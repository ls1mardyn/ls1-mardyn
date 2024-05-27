//
// Created by alex on 29.12.23.
//

#include "DensityProfile3D.h"
#include <cstring>
#include <array>
#include <vector>

void DensityProfile3D::init(double binWidth, Domain *domain, double smoothingFactor, double gridRadius) {
    _binWidth = binWidth;
	_smoothingFactor = smoothingFactor;
    double tmpMult = domain->getGlobalLength(0) * domain->getGlobalLength(1) * domain->getGlobalLength(2);
    for (int d = 0; d < 3; ++d) {
        _binDims[d] = static_cast<unsigned long>(domain->getGlobalLength(d) / binWidth);
        _binVolumes[d] = tmpMult / domain->getGlobalLength(d) * binWidth;
		_gmmDensities[d].n = 0;
    }

    Interpolation::createGaussianMatrix(0.0, domain->getGlobalLength(0), binWidth, smoothingFactor, _smoothingFilter);

	_gridSampler.SetMeasureRadius(gridRadius);
	_grid.MeshAllDomain();
	_grid.StartGrid(_binDims[0], _binDims[1], _binDims[2]);

	_gridSampler.GetGridHandler().SetGridBoundarySubsets(&_grid);
	_gridSampler.init(nullptr);
	_averager.setDataSize(_gridSampler.GetSamples().material_density);

	//ofstream mesh("mesh.txt");
	//mesh << grid;
	//mesh.close();
}

void DensityProfile3D::resetBuffers() {
	for (int d = 0; d < 3; ++d) {
		_localDensities[d].clear();
		_localDensities[d].resize(_binDims[d], 0.0);

		_globalDensities[d].clear();
		_globalDensities[d].resize(_binDims[d], 0.0);
	}
}

void DensityProfile3D::sampleDensities(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
                                       Domain *domain) {
    resetBuffers();
    //we first count the molecule counts
    //only after the global reduction, we divide by bin volume
    //gather local densities first
    #if defined(_OPENMP)
    #pragma omp declare reduction(vec_plus : std::vector<double> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<>())) \
        initializer(omp_priv = decltype(omp_orig)(omp_orig.size(), 0.0))
    auto* raw_array = _localDensities.data();
    #pragma omp parallel reduction(vec_plus: raw_array[:3])
    #endif
    for(auto itM = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); itM.isValid(); ++itM) {
        std::array<double, 3> R = (*itM).r_arr();
        for(int d = 0; d < 3; d++) {
            auto bin_pos = static_cast<unsigned long>(R[d] / _binWidth);
            bin_pos = std::clamp(bin_pos, 0UL, _binDims[d]-1UL);
            raw_array[d][bin_pos] += 1;
        }
    }

#if defined(ENABLE_MPI)
    MPI_Allreduce(_localDensities[0].data(), _globalDensities[0].data(), _localDensities[0].size(), MPI_DOUBLE, MPI_SUM, domainDecomp->getCommunicator());
    MPI_Allreduce(_localDensities[1].data(), _globalDensities[1].data(), _localDensities[1].size(), MPI_DOUBLE, MPI_SUM, domainDecomp->getCommunicator());
    MPI_Allreduce(_localDensities[2].data(), _globalDensities[2].data(), _localDensities[2].size(), MPI_DOUBLE, MPI_SUM, domainDecomp->getCommunicator());
#else
    std::memcpy(_localDensities[0].data(), _globalDensities[0].data(), _localDensities[0].size());
    std::memcpy(_localDensities[1].data(), _globalDensities[1].data(), _localDensities[1].size());
    std::memcpy(_localDensities[2].data(), _globalDensities[2].data(), _localDensities[2].size());
#endif

    /*//acquire global information
    domainDecomp->collCommInit(_localDensities[0].size() + _localDensities[1].size() + _localDensities[2].size());
    for(int d = 0; d < 3; d++) {
        for(int i : _localDensities[d]) {
            domainDecomp->collCommAppendInt(i);
        }
    }

    domainDecomp->collCommAllreduceSum();
    for(int d = 0; d < 3; d++) {
        for(int i = 0; i < _globalDensities[d].size(); i++) {
            _globalDensities[d][i] = domainDecomp->collCommGetInt();
        }
    }
    domainDecomp->collCommFinalize();*/
	//TODO: why is densities multidimensional? One per bin aint it?
    //divide global bins by volume
    for(int d = 0; d < 3; d++) {
        for(int i = 0; i < _globalDensities[d].size(); i++) {
            _globalDensities[d][i] /= _binVolumes[d];
        }
    }
}

void DensityProfile3D::computeGMMDensities(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
    std::array<std::vector<double>,3> global_mol_pos = getGlobalMolPos(particleContainer, domainDecomp, domain);
	for(int d = 0; d < 3; d++) {
		Interpolation::createGMM(0.0, domain->getGlobalLength(d), _binDims[d], _smoothingFactor, global_mol_pos[d], _gmmDensities[d]);
	}
}

void DensityProfile3D::computeFTDensities(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
										  Domain *domain) {
	std::array<std::vector<double>,3> global_mol_pos = getGlobalMolPos(particleContainer, domainDecomp, domain);

	for(int d = 0; d < 3; d++) {
		const auto dom_size = domain->getGlobalLength(d);
		std::vector<std::complex<double>> F;
		Interpolation::realFT(global_mol_pos[d], 1000, dom_size, F);
		Interpolation::filterFT(F);
		Interpolation::ift(F, 1000, dom_size, 0, dom_size, _binDims[d], _ftDensities[d]);

		for(unsigned long i = 0; i < _ftDensities[d].n; i++) {
			_ftDensities[d].function_values[i] /= _binVolumes[d];
		}
	}
}

const std::vector<double> &DensityProfile3D::getDensity(int dim) const {
	mardyn_assert(((dim>=0) && (dim<=3)));
	return _globalDensities[dim];
}

std::vector<double> DensityProfile3D::getDensitySmoothed(int dim) const {
	mardyn_assert(((dim>=0) && (dim<=3)));
	mardyn_assert((dim==0)); // TODO make for all dims a filter
	return _smoothingFilter * _globalDensities[dim];
}

const Interpolation::Function &DensityProfile3D::getGMMDensity(int dim) const {
    return _gmmDensities[dim];
}

const Interpolation::Function &DensityProfile3D::getFTDensity(int dim) const {
	return _ftDensities[dim];
}

void DensityProfile3D::writeDensity(const std::string &filename, const std::string &separator, int dim, Type type) {
#ifdef ENABLE_MPI
	if (_simulation.domainDecomposition().getRank() != 0) return;
#endif
	if(type == GMM) {
		getGMMDensity(dim).writeXML(filename);
		return;
	}
	if(type == FT) {
		getFTDensity(dim).writeXML(filename);
		return;
	}
	if(type == GRID){
		std::ofstream ppc2(filename);
		_gridSampler.WritePlaneSamples<std::vector<double>>(ppc2, _averager.getAveragedDataCopy());
		ppc2.close();
		return;
	}

	std::vector<double> densities;
	if(type == SMOOTH) densities = getDensitySmoothed(dim);
	else if(type == SAMPLE) densities = getDensity(dim);
	try {
		std::ofstream file {filename};
		for(unsigned long i = 0; i < densities.size()-1; i++) {
			file << densities[i] << separator;
		}
		file << densities[densities.size()-1] << std::endl;
	} catch (std::ifstream::failure& e) {
		Log::global_log->error() << "[AdResS] Failed to write densities to file.\n" << e.what() << std::endl;
		_simulation.exit(-1);
	}
}

void DensityProfile3D::step(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
	std::array<std::vector<double>,3> global_mol_pos = getGlobalMolPos(particleContainer, domainDecomp, domain);
	_gridSampler.SampleAtNodes(global_mol_pos);
	_gridSampler.SetSubsetMaterialDensityValues();
	_averager.averageData(_gridSampler.GetSamples().material_density);
}

std::array<std::vector<double>,3> DensityProfile3D::getGlobalMolPos(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
	std::array<std::vector<double>,3> local_mol_pos;
	std::array<std::vector<double>,3> global_mol_pos;
	local_mol_pos[0].resize(particleContainer->getNumberOfParticles(ParticleIterator::ONLY_INNER_AND_BOUNDARY), 0.0);
	local_mol_pos[1].resize(particleContainer->getNumberOfParticles(ParticleIterator::ONLY_INNER_AND_BOUNDARY), 0.0);
	local_mol_pos[2].resize(particleContainer->getNumberOfParticles(ParticleIterator::ONLY_INNER_AND_BOUNDARY), 0.0);
	global_mol_pos[0].resize(domain->getglobalNumMolecules(true, particleContainer, domainDecomp), 0.0);
	global_mol_pos[1].resize(domain->getglobalNumMolecules(true, particleContainer, domainDecomp), 0.0);
	global_mol_pos[2].resize(domain->getglobalNumMolecules(true, particleContainer, domainDecomp), 0.0);

	//Load all local positions
	auto begin = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
#if defined(_OPENMP)
#pragma omp parallel
#endif
	for(auto itM = begin; itM.isValid(); ++itM) {
		const std::array<double, 3> R = (*itM).r_arr();
		local_mol_pos[0][itM->getID()] = R[0];
		local_mol_pos[1][itM->getID()] = R[1];
		local_mol_pos[2][itM->getID()] = R[2];
	}

	//share maps across all ranks
#if defined(ENABLE_MPI)
	int num_ranks = 0;
	int rank = domainDecomp->getRank();
	MPI_Comm_size(domainDecomp->getCommunicator(), &num_ranks);

	std::vector<std::size_t> rank_mol_counts;
	rank_mol_counts.resize(num_ranks, 0UL);
	rank_mol_counts[rank] = local_mol_pos[0].size();

	MPI_Allreduce(MPI_IN_PLACE, rank_mol_counts.data(), num_ranks, MPI_UNSIGNED_LONG, MPI_SUM, domainDecomp->getCommunicator());

	std::size_t prior_mol_count = 0;
	for(int r = 0; r < rank; r++) {
		prior_mol_count += rank_mol_counts[r];
	}

	// write data to correct positions
	std::size_t offset = prior_mol_count;
	std::memcpy(global_mol_pos[0].data() + offset, local_mol_pos[0].data(), rank_mol_counts[rank]*sizeof(double));
	std::memcpy(global_mol_pos[1].data() + offset, local_mol_pos[1].data(), rank_mol_counts[rank]*sizeof(double));
	std::memcpy(global_mol_pos[2].data() + offset, local_mol_pos[2].data(), rank_mol_counts[rank]*sizeof(double));

	// share across ranks
	MPI_Allreduce(MPI_IN_PLACE, global_mol_pos[0].data(), global_mol_pos[0].size(), MPI_DOUBLE, MPI_SUM, domainDecomp->getCommunicator());
	MPI_Allreduce(MPI_IN_PLACE, global_mol_pos[1].data(), global_mol_pos[1].size(), MPI_DOUBLE, MPI_SUM, domainDecomp->getCommunicator());
	MPI_Allreduce(MPI_IN_PLACE, global_mol_pos[2].data(), global_mol_pos[2].size(), MPI_DOUBLE, MPI_SUM, domainDecomp->getCommunicator());
#else
	global_mol_pos[0] = std::move(local_mol_pos[0]);
    global_mol_pos[1] = std::move(local_mol_pos[1]);
    global_mol_pos[2] = std::move(local_mol_pos[2]);
#endif

	return global_mol_pos;
}
