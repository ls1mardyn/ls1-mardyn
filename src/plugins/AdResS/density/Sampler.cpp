/*
 * Created on Wed May 08 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */


#include "Sampler.h"

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"

/*************************************************************
 *						GRID_SAMPLER
 *********************************************************** */
GridSampler::GridSampler(FTH::grid_t *grid, double rad): _grid(grid), _measure_radius(rad) { }

void GridSampler::init(Domain *domain) {
	unsigned long total_nodes = _grid->getNodes().size();
	_sampled_data.resize(total_nodes, 0.0);
}

void GridSampler::sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) {
	sampleAtNodes(pc);
	_was_sampled = true;
}

void GridSampler::sampleAtNodes(ParticleContainer* pc){
	std::fill(_sampled_data.begin(), _sampled_data.end(), 0);
	auto& nodes = _grid->getNodes();
    auto* domain = _simulation.getDomain();
    auto& decomp = _simulation.domainDecomposition();
    const double r_cutoff = _simulation.getcutoffRadius();
    const std::array<double, 3> lower {
        decomp.getBoundingBoxMin(0, domain) - r_cutoff,
        decomp.getBoundingBoxMin(1, domain) - r_cutoff,
        decomp.getBoundingBoxMin(2, domain) - r_cutoff
    };

	// put mol pos into grid of spacing 2x measure_radius
	const auto box_dim = 2 * _measure_radius;
	const std::array<int, 3> bins_per_dim = { static_cast<int>(std::ceil((decomp.getBoundingBoxMax(0, domain) + r_cutoff - lower[0]) / box_dim)),
											  static_cast<int>(std::ceil((decomp.getBoundingBoxMax(1, domain) + r_cutoff - lower[1]) / box_dim)),
											  static_cast<int>(std::ceil((decomp.getBoundingBoxMax(2, domain) + r_cutoff - lower[2]) / box_dim))};
	const auto bins = bins_per_dim[0]*bins_per_dim[1]*bins_per_dim[2];
	std::vector<std::vector<std::array<double, 3>>> mol_pos_binned;
	mol_pos_binned.resize(bins);
	for(auto it = pc->iterator(ParticleIterator::ALL_CELLS); it.isValid(); ++it) {
		std::array<double, 3> part_pos {it->r_arr()};
		int bin_x = std::min(static_cast<int>((part_pos[0] - lower[0]) / box_dim), bins_per_dim[0]-1);
		int bin_y = std::min(static_cast<int>((part_pos[1] - lower[1]) / box_dim), bins_per_dim[1]-1);
		int bin_z = std::min(static_cast<int>((part_pos[2] - lower[2]) / box_dim), bins_per_dim[2]-1);

		mol_pos_binned[bins_per_dim[0]*bins_per_dim[1]*bin_z + bins_per_dim[0]*bin_y + bin_x].push_back(part_pos);
	}

	#if defined(_OPENMP)
	#pragma omp parallel for
	#endif
	for(int nidx = 0; nidx < _sampled_data.size(); nidx++) {
		auto nodal_pos=nodes[nidx].getPos();
		int bin_x = static_cast<int>((nodal_pos[0] - lower[0]) / box_dim);
		int bin_y = static_cast<int>((nodal_pos[1] - lower[1]) / box_dim);
		int bin_z = static_cast<int>((nodal_pos[2] - lower[2]) / box_dim);

		for(int offset_z = -1; offset_z <= 1; offset_z++) {
			int bz = offset_z + bin_z;
			if(bz < 0 || bz >= bins_per_dim[2]) continue;
			for(int offset_y = -1; offset_y <= 1; offset_y++) {
				int by = offset_y + bin_y;
				if(by < 0 || by >= bins_per_dim[1]) continue;
				for(int offset_x = -1; offset_x <= 1; offset_x++) {
					int bx = offset_x + bin_x;
					if(bx < 0 || bx >= bins_per_dim[0]) continue;

					for(auto& mol_pos : mol_pos_binned[bins_per_dim[0]*bins_per_dim[1]*bz + bins_per_dim[0]*by + bx]) {
						if(ParticleInsideMeasuringSpace(nodal_pos, mol_pos)){
							_sampled_data[nidx] +=1;
						}
					}
				}
			}
		}
	}
	//Convert to material density values
	const double sphere_volume = 4.0/3.0 * M_PI * std::pow(_measure_radius, 3);
	#if defined(_OPENMP)
	#pragma omp parallel for
	#endif
	for(int nidx = 0; nidx < _sampled_data.size(); nidx++){
		_grid->getNodes()[nidx].data().particles = static_cast<int>(_sampled_data[nidx]);
		const auto mat_density = (double)_sampled_data[nidx]/sphere_volume;
		_sampled_data[nidx] = mat_density;
		_grid->getNodes()[nidx].data().density = mat_density;
	}
}

bool GridSampler::ParticleInsideMeasuringSpace(const std::array<double,3>& nodal_pos, const std::array<double,3>& par_pos) const {
	bool is_inside = false;

	std::array<double, 3> distance { };
	distance[0] = par_pos[0]-nodal_pos[0];
	distance[1] = par_pos[1]-nodal_pos[1];
	distance[2] = par_pos[2]-nodal_pos[2];

	double norm2 = distance[0]*distance[0] + distance[1]*distance[1] + distance[2]*distance[2];

	if(norm2 < _measure_radius*_measure_radius){
		is_inside = true;
	}

	return is_inside;
}

void GridSampler::writeSample(const std::string &filename, int simstep) {
	std::stringstream ss;
	ss << filename << "_" << simstep << ".txt";
    std::ofstream file(ss.str());
    file << "#Total nodes sampled: " << _sampled_data.size() << "\n";
    auto& nodes = _grid->getNodes();

    for(int i = 0; i < _sampled_data.size(); i++) {
        auto& pos = nodes[i].getPos();
        file << i << "\t" << pos[0] << "\t" << pos[1] << "\t" << pos[2] << "\t" << _sampled_data[i] << "\n";
    }
}

/*************************************************************
 *						AVG_GRID_SAMPLER
 *********************************************************** */
AveragedGridSampler::AveragedGridSampler(FTH::grid_t *grid, double rad): GridSampler(grid, rad) { }

void AveragedGridSampler::init(Domain *domain) {
	GridSampler::init(nullptr);
	_averager.setDataSize(_sampled_data);
}

void AveragedGridSampler::sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) {
	GridSampler::sampleData(pc, domainDecomp, domain);
	_averager.averageData(_sampled_data);
	_averager.getAveragedData(_sampled_data);

	//write averages into grid
	#if defined(_OPENMP)
	#pragma omp parallel for
	#endif
	for(int nidx = 0; nidx < _sampled_data.size(); nidx++){
		_grid->getNodes()[nidx].data().density = _sampled_data[nidx];
	}
}

void AveragedGridSampler::writeSample(const std::string &filename, int simstep) {
	GridSampler::writeSample(filename, simstep);

	std::stringstream ss;
	ss << "AVG_" << filename << "_" << simstep << ".txt";
	std::ofstream avg_file(ss.str());
	_averager.writeAverage(avg_file, _averager.getSumData());
	avg_file.close();
}

Averager<std::vector<double>> &AveragedGridSampler::getAverager() {
	return _averager;
}

/*************************************************************
 *						PROJ_SAMPLER
 *********************************************************** */
void ProjectedSampler::loadGlobalMolPos(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
										Domain *domain) {
	_local_mol_pos.resize(particleContainer->getNumberOfParticles(ParticleIterator::ONLY_INNER_AND_BOUNDARY), 0.0);
	_global_mol_pos.resize(domain->getglobalNumMolecules(true, particleContainer, domainDecomp), 0.0);
	std::fill(_global_mol_pos.begin(), _global_mol_pos.end(), 0.0);

	//Load all local positions
	auto begin = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	for(struct {int idx; ParticleIterator it; } s = {0, begin}; s.it.isValid(); ++s.it, ++s.idx) {
		const std::array<double, 3> R = (*s.it).r_arr();
		_local_mol_pos[s.idx] = R[_dim];
	}

	//share maps across all ranks
#if defined(ENABLE_MPI)
	int num_ranks = 0;
	int rank = domainDecomp->getRank();
	MPI_Comm_size(domainDecomp->getCommunicator(), &num_ranks);

	std::vector<std::size_t> rank_mol_counts;
	rank_mol_counts.resize(num_ranks, 0UL);
	rank_mol_counts[rank] = _local_mol_pos.size();

	MPI_Allreduce(MPI_IN_PLACE, rank_mol_counts.data(), num_ranks, MPI_UNSIGNED_LONG, MPI_SUM, domainDecomp->getCommunicator());

	std::size_t prior_mol_count = 0;
	for(int r = 0; r < rank; r++) {
		prior_mol_count += rank_mol_counts[r];
	}

	// write data to correct positions
	std::size_t offset = prior_mol_count;
	std::memcpy(_global_mol_pos.data() + offset, _local_mol_pos.data(), rank_mol_counts[rank]*sizeof(double));

	// share across ranks
	MPI_Allreduce(MPI_IN_PLACE, _global_mol_pos.data(), _global_mol_pos.size(), MPI_DOUBLE, MPI_SUM, domainDecomp->getCommunicator());
#else
	std::copy(_local_mol_pos.begin(), _local_mol_pos.end(), _global_mol_pos.begin());
#endif
}

/*************************************************************
 *						DIR_PROJ_SAMPLER
 *********************************************************** */
DirectProjectedSampler::DirectProjectedSampler(int dim, double bin_width) : ProjectedSampler(dim), _bins(0), _binWidth(bin_width), _binVolume(0) { }

void DirectProjectedSampler::init(Domain *domain) {
	double full_volume = domain->getGlobalLength(0) * domain->getGlobalLength(1) * domain->getGlobalLength(2);
	_bins = static_cast<unsigned long>(domain->getGlobalLength(_dim) / _binWidth);
	_binVolume = full_volume / domain->getGlobalLength(_dim) * _binWidth;
}

void DirectProjectedSampler::sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) {
	_localDensities.clear();
	_localDensities.resize(_bins, 0.0);

	_sampled_data.clear();
	_sampled_data.resize(_bins, 0.0);

	//we first count the molecule counts
	//only after the global reduction, we divide by bin volume
	//gather local densities first
	#if defined(_OPENMP)
	auto* raw_array = _localDensities.data();
	#pragma omp parallel reduction(+: raw_array[:_localDensities.size()])
	#endif
	for(auto itM = pc->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); itM.isValid(); ++itM) {
		std::array<double, 3> R = (*itM).r_arr();
		auto bin_pos = static_cast<unsigned long>(R[_dim] / _binWidth);
		bin_pos = std::clamp(bin_pos, 0UL, _bins-1UL);
		raw_array[bin_pos] += 1;
	}

#if defined(ENABLE_MPI)
	MPI_Allreduce(_localDensities.data(), _sampled_data.data(), _localDensities.size(), MPI_DOUBLE, MPI_SUM, domainDecomp->getCommunicator());
#else
	std::memcpy(_localDensities.data(), _sampled_data.data(), _localDensities.size());
#endif

	//divide global bins by volume
	for(int i = 0; i < _sampled_data.size(); i++) {
		_sampled_data[i] /= _binVolume;
	}

	_was_sampled = true;
}

void DirectProjectedSampler::writeSample(const std::string &filename, int simstep) {
	std::stringstream ss;
	ss << filename << "_" << simstep << ".txt";
	std::ofstream file(ss.str());
	file << "#Total bins sampled: " << _sampled_data.size() << "\n";

	for(int i = 0; i < _sampled_data.size(); i++) {
		file << i << "\t" << _sampled_data[i] << "\n";
	}
}

/*************************************************************
 *						SMOOTH_PROJ_SAMPLER
 *********************************************************** */
SmoothedProjectedSampler::SmoothedProjectedSampler(int dim, double bin_width, double smoothing_strength) : DirectProjectedSampler(dim, bin_width), _filterStrength(smoothing_strength) { }

void SmoothedProjectedSampler::init(Domain *domain) {
	DirectProjectedSampler::init(domain);
	Interpolation::createGaussianMatrix(0.0, domain->getGlobalLength(0), _binWidth, _filterStrength, _smoothingFilter);
}

void SmoothedProjectedSampler::sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) {
	DirectProjectedSampler::sampleData(pc, domainDecomp, domain);
	std::vector<double> buffer = _smoothingFilter * _sampled_data;
	std::copy(buffer.begin(), buffer.end(), _sampled_data.begin());
}

/*************************************************************
 *						FT_PROJ_SAMPLER
 *********************************************************** */
FTProjectedSampler::FTProjectedSampler(int dim, int frequencies, int samples) : ProjectedSampler(dim), _frequencies(frequencies), _samples(samples), _realFT(), _fun(), _binVolume(0) { }

void FTProjectedSampler::init(Domain *domain) {
	double binWidth = domain->getGlobalLength(_dim) / _samples;
	_binVolume = domain->getGlobalLength(0) * domain->getGlobalLength(1) * domain->getGlobalLength(2);
	_binVolume = _binVolume / domain->getGlobalLength(_dim) * binWidth;
}

void FTProjectedSampler::sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) {
	ProjectedSampler::loadGlobalMolPos(pc, domainDecomp, domain);
	const auto dom_size = domain->getGlobalLength(_dim);
	Interpolation::realFT(_global_mol_pos, _frequencies, dom_size, _realFT);
	Interpolation::filterFT(_realFT);
	Interpolation::ift(_realFT, _frequencies, dom_size, 0, dom_size, _samples, _fun);

	for(unsigned long i = 0; i < _fun.n; i++) {
		_fun.function_values[i] /= _binVolume;
	}

	_was_sampled = true;
}

void FTProjectedSampler::writeSample(const std::string &filename, int simstep) {
	std::stringstream ss;
	ss << filename << "_" << simstep << ".txt";
	_fun.writeTXT(ss.str());
}

void FTProjectedSampler::getSampleFunction(Interpolation::Function &fun) {
	fun = _fun;
}

/*************************************************************
 *						GMM_PROJ_SAMPLER
 *********************************************************** */
GMMProjectedSampler::GMMProjectedSampler(int dim, int samples, double smoothing_strength) : ProjectedSampler(dim), _samples(samples), _filterStrength(smoothing_strength) {}

void GMMProjectedSampler::sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) {
	ProjectedSampler::loadGlobalMolPos(pc, domainDecomp, domain);
	Interpolation::createGMM(0.0, domain->getGlobalLength(_dim), _samples, _filterStrength, _global_mol_pos, _fun);

	_was_sampled = true;
}

void GMMProjectedSampler::writeSample(const std::string &filename, int simstep) {
	std::stringstream ss;
	ss << filename << "_" << simstep << ".txt";
	_fun.writeTXT(ss.str());
}

void GMMProjectedSampler::getSampleFunction(Interpolation::Function &fun) {
	fun = _fun;
}