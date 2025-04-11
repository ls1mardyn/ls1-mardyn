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

void SamplerBase::checkReset(unsigned long simstep) {
    if (shouldReset(simstep)) _averager.reset();
}

void SamplerBase::writeAverage(const std::string &filename, int simstep) {
    auto avg = getAverage();
    const auto& sum_counts = _averager.getSumData();

    std::stringstream ss;
    ss << filename << "_" << simstep << ".txt";
    std::ofstream avg_file(ss.str());

    std::string prefix ="//[TimeAverage]: ";
    avg_file << prefix << "data average after: " << _averager.getStepCount() << " steps" << "\n";
    avg_file << prefix << "data structure with size: " << avg.size() << "\n";

    for (int i = 0; i < avg.size(); i++) {
        avg_file << i << "\t" << sum_counts[i] << "\t" << avg[i] << "\n";
    }

    avg_file.close();
}

bool SamplerBase::shouldReset(unsigned long simstep) {
    if (!_useReset) return false;
    if (_averager.getStepCount() > 0 && simstep > 0 && simstep % _window == 0) return true;
    return false;
}

/*************************************************************
 *						GRID_SAMPLER
 *********************************************************** */
GridSampler::GridSampler(int window, bool useReset, FTH::grid_t *grid, double rad): SamplerBase(window, useReset), _grid(grid), _measure_radius(rad) { }

void GridSampler::init(Domain *domain) {
	unsigned long total_nodes = _grid->getNodes().size();
	_sampled_data.resize(total_nodes, 0.0);
    _averager.setDataSize(_sampled_data);
}

void GridSampler::sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) {
	sampleAtNodes(pc);
	_was_sampled = true;

    _averager.averageData(_sampled_data); // _sample_data buffer is now free to use
    _averager.getAveragedData(_sampled_data);

    //write averages into grid
    const double sphere_volume = 4.0/3.0 * M_PI * std::pow(_measure_radius, 3);
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for(int nidx = 0; nidx < _sampled_data.size(); nidx++){
        _grid->getNodes()[nidx].data().density = _sampled_data[nidx] / sphere_volume;
    }
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

std::vector<double> GridSampler::getAverage() {
    auto avg = _averager.getAveragedDataCopy();

    const double sphere_volume = 4.0/3.0 * M_PI * std::pow(_measure_radius, 3);
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for (int idx = 0; idx < avg.size(); idx++) {
        avg[idx] /= sphere_volume;
    }
    return avg;
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
DirectProjectedSampler::DirectProjectedSampler(int window, bool useReset, int dim, int bins) : ProjectedSampler(window, useReset, dim), _bins(bins), _binWidth(0), _binVolume(0) { }

void DirectProjectedSampler::init(Domain *domain) {
	double full_volume = domain->getGlobalLength(0) * domain->getGlobalLength(1) * domain->getGlobalLength(2);
    _binWidth = domain->getGlobalLength(_dim) / _bins;
	_binVolume = full_volume / domain->getGlobalLength(_dim) * _binWidth;

    _sampled_data.resize(_bins, 0.0);
}

void DirectProjectedSampler::sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) {
    ProjectedSampler::loadGlobalMolPos(pc, domainDecomp, domain);
    std::fill(_sampled_data.begin(), _sampled_data.end(), 0.0);

	#if defined(_OPENMP)
    #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
			std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
			initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
	#pragma omp parallel for reduction(vec_double_plus: _sampled_data)
	#endif
	for(int idx = 0; idx < _global_mol_pos.size(); idx++) {
		auto bin_pos = static_cast<unsigned long>(_global_mol_pos[idx] / _binWidth);
		bin_pos = std::clamp(bin_pos, 0UL, _bins-1UL);
		_sampled_data[bin_pos] += 1.0;
	}

    _averager.averageData(_sampled_data);
	_was_sampled = true;
}

std::vector<double> DirectProjectedSampler::getAverage() {
    auto avg = _averager.getAveragedDataCopy();

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for(int i = 0; i < avg.size(); i++) {
        avg[i] /= _binVolume;
    }
    return avg;
}

/*************************************************************
 *						SMOOTH_PROJ_SAMPLER
 *********************************************************** */
SmoothedProjectedSampler::SmoothedProjectedSampler(int window, bool useReset, int dim, int bins, double smoothing_strength) : DirectProjectedSampler(window, useReset, dim, bins), _filterStrength(smoothing_strength) { }

void SmoothedProjectedSampler::init(Domain *domain) {
	DirectProjectedSampler::init(domain);
	Interpolation::createGaussianMatrix(0.0, domain->getGlobalLength(0), _binWidth, _filterStrength, _smoothingFilter);
}

std::vector<double> SmoothedProjectedSampler::getAverage() {
    return _smoothingFilter * DirectProjectedSampler::getAverage();
}

/*************************************************************
 *						FT_PROJ_SAMPLER
 *********************************************************** */
FTProjectedSampler::FTProjectedSampler(int window, bool useReset, int dim, int bins, int frequencies) : DirectProjectedSampler(window, useReset, dim, bins), _frequencies(frequencies), _realFT(), _fun(), _binVolume(0) { }

void FTProjectedSampler::init(Domain *domain) {
    DirectProjectedSampler::init(domain);
}

std::vector<double> FTProjectedSampler::getAverage() {
    auto avg = DirectProjectedSampler::getAverage();
    auto freq = Interpolation::realFT(avg, _frequencies);
    return Interpolation::irft(freq, static_cast<int>(_bins));
}

/*************************************************************
 *						GMM_PROJ_SAMPLER
 *********************************************************** */
GMMProjectedSampler::GMMProjectedSampler(int window, bool useReset, int dim, int bins, double smoothing_strength) : ProjectedSampler(window, useReset, dim), _bins(bins), _filterStrength(smoothing_strength) {}

void GMMProjectedSampler::init(Domain *domain) {
    _sampled_data.resize(_bins, 0.0);
}

void GMMProjectedSampler::sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) {
	ProjectedSampler::loadGlobalMolPos(pc, domainDecomp, domain);
	Interpolation::createGMM(0.0, domain->getGlobalLength(_dim), _bins, _filterStrength, _global_mol_pos, _fun);

    std::copy(_fun.function_values.begin(), _fun.function_values.end(), _sampled_data.begin());
    _averager.averageData(_sampled_data);
	_was_sampled = true;
}

std::vector<double> GMMProjectedSampler::getAverage() {
    return _averager.getAveragedDataCopy();
}
