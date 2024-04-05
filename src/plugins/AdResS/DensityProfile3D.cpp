//
// Created by alex on 29.12.23.
//

#include "DensityProfile3D.h"
#include <cstring>
#include <array>
#include <map>
#include <unordered_map>

void DensityProfile3D::init(double binWidth, Domain *domain, double rho0, double smoothingFactor) {
    _binWidth = binWidth;
    _rho0 = rho0;
    double tmpMult = domain->getGlobalLength(0) * domain->getGlobalLength(1) * domain->getGlobalLength(2);
    for (int d = 0; d < 3; ++d) {
        _binDims[d] = static_cast<unsigned long>(domain->getGlobalLength(d) / binWidth);
        _binVolumes[d] = tmpMult / domain->getGlobalLength(d) * binWidth;
    }

    Interpolation::createGaussianMatrix(0.0, domain->getGlobalLength(0), binWidth, smoothingFactor, _smoothingFilter);
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

    //divide global bins by volume
    for(int d = 0; d < 3; d++) {
        for(int i = 0; i < _globalDensities[d].size(); i++) {
            _globalDensities[d][i] /= _binVolumes[d];
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

void DensityProfile3D::resetBuffers() {
    for (int d = 0; d < 3; ++d) {
        _localDensities[d].clear();
        _localDensities[d].resize(_binDims[d], 0.0);

        _globalDensities[d].clear();
        _globalDensities[d].resize(_binDims[d], 0.0);
    }
}

static void map_add(std::map<double, double> &in, std::map<double, double> &out) {
    for(auto [pos, val] : in ) {
        out[pos] += val;
    }
}

static void umap_add(std::unordered_map<double, double> &in, std::unordered_map<double, double> &out) {
    for(auto [pos, val] : in ) {
        out[pos] += val;
    }
}

static std::unordered_map<double,double> create_map() {
    return {};
}

void DensityProfile3D::computeDensities(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
                                        Domain *domain) {
    std::array<std::unordered_map<double,double>,3> buffer;
    std::array<std::map<double,double>,3> global_forces;

    //Load all forces
    #if defined(_OPENMP)
    #pragma omp declare reduction(map_plus : std::unordered_map<double, double> : umap_add(omp_in, omp_out)) \
        initializer (omp_priv=create_map())
    auto* raw_array = buffer.data();
    #pragma omp parallel reduction(map_plus: raw_array[:3])
    #endif
    for(auto itM = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); itM.isValid(); ++itM) {
        std::array<double, 3> R = (*itM).r_arr();
        for(int d = 0; d < 3; d++) {
            raw_array[d][R[d]] += (*itM).F(d);
        }
    }

    //share maps across all ranks
#if defined(ENABLE_MPI)
    int num_ranks = 0;
    int rank = domainDecomp->getRank();
    MPI_Comm_size(domainDecomp->getCommunicator(), &num_ranks);

    std::vector<std::size_t> local_entry_counts;
    local_entry_counts.resize(num_ranks * 3, 0UL);
    for(int d = 0; d < 3; d++) local_entry_counts[3 * rank + d] = buffer[d].size();

    std::vector<std::size_t> global_entry_counts;
    global_entry_counts.resize(num_ranks * 3, 0UL);
    MPI_Allreduce(local_entry_counts.data(), global_entry_counts.data(), num_ranks * 3, MPI_UNSIGNED_LONG, MPI_SUM, domainDecomp->getCommunicator());

    std::array<std::size_t, 3> entries_per_dim {0,0,0};
    std::array<std::size_t, 3> prior_entries_per_dim {0,0,0};
    for(int r = 0; r < num_ranks; r++) {
        for(int d = 0; d < 3; d++) {
            entries_per_dim[d] += global_entry_counts[3 * r + d];
            if(r < rank) prior_entries_per_dim[d] += global_entry_counts[3 * r + d];
        }
    }

    // prepare data send buffers, we need twice as much space, since we are also sending the keys
    std::array<std::vector<double>, 3> local_entries;
    for(int d = 0; d < 3; d++) local_entries[d].resize(2 * entries_per_dim[d], 0.0);
    for(int d = 0; d < 3; d++) {
        unsigned long offset = 2 * prior_entries_per_dim[d];
        for(auto [pos, val] : buffer[d]) {
            local_entries[d][offset + 0] = pos;
            local_entries[d][offset + 1] = val;
            offset += 2;
        }
    }

    std::array<std::vector<double>, 3> global_entries;
    for(int d = 0; d < 3; d++) global_entries[d].resize(2 * entries_per_dim[d], 0.0);
    for(int d = 0; d < 3; d++) {
        MPI_Allreduce(local_entries[d].data(), global_entries[d].data(), 2 * entries_per_dim[d], MPI_DOUBLE, MPI_SUM, domainDecomp->getCommunicator());
    }

    for(int d = 0; d < 3; d++) {
        for(unsigned long i = 0; i < 2 * entries_per_dim[d]; i+=2) {
            global_forces[d][global_entries[d][i+0]] += global_entries[d][i+1];
        }
    }
#else
    for(int d = 0; d < 3; d++) {
        global_forces[d].insert(buffer[d].begin(), buffer[d].end());
    }
#endif

    //normalize data
    {
        unsigned long num_mols = _simulation.getTotalNumberOfMolecules();
        for(int d = 0; d < 3; d++) {
            for(auto it = global_forces[d].begin(); it != global_forces[d].end(); ++it) {
                global_forces[d][it->first] /= num_mols;
            }
        }
    }

    std::array<std::vector<double>,3> forces;
    std::array<std::vector<double>,3> steps;
    for(int d = 0; d < 3; d++) {
        forces[d].reserve(global_forces[d].size());
        steps[d].reserve(global_forces[d].size()-1);
        double last_pos = 0.0;
        for(auto [pos, val] : global_forces[d]) {
            steps[d].emplace_back(pos - last_pos);
            forces[d].emplace_back(val);

            last_pos = pos;
        }
    }

    std::array<Interpolation::Function, 3> force_functions;
    for(int d = 0; d < 3; d++) Interpolation::computeHermite(0.0, forces[d], steps[d], forces[d].size(), force_functions[d]);
    std::array<Interpolation::Function, 3> force_function_ints;
    for(int d = 0; d < 3; d++) Interpolation::computeIntegral(force_functions[d], force_function_ints[d]);

    double T = domain->getGlobalCurrentTemperature();
    const double kB = 1.380649e-23;

    double scaling = 1 / (T * kB);

    for(int d = 0; d < 3; d++) {
        for(unsigned long i = 0; i < force_function_ints[d].n; i++) {
            force_function_ints[d].function_values[i] *= scaling;
            force_function_ints[d].gradients[i] *= scaling;

            force_function_ints[d].function_values[i] += _rho0;
        }
    }

    _histDensities = std::move(force_function_ints);
}

const Interpolation::Function &DensityProfile3D::getHistDensity(int dim) const {
    return _histDensities[dim];
}

void DensityProfile3D::writeDensity(const std::string &filename, const std::string &separator, int dim, bool smoothed) {
	std::vector<double> densities;
	if(smoothed) densities = getDensitySmoothed(dim);
	else densities = getDensity(dim);

#ifdef ENABLE_MPI
	if (_simulation.domainDecomposition().getRank() != 0) return;
#endif
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
