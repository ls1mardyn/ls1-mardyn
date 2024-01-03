//
// Created by alex on 29.12.23.
//

#include "DensityProfile3D.h"
#include <cstring>

void DensityProfile3D::init(double binWidth, Domain *domain) {
    _binWidth = binWidth;
    double tmpMult = domain->getGlobalLength(0) * domain->getGlobalLength(1) * domain->getGlobalLength(2);
    for (int d = 0; d < 3; ++d) {
        _binDims[d] = static_cast<unsigned long>(domain->getGlobalLength(d) / binWidth);
        _binVolumes[d] = tmpMult / domain->getGlobalLength(d) * binWidth;
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

void DensityProfile3D::resetBuffers() {
    for (int d = 0; d < 3; ++d) {
        _localDensities[d].clear();
        _localDensities[d].resize(_binDims[d], 0.0);

        _globalDensities[d].clear();
        _globalDensities[d].resize(_binDims[d], 0.0);
    }
}
