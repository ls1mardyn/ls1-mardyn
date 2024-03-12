//
// Created by alex on 25.02.24.
//

#include "PressureProfiler.h"
#include "Simulation.h"
#include "ensemble/EnsembleBase.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"

#include <string>
/********************************************
 * Plugin Code
 * ******************************************/
void PressureProfiler::init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {
    _profiler->init(particleContainer, domainDecomp, domain);
}

void PressureProfiler::readXML(XMLfileUnits &xmlconfig) {
    std::array<std::string, ModeCount> mode_strings{"MOP", "Direct"};
    std::string mode_str = mode_strings[1];
    xmlconfig.getNodeValue("mode", mode_str);
    unsigned long index = std::distance(mode_strings.begin(),
                                        std::find(mode_strings.begin(), mode_strings.end(), mode_str));
    if (index >= ModeCount) {
        Log::global_log->warning() << "[Pressure Profiler] Invalid mode specified. Falling back to: " << mode_strings[1]
                              << std::endl;
        index = 1;
    }
    _mode = static_cast<PressureMode>(index);

    if(_mode == Direct) {
        _profiler = std::make_unique<DirectProfiler>();
    }
    else if (_mode == MOP) {
        _profiler = std::make_unique<MOPProfiler>();
    }
    else {
        Log::global_log->error() << "[Pressure Profiler] Unknown mode." << std::endl;
        _simulation.exit(-1);
    }

    _profiler->_bins_per_dim = {};
    xmlconfig.getNodeValue("bins/x", _profiler->_bins_per_dim[0]);
    xmlconfig.getNodeValue("bins/y", _profiler->_bins_per_dim[1]);
    xmlconfig.getNodeValue("bins/z", _profiler->_bins_per_dim[2]);

    std::string outdims = "xyz";
    xmlconfig.getNodeValue("outdims", outdims);
    if (outdims.find('x') == std::string::npos) _profiler->_reduce_dim[0] = true;
    if (outdims.find('y') == std::string::npos) _profiler->_reduce_dim[1] = true;
    if (outdims.find('z') == std::string::npos) _profiler->_reduce_dim[2] = true;



    _profiler->_max_samples = 100;
    xmlconfig.getNodeValue("samples", _profiler->_max_samples);
    _profiler->_sample_gap = 1;
    xmlconfig.getNodeValue("sampleFreq", _profiler->_sample_gap);
    _profiler->_write_gap = 100;
    xmlconfig.getNodeValue("writeFreq", _profiler->_write_gap);
    _profiler->_output_prefix = "pressures";
    xmlconfig.getNodeValue("filename", _profiler->_output_prefix);

	_start = 0;
	xmlconfig.getNodeValue("start", _start);
}

void PressureProfiler::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                               unsigned long simstep) {
	if(simstep < _start) return;
    _profiler->step(particleContainer, domainDecomp, domain, simstep);
}

void PressureProfiler::finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {

}

std::string PressureProfiler::getPluginName() {
    return {"PressureProfiler"};
}

/********************************************
 * Base Profiler code
 * ******************************************/
void ProfilerBase::writePressure(const std::string &filename, Vec3D<double> &pressure_data) {
#ifdef ENABLE_MPI
    if (_simulation.domainDecomposition().getRank() != 0) return;
#endif

    std::vector<double> output((_reduce_dim[0] ? 1 : _bins_per_dim[0]) *
                               (_reduce_dim[1] ? 1 : _bins_per_dim[1]) *
                               (_reduce_dim[2] ? 1 : _bins_per_dim[2]), 0
    );
    std::size_t start_x = 0;
    std::size_t start_y = 0;
    std::size_t start_z = 0;
    for(std::size_t index = 0; index < output.size(); index++) {
        std::size_t count = 0;
        double pressure = 0;
        for (std::size_t pos_z = start_z; pos_z < _bins_per_dim[2]; pos_z += _reduce_dim[2] ? 1 : _bins_per_dim[2]) {
            for (std::size_t pos_y = start_y; pos_y < _bins_per_dim[1]; pos_y += _reduce_dim[1] ? 1 : _bins_per_dim[1]) {
                for (std::size_t pos_x = start_x; pos_x < _bins_per_dim[0]; pos_x += _reduce_dim[0] ? 1 : _bins_per_dim[0]) {
                    count++;
                    pressure += pressure_data.at(pos_x, pos_y, pos_z);
                }
            }
        }
        pressure /= static_cast<double>(count);
        output[index] = pressure;

        bool carry = false;
        if(!_reduce_dim[0]) {
            if(start_x == _bins_per_dim[0]-1) carry = true;
            start_x = (start_x+1) % _bins_per_dim[0];
        }
        else carry = true;
        if(!_reduce_dim[1] && carry) {
            carry = false;
            if(start_y == _bins_per_dim[1]-1) carry = true;
            start_y = (start_y+1) % _bins_per_dim[1];
        }
        if(!_reduce_dim[2] && carry) {
            start_z = (start_z+1) % _bins_per_dim[2];
        }
    }

    try {
        std::ofstream file{filename};
        for(double& d : output) {
            file << d << " ";
        }
    } catch (std::ofstream::failure& e) {
        Log::global_log->error() << "[Pressure Profiler] Failed to write pressure measurements to file.\n" << e.what() << std::endl;
        _simulation.exit(-1);
    }
}

/********************************************
 * Direct Profiler code
 * ******************************************/
void DirectProfiler::init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {
    for (int i = 0; i < 3; i++) {
        _bin_widths[i] = domain->getGlobalLength(i) / (double) _bins_per_dim[i];
    }
    _v = _bin_widths[0] * _bin_widths[1] * _bin_widths[2];

    _bin_data.init(_bins_per_dim, _max_samples);

    _sample_counter = 0;
    _write_counter = 0;
}

void DirectProfiler::step(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                          unsigned long simstep) {
    _write_counter = (_write_counter + 1) % _write_gap;
    _sample_counter = (_sample_counter + 1) % _sample_gap;

    if(_write_counter == 0) {
        Vec3D<double> pressureData;
        profilePressure(particleContainer, pressureData);
        std::stringstream name {};
        name << _output_prefix << simstep << ".pm";

        writePressure(name.str(), pressureData);
    }
    else if (_sample_counter == 0) {
        updateTempBuffers(particleContainer);
    }
}

double DirectProfiler::measureTemp(const std::array<double, 3> &low, const std::array<double, 3> &high,
                                     ParticleContainer *particleContainer, BinData& binData) {
    std::array<double, 3> v_mean = computeMeanVelocityStep(low, high, particleContainer, binData);
    double m = _simulation.getEnsemble()->getComponent(0)->m();
    double v = 0;
    #if defined(_OPENMP)
    #pragma omp parallel reduction(+:v)
    #endif
    for (auto it = particleContainer->regionIterator(std::data(low), std::data(high),
                                                     ParticleIterator::ALL_CELLS); it.isValid(); ++it) {
        for(int d = 0; d < 3; d++) {
            v += std::pow(it->v(d) - v_mean[d], 2);
        }
    }

    return v * m / 3;
}

std::array<double, 3> DirectProfiler::computeMeanVelocityStep(const std::array<double, 3> &low, const std::array<double, 3> &high,
                                               ParticleContainer *particleContainer, BinData& binData) {
    //compute mean velocity per dim using SMC method, see: Karimian et al. 2011
    std::array<double, 3> v_avg{0};
    double* v_raw = std::data(v_avg);
    std::size_t count = 0;
    #if defined(_OPENMP)
    #pragma omp parallel reduction(+:v_raw[:3], count)
    #endif
    for (auto it = particleContainer->regionIterator(std::data(low), std::data(high),
                                                     ParticleIterator::ALL_CELLS); it.isValid(); ++it) {
        for(int d = 0; d < 3; d++) {
            v_raw[d] += it->v(d);
        }
        count += 1;
    }
    //TODO check this
    if(count == 0) count = 1;
    binData.mol_counts.insert(count);
    binData.velocities.insert(v_avg);

    //compute CAM average first
    {
        std::size_t steps = binData.mol_counts.size();
        auto d = static_cast<double>(steps<=1?1:steps-1);

        std::array<double, 3> v_cam{ 0 };
        for(std::size_t i = 0; i < steps; i++) {
            const auto& vec = binData.velocities.get(i);
            for(int dim = 0; dim < 3; dim++) {
                v_cam[dim] += vec[dim];
            }
        }
        for(int dim = 0; dim < 3; dim++) {
            v_cam[dim] /= d;
        }

        std::size_t nks = 0;
        for(std::size_t i = 0; i < steps; i++) {
            nks += binData.mol_counts.get(i);
        }
        double nks_d = static_cast<double>(nks) / d;

        for(int dim = 0; dim < 3; dim++) {
            v_cam[dim] /= nks_d;
        }
        binData.v_cam.insert(v_cam);
    }

    //compute SAM over CAM
    {
        std::size_t steps = binData.v_cam.size();
        auto d = static_cast<double>(steps<=1?1:steps-1);

        std::array<double, 3> v_scm{ 0 };
        for(std::size_t i = 0; i < steps; i++) {
            const auto& vec = binData.v_cam.get(i);
            for(int dim = 0; dim < 3; dim++) {
                v_scm[dim] += vec[dim];
            }
        }
        for(int dim = 0; dim < 3; dim++) {
            v_scm[dim] /= d;
        }

        return v_scm;
    }
}

double DirectProfiler::computePotential(const std::array<double, 3> &low, const std::array<double, 3> &high,
                                          ParticleContainer *particleContainer) {
    ParticlePairsHandler* pairsHandler = _simulation.getParticlePairsHandler();
    double r_c = _simulation.getcutoffRadius();
    double r_c2 = r_c * r_c;
    double lj_r_c2 = _simulation.getLJCutoff();
    lj_r_c2 = lj_r_c2 * lj_r_c2;

    double u = 0;
    // Compute inner interactions, i.e., all within the region to be measured
    pairsHandler->init();
#if defined(_OPENMP)
#pragma omp parallel reduction(+:u)
#endif
    for (auto it1 = particleContainer->regionIterator(std::data(low), std::data(high),
                                                      ParticleIterator::ALL_CELLS); it1.isValid(); ++it1) {
        std::array<double, 3> distanceVector{ 0 };

        Molecule& m1 = *it1;
        auto it2 = it1;
        ++it2;
        for (; it2.isValid(); ++it2) {
            Molecule& m2 = *it2;

            mardyn_assert(&m1 != &m2);
            double dd = m2.dist2(m1, std::data(distanceVector));

            if (dd < r_c2) {
                Molecule m1_copy = m1;
                Molecule m2_copy = m2;

                m1_copy.buildOwnSoA();
                m2_copy.buildOwnSoA();

                pairsHandler->processPair(m1_copy, m2_copy, std::data(distanceVector), MOLECULE_MOLECULE, dd, (dd < lj_r_c2));
                m1_copy.calcFM();
                m2_copy.calcFM();

                m1_copy.releaseOwnSoA();
                m2_copy.releaseOwnSoA();

                std::array<double, 3> f = m1_copy.F_arr();
                u += f[0]*distanceVector[0] + f[1]*distanceVector[1] + f[2]*distanceVector[2];
            }
            else u += 0;
        }
    }
    pairsHandler->finish();

    // Compute interactions of molecules within the region with molecules outside, but within r_c
    std::array<double, 3> low_outer = low;
    std::array<double, 3> high_outer = high;
    for(int d = 0; d < 3; d++) {
        low_outer[d] -= r_c;
        high_outer[d] += r_c;
    }
    pairsHandler->init();
#if defined(_OPENMP)
#pragma omp parallel reduction(+:u)
#endif
    for (auto it1 = particleContainer->regionIterator(std::data(low), std::data(high),
                                                      ParticleIterator::ALL_CELLS); it1.isValid(); ++it1) {
        std::array<double, 3> distanceVector{ 0 };

        Molecule& m1 = *it1;
        for (auto it2 = particleContainer->regionIterator(std::data(low_outer), std::data(high_outer),
                                                          ParticleIterator::ALL_CELLS); it2.isValid(); ++it2) {
            Molecule& m2 = *it2;
            if(m2.inBox(std::data(low), std::data(high))) continue;
            mardyn_assert(&m1 != &m2);
            double dd = m2.dist2(m1, std::data(distanceVector));

            if (dd < r_c2) {
                Molecule m1_copy = m1;
                Molecule m2_copy = m2;

                m1_copy.buildOwnSoA();
                m2_copy.buildOwnSoA();

                pairsHandler->processPair(m1_copy, m2_copy, std::data(distanceVector), MOLECULE_MOLECULE, dd, (dd < lj_r_c2));
                m1_copy.calcFM();
                m2_copy.calcFM();

                m1_copy.releaseOwnSoA();
                m2_copy.releaseOwnSoA();

                std::array<double, 3> f = m1_copy.F_arr();
                u += f[0]*distanceVector[0] + f[1]*distanceVector[1] + f[2]*distanceVector[2];
            }
            else u += 0;
        }
    }
    pairsHandler->finish();
    return u;
}

double DirectProfiler::computePressure(const std::array<double, 3> &low, const std::array<double, 3> &high,
                                         ParticleContainer *particleContainer, BinData& binData) {
    double p = 0;
    //kinetic term
    p += measureTemp(low, high, particleContainer, binData) / _v;
    //potential term
    p += 1.0/(6.0*_v) * computePotential(low, high, particleContainer);
    return p;
}

void DirectProfiler::profilePressure(ParticleContainer *particleContainer, Vec3D<double>& pressure_data) {
    pressure_data.dims = _bins_per_dim;
    pressure_data.data.resize(_bins_per_dim[0]*_bins_per_dim[1]*_bins_per_dim[2]);

    std::array<double, 3> low { 0 };
    std::array<double, 3> high { 0 };
    for(std::size_t bin_z = 0; bin_z < _bins_per_dim[2]; bin_z++) {
        low[2] = _bin_widths[2] * static_cast<double>(bin_z);
        high[2] = _bin_widths[2] * static_cast<double>(bin_z+1);
        for(std::size_t bin_y = 0; bin_y < _bins_per_dim[1]; bin_y++) {
            low[1] = _bin_widths[1] * static_cast<double>(bin_y);
            high[1] = _bin_widths[1] * static_cast<double>(bin_y+1);
            for(std::size_t bin_x = 0; bin_x < _bins_per_dim[0]; bin_x++) {
                low[0] = _bin_widths[0] * static_cast<double>(bin_x);
                high[0] = _bin_widths[0] * static_cast<double>(bin_x+1);

                double p = computePressure(low, high, particleContainer, _bin_data.at(bin_x, bin_y, bin_z));
                pressure_data.at(bin_x, bin_y, bin_z) = p;
            }
        }
    }
}

void DirectProfiler::updateTempBuffers(ParticleContainer *particleContainer) {
    std::array<double, 3> low { 0 };
    std::array<double, 3> high { 0 };
    for(std::size_t bin_z = 0; bin_z < _bins_per_dim[2]; bin_z++) {
        low[2] = _bin_widths[2] * static_cast<double>(bin_z);
        high[2] = _bin_widths[2] * static_cast<double>(bin_z+1);
        for(std::size_t bin_y = 0; bin_y < _bins_per_dim[1]; bin_y++) {
            low[1] = _bin_widths[1] * static_cast<double>(bin_y);
            high[1] = _bin_widths[1] * static_cast<double>(bin_y+1);
            for(std::size_t bin_x = 0; bin_x < _bins_per_dim[0]; bin_x++) {
                low[0] = _bin_widths[0] * static_cast<double>(bin_x);
                high[0] = _bin_widths[0] * static_cast<double>(bin_x+1);

                (void) measureTemp(low, high, particleContainer, _bin_data.at(bin_x, bin_y, bin_z));
            }
        }
    }
}

/********************************************
 * MOP Profiler code
 * ******************************************/
void MOPProfiler::init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {
    for (int i = 0; i < 3; i++) {
        _bin_widths[i] = domain->getGlobalLength(i) / (double) _bins_per_dim[i];
    }
    _v = _bin_widths[0] * _bin_widths[1] * _bin_widths[2];

    _bin_data.dims = _bins_per_dim;
    _bin_data.data.resize(_bins_per_dim[0]*_bins_per_dim[1]*_bins_per_dim[2]);
    for(auto& d : _bin_data.data) d.init(_max_samples);

    _sample_counter = 0;
    _write_counter = 0;
}

void MOPProfiler::step(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                       unsigned long simstep) {
    _write_counter = (_write_counter + 1) % _write_gap;
    _sample_counter = (_sample_counter + 1) % _sample_gap;

    if(_write_counter == 0) writeAvg(simstep);
    if(_sample_counter == 0) updatePBuffers(particleContainer);
}

void MOPProfiler::writeAvg(unsigned long simstep) {
    Vec3D<double> p_avg;
    p_avg.dims = _bin_data.dims;
    p_avg.data.resize(_bin_data.data.size(), 0);

    std::size_t steps = _bin_data.data[0].size();
    if(steps <= 0) return;
    for(std::size_t pos_z = 0; pos_z < _bins_per_dim[2]; pos_z++) {
        for(std::size_t pos_y = 0; pos_y < _bins_per_dim[1]; pos_y++) {
            for(std::size_t pos_x = 0; pos_x < _bins_per_dim[0]; pos_x++) {
                double avg = 0;
                auto& bin = _bin_data.at(pos_x,pos_y,pos_z);
                for(std::size_t step = 0; step < steps; step++) {
                    avg += bin.get(step);
                }
                avg /= static_cast<double>(steps > 1 ? (steps-1) : 1);
                p_avg.at(pos_x,pos_y,pos_z) = avg / (3*_v);
            }
        }
    }

    std::stringstream name {};
    name << _output_prefix << simstep << ".pm";
    writePressure(name.str(), p_avg);
}

void MOPProfiler::updatePBuffers(ParticleContainer *particleContainer) {
    std::array<double, 3> low { 0 };
    std::array<double, 3> high { 0 };
    for(std::size_t bin_z = 0; bin_z < _bins_per_dim[2]; bin_z++) {
        low[2] = _bin_widths[2] * static_cast<double>(bin_z);
        high[2] = _bin_widths[2] * static_cast<double>(bin_z+1);
        for(std::size_t bin_y = 0; bin_y < _bins_per_dim[1]; bin_y++) {
            low[1] = _bin_widths[1] * static_cast<double>(bin_y);
            high[1] = _bin_widths[1] * static_cast<double>(bin_y+1);
            for(std::size_t bin_x = 0; bin_x < _bins_per_dim[0]; bin_x++) {
                low[0] = _bin_widths[0] * static_cast<double>(bin_x);
                high[0] = _bin_widths[0] * static_cast<double>(bin_x+1);

                double p = 0;
                p += computeMomentum(low, high, particleContainer);
                p += computePotential(low, high, particleContainer);

                _bin_data.at(bin_x, bin_y, bin_z).insert(p);
            }
        }
    }
}

double MOPProfiler::computeMomentum(const std::array<double, 3> &low, const std::array<double, 3> &high,
                                    ParticleContainer *particleContainer) {
    double p = 0;
#if defined(_OPENMP)
#pragma omp parallel reduction(+:p)
#endif
    for (auto it = particleContainer->regionIterator(std::data(low), std::data(high),
                                                     ParticleIterator::ALL_CELLS); it.isValid(); ++it) {
        double tmp = 0;
        for(int d = 0; d < 3; d++) {
            tmp += it->Vi(d) + it->Vi(d);
        }
        tmp /= it->mass();
        p += tmp;
    }
    return p;
}

double MOPProfiler::computePotential(const std::array<double, 3> &low, const std::array<double, 3> &high,
                                     ParticleContainer *particleContainer) {
    double r_c = _simulation.getcutoffRadius();
    double r_c2 = r_c * r_c;
    double lj_r_c2 = _simulation.getLJCutoff();
    lj_r_c2 = lj_r_c2 * lj_r_c2;
    ParticlePairsHandler* pairsHandler = _simulation.getParticlePairsHandler();
    double u = 0;

    std::array<double, 3> low_outer = low;
    std::array<double, 3> high_outer = high;
    for(int d = 0; d < 3; d++) {
        low_outer[d] -= r_c;
        high_outer[d] += r_c;
    }

    pairsHandler->init();
    #if defined(_OPENMP)
    #pragma omp parallel reduction(+:u)
    #endif
    for (auto it1 = particleContainer->regionIterator(std::data(low_outer), std::data(high_outer),
                                                     ParticleIterator::ALL_CELLS); it1.isValid(); ++it1) {
        std::array<double, 3> distanceVector{ 0 };
        Molecule& m1 = *it1;

        auto it2 = it1;
        ++it2;
        for (; it2.isValid(); ++it2) {
            Molecule& m2 = *it2;
            double tmp = 0;

            mardyn_assert(&m1 != &m2);
            double dd = m2.dist2(m1, std::data(distanceVector));

            if (dd < r_c2) {
                Molecule m1_copy = m1;
                Molecule m2_copy = m2;

                m1_copy.buildOwnSoA();
                m2_copy.buildOwnSoA();

                pairsHandler->processPair(m1_copy, m2_copy, std::data(distanceVector), MOLECULE_MOLECULE, dd, (dd < lj_r_c2));
                m1_copy.calcFM();
                m2_copy.calcFM();

                m1_copy.releaseOwnSoA();
                m2_copy.releaseOwnSoA();

                std::array<double, 3> f_ij = m1_copy.F_arr();
                for(int d = 0; d < 3; d++) {
                    tmp += f_ij[d] * distanceVector[d];
                }

                tmp *= computeIntersection(low, high, m1_copy.r_arr(), m2.r_arr());
                u += tmp;
            }
            else u += 0;
        }
    }
    pairsHandler->finish();

    return u;
}

static inline bool inBox(const std::array<double, 3> &low, const std::array<double, 3> &high, const std::array<double, 3> &p) {
    constexpr double tolerance = 1e-12;
    bool result = true;
    for(int d = 0; d < 3; d++) {
        result &= ((p[d] + tolerance >=  low[d]) || (p[d] - tolerance >=  low[d])) &&
                  ((p[d] + tolerance <= high[d]) || (p[d] - tolerance <= high[d]));
    }
    return result;
}

static inline double dot(const std::array<double, 3> &p1, const std::array<double, 3> &p2) {
    double result = 0;
    for(int d = 0; d < 3; d++) result += p1[d]*p2[d];
    return result;
}

static inline void intersectLinePlane(const std::array<double, 3>& low, const std::array<double, 3>& high,
                                      const std::array<double, 3>& center, const std::array<double, 3>& p,
                                      const std::array<double, 3>& v, std::array<double, 6>& lambdas) {
    using d3 = std::array<double,3>;
    constexpr std::array<d3, 6> normals {d3{1,0,0}, d3{-1,0,0}, d3{0,1,0}, d3{0,-1,0}, d3{0,0,1}, d3{0,0,-1}};
    std::array<d3, 6> bases {center, center, center, center, center, center};
    bases[0][0]=high[0]; bases[1][0]=low[0]; bases[2][1]=high[1];
    bases[3][1]=low[1]; bases[4][2]=high[2]; bases[5][2]=low[2];

    for(int i = 0; i < 6; i++) {
        lambdas[i] = (dot(normals[i], bases[i]) - dot(normals[i], p)) / dot(normals[i], v);
        if(lambdas[i] < 0) lambdas[i] = std::numeric_limits<double>::infinity();
    }
}

static inline std::array<double, 3> projectPointVector(const std::array<double, 3>& p, const std::array<double, 3>& v, double lambda) {
    using d3 = std::array<double,3>;
    d3 x = p;
    for(int d = 0; d < 3; d++) x[d] += lambda * v[d];
    return x;
}

double MOPProfiler::computeIntersection(const std::array<double, 3> &low, const std::array<double, 3> &high,
                                        const std::array<double, 3> &p1, const std::array<double, 3> &p2) {
    using d3 = std::array<double,3>;
    if(inBox(low, high, p1) && inBox(low, high, p2)) return 1;

    d3 center = high;
    for(int d = 0; d<3; d++) center[d] = (center[d] - low[d]) / 2;

    //one in box
    if(inBox(low, high, p1) ^ inBox(low, high, p2))
    {
        d3 p {}, v {};
        if(inBox(low,high,p1)) for(int d = 0; d < 3; d++) { v[d] = p2[d] - p1[d]; p[d] = p1[d]; }
        else for(int d = 0; d < 3; d++) { v[d] = p1[d] - p2[d]; p[d] = p2[d]; }

        std::array<double, 6> lambdas {};
        intersectLinePlane(low, high, center, p, v, lambdas);
        double lambda = *std::min_element(lambdas.begin(), lambdas.end());

#ifndef NDEBUG
        d3 x = projectPointVector(p, v, lambda);
        if(!inBox(low, high, x))
            std::cout << "yeet" << std::endl;
        mardyn_assert(inBox(low, high, x));
        mardyn_assert(lambda > 0);
        mardyn_assert(lambda <= 1);
#endif
        return lambda;
    }

    //both outside box
    {
        //project p1 and p2 along line based on p1p2 or p2p1 vector respectively onto all region
        //get the smallest resulting lambda, while resulting point is in region
        d3 p {}, v{}, x1 {}, x2 {};
        std::array<double, 6> lambdas {};


        //p1 to p2
        for(int d = 0; d < 3; d++) { v[d] = p2[d] - p1[d]; p[d] = p1[d]; }
        intersectLinePlane(low, high, center, p, v, lambdas);

        double min = std::numeric_limits<double>::infinity();
        for(double lambda : lambdas) {
            if(lambda > min) continue;
            d3 x = projectPointVector(p, v, lambda);
            if(!inBox(low, high, x)) continue;
            min = lambda;
        }
        //if we fail to project the point, then intersection with region is 0, leading to 0 contribution
        if(std::isinf(min)) return 0;
        x1 = projectPointVector(p, v, min);

        //p2 to p1
        for(int d = 0; d < 3; d++) { v[d] = p1[d] - p2[d]; p[d] = p2[d]; }
        intersectLinePlane(low, high, center, p, v, lambdas);

        min = std::numeric_limits<double>::infinity();
        for(double lambda : lambdas) {
            if(lambda > min) continue;
            d3 x = projectPointVector(p, v, lambda);
            if(!inBox(low, high, x)) continue;
            min = lambda;
        }
        //if we fail to project the point, then intersection with region is 0, leading to 0 contribution
        if(std::isinf(min)) return 0;
        x2 = projectPointVector(p, v, min);


        //l2 distance between x1 and x2 is contribution
        double result = 0;
        for(int d = 0; d < 3; d++) {
            result += std::pow(x1[d] - x2[d], 2);
        }
        return std::sqrt(result);
    }

}
