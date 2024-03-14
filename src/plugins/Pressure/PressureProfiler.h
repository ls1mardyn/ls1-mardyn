//
// Created by alex on 25.02.24.
//

#ifndef MARDYN_PRESSUREPROFILER_H
#define MARDYN_PRESSUREPROFILER_H

#include "plugins/PluginBase.h"
#include "utils/mardyn_assert.h"
#include "DataStructures.h"

#include <array>
#include <vector>

enum PressureMode {
    MOP = 0, Direct, ModeCount
};

struct ProfilerBase {
    virtual ~ProfilerBase() = default;
    virtual void init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) = 0;
    virtual void step(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain, unsigned long simstep) = 0;

    std::array<std::size_t, 3> _bins_per_dim;
    std::array<bool, 3> _reduce_dim;
    std::array<double, 3> _bin_widths;
    double _v;
    std::size_t _max_samples;
    //! @brief Amount of simulation steps until next sampling
    std::size_t _sample_gap;
    std::size_t _sample_counter;
    //! @brief Amount of simulation steps until output file
    std::size_t _write_gap;
    std::size_t _write_counter;
    std::string _output_prefix;

    void writePressure(const std::string &filename, Vec3D<double>& pressure_data);
};

/**
 * Based on: Details about pressure calculation in molecular dynamic analysis by Karimian et al. 2014
 * */
struct DirectProfiler : public ProfilerBase {
    void init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) override;

    void step(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
              unsigned long simstep) override;

    BinVec3D _bin_data;
    /**
     * returns T/(k_b * N_bin)
     * */
    double measureTemp(const std::array<double,3>& low, const std::array<double, 3>& high, ParticleContainer *particleContainer, BinData& binData);
    std::array<double, 3> computeMeanVelocityStep(const std::array<double,3>& low, const std::array<double, 3>& high, ParticleContainer *particleContainer, BinData& binData);
    double computePotential(const std::array<double,3>& low, const std::array<double, 3>& high, ParticleContainer *particleContainer);
    double computePressure(const std::array<double,3>& low, const std::array<double, 3>& high, ParticleContainer *particleContainer, BinData& binData);
    void profilePressure(ParticleContainer *particleContainer, Vec3D<double>& pressure_data);
    void updateTempBuffers(ParticleContainer *particleContainer);
};

/**
 * Based on DOI: 10.1088/0953-8984/24/28/284133
 * */
struct MOPProfiler : public ProfilerBase {
    void init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) override;

    void step(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
              unsigned long simstep) override;

    Vec3D<NBuffer<double>> _bin_data;

    void writeAvg(unsigned long simstep);
    void updatePBuffers(ParticleContainer* particleContainer);
    double computeMomentum(const std::array<double,3>& low, const std::array<double, 3>& high, ParticleContainer *particleContainer);
    double computePotential(const std::array<double,3>& low, const std::array<double, 3>& high, ParticleContainer *particleContainer);
    double computeIntersection(const std::array<double,3>& low, const std::array<double, 3>& high, const std::array<double,3>& p1, const std::array<double, 3>& p2);
};

//! todo We assume single phase for now
class PressureProfiler : public PluginBase {
public:
    void init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) override;

    void readXML(XMLfileUnits &xmlconfig) override;

    void endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                 unsigned long simstep) override;

    void finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) override;

    std::string getPluginName() override;

    static PluginBase* createInstance() {
        return new PressureProfiler();
    }

private:
	std::size_t _start;
    PressureMode _mode;
    std::unique_ptr<ProfilerBase> _profiler;
};


#endif //MARDYN_PRESSUREPROFILER_H
