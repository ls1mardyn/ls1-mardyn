#pragma once

class VelocityExchangeTest;
#include <string>
#include <map>
#include <unordered_map>
#include <list>
#include <cstdint>
#include <vector>
#include <array>

#include "plugins/PluginBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/CommVar.h"

/** @brief Read in XML configuration for VelocityExchange
 *
 * This plugin can be used to create a temperature gradient by exchanging the velocities of
 * the coldest particle in the specified warm region and the warmest particle in the specified
 * cold region. As temperature gradients are often investigated using a symmetric system, this
 * plugin can be configured to set up a second warm region at the opposite side of the simulation
 * box in y-direction.
 *
 * The following XML object structure is handled by this method:
 * \code{.xml}
    <plugin name="VelocityExchange">
        <control>
            <start>0</start>                     <!-- step to start; default 0 -->
            <frequency>10000</frequency>         <!-- frequency to conduct velo swap; default 10000 -->
            <stop>200000000</stop>               <!-- step to stop; default 200000000 -->
        </control>
        <coldregion>                             <!-- region with lower temperature -->
            <xmin>0</xmin> <xmax>box</xmax>      <!-- range x-axis; default: 0 to box size -->
            <ymin>140</ymin> <ymax>160</ymax>    <!-- range y-axis -->
            <zmin>0</zmin> <zmax>box</zmax>      <!-- range z-axis; default: 0 to box size -->
        </coldregion>
        <warmregion>                             <!-- region with higher temperature -->
            <symmetric>1</symmetric>             <!-- if set to 1 (default), a second warm region is inserted at the opposite side (in y-direction) of the
                                                      simulation box (from yBoxsize-ymax to yBoxsize-ymin) -->
            <xmin>0</xmin> <xmax>box</xmax>      <!-- range x-axis; default: 0 to box size -->
            <ymin>20</ymin> <ymax>30</ymax>      <!-- range y-axis -->
            <zmin>0</zmin> <zmax>box</zmax>      <!-- range z-axis; default: 0 to box size -->
        </warmregion>
    </plugin>
 * \endcode
 */

class VelocityExchange: public PluginBase {
 private:
    friend VelocityExchangeTest;

    struct TimestepControl {
        uint64_t start {0ul};
        uint64_t freq {10000ul};
        uint64_t stop {200000000ul};
    } _control;

    struct Region {
        double min[3] {0.0f, 0.0f, 0.0f};
        double max[3] {0.0f, 0.0f, 0.0f};
    };

    Region _cold_region;
    Region _warm_region;

    bool _symmetry {true};

    double _boxLength[3];

    uint32_t _numComp;

    void exchangeVelocities(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp);

    // find molecule with extreme absolute velocity (warmest or coldest) in the region "begin_iterator"
    void findExtremeMols(const RegionParticleIterator& begin_iterator, const bool flgColdRegion,
                         CommVar<std::vector<double>>& velocity_abs, std::vector<Molecule*>& mol_ptr);

    // assign the (rotational) velocities to the molecule with "molID" if in the region "begin_iterator"
    void assignVelocities(const RegionParticleIterator& begin_iterator, const CommVar<std::vector<unsigned long>>& molID,
                          const CommVar<std::array<std::vector<double>, 3>>& velocity, const CommVar<std::array<std::vector<double>, 3>>& rotVelo);

    // Used in unit test
    // gets minimum and maximum value of specified region in direction d
    void getRegionCoords(double& min_val, double& max_val, unsigned short d, const bool flgColdRegion) {
        if (flgColdRegion) {
            min_val = _cold_region.min[d];
            max_val = _cold_region.max[d];
        } else {
            min_val = _warm_region.min[d];
            max_val = _warm_region.max[d];
        }
    }

 public:
    VelocityExchange();
    ~VelocityExchange() override = default;

    void readXML(XMLfileUnits& xmlconfig) override;

    void init(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */) override;

    void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;

    void endStep(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */, unsigned long /* simstep */) override {}

    void finish(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */) override {}

    std::string getPluginName() override { return std::string("VelocityExchange"); }

    static PluginBase* createInstance() { return new VelocityExchange(); }
};
