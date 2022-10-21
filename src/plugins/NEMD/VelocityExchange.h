#ifndef VELOCITYEXCHANGE_H_
#define VELOCITYEXCHANGE_H_

#include <string>
#include <map>
#include <unordered_map>
#include <list>
#include <cstdint>
#include <vector>
#include <array>

#include "plugins/PluginBase.h"
#include "utils/CommVar.h"

/** @brief Read in XML configuration for VelocityExchange
 *
 * The following XML object structure is handled by this method:
 * \code{.xml}
    <plugin name="VelocityExchange">
        <control>
            <start>0</start>                     <!-- step to start; default 0 -->
            <frequency>10000</frequency>         <!-- frequency to conduct velo swap; default 10000 -->
            <stop>200000000</stop>               <!-- step to stop; default 200000000 -->
        </control>
        <coldrange>                              <!-- region with lower temperature -->
            <xmin>0</xmin> <xmax>100</xmax>      <!-- range x-axis -->
            <ymin>140</ymin> <ymax>160</ymax>    <!-- range y-axis -->
            <zmin>0</zmin> <zmax>100</zmax>      <!-- range z-axis -->
        </coldrange>
        <warmrange>                              <!-- region with higher temperature -->
            <symmetric>1</symmetric>             <!-- 0: no symmetry; 1: symmetry in y direction -->
            <xmin>0</xmin> <xmax>100</xmax>      <!-- range x-axis -->
            <ymin>20</ymin> <ymax>30</ymax>      <!-- range y-axis -->
            <zmin>0</zmin> <zmax>100</zmax>      <!-- range z-axis -->
        </warmrange>
    </plugin>
 * \endcode
 */

class VelocityExchange: public PluginBase {
 private:
    struct TimestepControl {
        uint64_t start {0ul};
        uint64_t freq {10000ul};
        uint64_t stop {200000000ul};
    } _control;

    struct C_Range {
        double xmin{0.0f};
        double xmax{0.0f};
        double ymin{0.0f};
        double ymax{0.0f};
        double zmin{0.0f};
        double zmax{0.0f};
    } _cold_range;

    struct W_Range {
        double xmin{0.0f};
        double xmax{0.0f};
        double ymin{0.0f};
        double ymax{0.0f};
        double zmin{0.0f};
        double zmax{0.0f};
    } _warm_range;

    bool _symmetry {true};

    double _boxLength[3];

    uint32_t _numComp;

    void exchangeVelocities(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp);

 public:
    VelocityExchange();
    ~VelocityExchange() override = default;

    void readXML(XMLfileUnits& xmlconfig) override;

    void init(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* domain) override;

    void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;

    void endStep(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */, unsigned long /* simstep */) override {}

    void finish(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */) override {}

    std::string getPluginName() override { return std::string("VelocityExchange"); }

    static PluginBase* createInstance() { return new VelocityExchange(); }
};

#endif /*VELOCITYEXCHANGE_H_*/
