/*
 * CoalescenceSampling.h
 *
 *  Created on: May 2024
 *      Author: niemann
 */


#ifndef MARDYN_TRUNK_COALESCENCESAMPLING_H
#define MARDYN_TRUNK_COALESCENCESAMPLING_H

#include <memory>
#include <vector>
#include <string>
#include <limits>

#include "PluginBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "utils/CommVar.h"
#include "parallel/DomainDecompBase.h"


/** @brief
* Plugin: can be enabled via config.xml <br>
*
* Samples values in cylinder binwise in y-direction
* There is NO differentiation between components; it is assumed that all (pseudo) components have the same mass (for correct calculation of temperature)
* \code{.xml}
  <plugin name="CoalescenceSampling">
            <binwidth>FLOAT</binwidth>                  <!-- Width of sampling bins; default 1.0 -->
            <radius>FLOAT</radius>                      <!-- radius of Cylinder -->
            <start>INT</start>                          <!-- Simstep to start sampling; default 0 -->
            <stop>INT</stop>                            <!-- Simstep to stop sampling; default std::numeric_limits<unsigned long>::max() -->
            <writefrequency>INT</writefrequency>        <!-- Simstep to write out result file; default 100 -->
            <yLower>FLOAT</yLower>                      <!-- start of Cylinder -->
            <yUpper>FLOAT</yUpper>                      <!-- end of Cylinder -->
            <outputPrefix>STRING</outputPrefix>                      <!-- outputPrefix -->
  </plugin>
* \endcode
*/

class CoalescenceSampling : public PluginBase{

 private: 
    // Control: general
    float _binwidth {.5f};
    unsigned long _startSampling {0ul};
    unsigned long _stopSampling {std::numeric_limits<unsigned long>::max()};
    unsigned long _writeFrequency {100ul};
    float _yLower {0.f};
    float _yUpper {std::numeric_limits<float>::max()};
    float _radius {4.f};
    std::string _outputPrefix = "coalSampling_";
    std::string _rhoOutputFilename = _outputPrefix +  "_coalSampling_r"+std::to_string(_radius)+".csv";

    // Auxiliary variables
    unsigned int _numBinsGlobal;
    std::array<double, 3> _globalBoxLength;
    
    // Accumulated quantities over _writeFrequency per bin
    // NOTE: Only the root process knows correct values (except number of molecules and temperature)
    // With only the root process writing data, MPI_Reduce instead of MPI_Allreduce could be used -> only root has correct data
    std::vector<unsigned long> _numMolecules_accum;             // Number of molecules in bin
    // std::vector<unsigned long> _doftotal_accum;                 // DOF in bin
    // std::vector<double> _mass_accum;                            // Mass
    // std::vector<double> _ekin_accum;                            // Kinetic energy including drift
    // std::vector<double> _virial_accum;                          // Virial
    // std::array<std::vector<double>, 3> _ekinVect_accum;         // Kinetic energy in each direction (drift corrected); radial, height (y), tangantial
    // std::array<std::vector<double>, 3> _velocityVect_accum;     // Drift velocity in each direction; radial, height (y), tangantial
    // std::array<std::vector<double>, 3> _virialVect_accum;       // Virial in each direction; radial, height (y), tangantial

    std::vector<unsigned long> _countSamples;                   // Number of samples; can vary from bin to bin as some bins could be empty

    void resizeVectors();  // Change size of accumulation vectors
    void resetVectors();   // Set accumulation vectors to zero

 public:
    CoalescenceSampling();
    ~CoalescenceSampling() override = default;

    void init(ParticleContainer* /* particleContainer */, DomainDecompBase* domainDecomp, Domain* domain) override;

    void readXML(XMLfileUnits& xmlconfig) override;

    void beforeForces(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, unsigned long /* simstep */) override {}

    // Needs to be called when halo cells are still existing (for sampling of chemical potential)
    void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
    
    void endStep(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */, unsigned long /* simstep */) override {}
    
    void finish(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */) override {}

    std::string getPluginName() override { return "CoalescenceSampling"; }

    static PluginBase* createInstance() { return new CoalescenceSampling(); }
};



#endif //MARDYN_TRUNK_COALESCENCESAMPLING_H

