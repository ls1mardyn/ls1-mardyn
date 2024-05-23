/*
 * CapiSampling.h
 *
 *  Created on: May 2024
 *      Author: jniemann
 */
#pragma once

#include <memory>
#include <vector>
#include <string>

#include "PluginBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "utils/CommVar.h"
#include "parallel/DomainDecompBase.h"

/** @brief
* Plugin: can be enabled via config.xml <br>
*
* Samples values binwise in all directions within given Box
* There is NO differentiation between components; it is assumed that all (pseudo) components have the same mass (for correct calculation of temperature)
* \code{.xml}
* <plugin name="CapiSampling">
*           <numBinsX>FLOAT</numBinsX>                  <!-- # of sampling bins in x direction; default 64 -->
*           <numBinsZ>FLOAT</numBinsZ>                  <!-- # of sampling bins in z direction; default 64 -->
            <start>INT</start>                          <!-- Simstep to start sampling; default 0 -->
            <writefrequency>INT</writefrequency>        <!-- Simstep to write out result file; default 10000 -->
            <stop>INT</stop>                            <!-- Simstep to stop sampling; default 1000000000 -->
            <yLower>FLOAT</yLower>                            <!-- lower bound in y-direction -->
            <yUpper>FLOAT</yUpper>                            <!-- upper bound in y-direction -->
* </plugin>
* \endcode
*/
class CapiSampling : public PluginBase {

 private:
    // Control: general
    unsigned int _numBinsX= 64;
    unsigned int _numBinsZ= 64;
    unsigned long _startSampling {0ul};
    unsigned long _writeFrequency {1000ul};
    unsigned long _stopSampling {std::numeric_limits<unsigned long>::max()};
    double _yLower = 0.;
    double _yUpper = {std::numeric_limits<double>::max()};

    // Auxiliary variables
    double _binwidthX;
    double _binwidthZ;
    unsigned long _lenVector;
    std::array<double, 3> _globalBoxLength;
    double _cellVolume;

    // Accumulated quantities over _writeFrequency per bin
    // NOTE: Only the root process knows correct values (except number of molecules and temperature)
    // With only the root process writing data, MPI_Reduce instead of MPI_Allreduce could be used -> only root has correct data
    std::vector<unsigned long> _numMolecules_accum;             // Number of molecules in bin
    std::vector<unsigned long> _doftotal_accum;                 // DOF in bin
    std::vector<double> _mass_accum;                            // Mass
    std::vector<double> _ekin_accum;                            // Kinetic energy including drift
    std::vector<double> _virial_accum;                          // Virial
    std::array<std::vector<double>, 3> _ekinVect_accum;         // Kinetic energy in each direction (drift corrected); radial, height (y), tangantial
    std::array<std::vector<double>, 3> _velocityVect_accum;     // Drift velocity in each direction; radial, height (y), tangantial
    std::array<std::vector<double>, 3> _virialVect_accum;       // Virial in each direction; radial, height (y), tangantial

    std::vector<unsigned long> _countSamples;                   // Number of samples; can vary from bin to bin as some bins could be empty

    void resizeVectors();  // Change size of accumulation vectors
    void resetVectors();   // Set accumulation vectors to zero

 public:
    CapiSampling();
    ~CapiSampling() override = default;

    void init(ParticleContainer* /* particleContainer */, DomainDecompBase* domainDecomp, Domain* domain) override;

    void readXML(XMLfileUnits& xmlconfig) override;

    void beforeForces(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, unsigned long /* simstep */) override {}

    // Needs to be called when halo cells are still existing (for sampling of chemical potential)
    void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
    
    void endStep(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */, unsigned long /* simstep */) override {}
    
    void finish(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */) override {}

    std::string getPluginName() override { return "CapiSampling"; }

    static PluginBase* createInstance() { return new CapiSampling(); }
};

