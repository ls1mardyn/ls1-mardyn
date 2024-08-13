/*
 * SphericalSampling.h
 *
 *  Created on: July 24
 *      Author: JakNiem
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
* Samples values binwise in radial and y-direction
* There is NO differentiation between components; it is assumed that all (pseudo) components have the same mass (for correct calculation of temperature)
* \code{.xml}
* <plugin name="SphericalSampling">
*           <nShells>INT</nShells>                  <!-- Width of sampling bins; default 1.0 -->
            <start>INT</start>                          <!-- Simstep to start sampling; default 0 -->
            <writefrequency>INT</writefrequency>        <!-- Simstep to write out result file; default 10000 -->
            <stop>INT</stop>                            <!-- Simstep to stop sampling; default MAX -->
* </plugin>
* \endcode
*/

class SphericalSampling : public PluginBase {
 private:
    // Control: general
    unsigned int _nShells {60};
    unsigned long _startSampling {0ul};
    unsigned long _writeFrequency {1000ul};
    unsigned long _stopSampling {std::numeric_limits<unsigned long>::max()};

    // Auxiliary variables
    double _shellMappingDistMax {1.0f}; // for using the squared dist from center instead of dist from center to assign shells
    unsigned long _lenVector {_nShells + 1}; // outer-most element for residual volume (cuboid minus sphere)
    std::array<double, 3> _globalBoxLength;
    std::array<double, 3> _globalCenter;
    double _distMax;

    // Accumulated quantities over _writeFrequency per bin
    // NOTE: Only the root process knows correct values (except number of molecules and temperature)
    // With only the root process writing data, MPI_Reduce instead of MPI_Allreduce could be used -> only root has correct data
    
    
    std::vector<unsigned long> _numMolecules_accum;             // Number of molecules in bin
    std::vector<unsigned long> _doftotal_accum;                 // DOF in bin
    std::vector<double> _mass_accum;                            // Mass
    std::vector<double> _ekin_accum;                            // Kinetic energy including drift
    std::vector<double> _virSph_accum;                          // Virial (sum over Spherical Coord)
    std::vector<double> _virN_accum;                            // Virial in normal diection
    std::vector<double> _virT_accum;                            // Virial in tangential diection
    std::vector<double> _velocityN_accum;                            // Virial in tangential diection
    std::array<std::vector<double>, 3> _virialVect_accum;       // Virial in x,y,z direction
    std::array<std::vector<double>, 3> _ekinVect_accum;         // Kinetic energy in each direction (drift corrected); radial, height (y), tangantial

    std::vector<unsigned long> _countSamples;                   // Number of samples; can vary from bin to bin as some bins could be empty

    void resizeVectors();  // Change size of accumulation vectors
    void resetVectors();   // Set accumulation vectors to zero

 public:
    SphericalSampling();
    ~SphericalSampling() override = default;

    void init(ParticleContainer* /* particleContainer */, DomainDecompBase* domainDecomp, Domain* domain) override;

    void readXML(XMLfileUnits& xmlconfig) override;

    void beforeForces(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, unsigned long /* simstep */) override {}

    // Needs to be called when halo cells are still existing (for sampling of chemical potential)
    void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
    
    void endStep(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */, unsigned long /* simstep */) override {}
    
    void finish(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */) override {}

    std::string getPluginName() override { return "SphericalSampling"; }

    static PluginBase* createInstance() { return new SphericalSampling(); }
};

