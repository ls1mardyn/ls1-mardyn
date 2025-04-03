/*
 * ExtendedProfileSampling.h
 *
 *  Created on: Feb 2022
 *      Author: homes
 */

#ifndef MARDYN_TRUNK_EXTENDEDPROFILESAMPLING_H
#define MARDYN_TRUNK_EXTENDEDPROFILESAMPLING_H

class ExtendedProfileSamplingTest;
#include <memory>
#include <vector>
#include <string>
#include <random>

#include "PluginBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "utils/CommVar.h"
#include "parallel/DomainDecompBase.h"

/** @brief
* Plugin: can be enabled via config.xml <br>
*
* Samples residual chemical potential binwise in y-direction for the LJTS fluid (pot. energy corrections not implemented!) using the Widom insertion method
* TODO for general case: Rotation in kin. Temperatur; Upot correction for getEnergy; multi component chem pot;
* \code{.xml}
* <plugin name="ExtendedProfileSampling">
*            <binwidth>FLOAT</binwidth>                  <!-- Width of sampling bins; default 1.0 -->
            <start>INT</start>                          <!-- Simstep to start sampling; default 0 -->
            <writefrequency>INT</writefrequency>        <!-- Simstep to write out result file; default 10000 -->
            <stop>INT</stop>                            <!-- Simstep to stop sampling; default 1000000000 -->
            <singlecomponent>BOOL</singlecomponent>     <!-- Ignore component of the particles; default false -->
            <highermoments>BOOL</highermoments>         <!-- Sample higher moments; default false -->
            <chemicalpotential enable=BOOL>             <!-- Enable or disable sampling of chemical potential; default false -->
                <lattice>BOOL</lattice>                 <!-- Choose if lattice or random insertion; Note: The random method does not take local density into account; default true -->
                <factorNumTest>FLOAT</factorNumTest>    <!-- Factor which specifies number of inserted test particles (numTest = factor*numPartsGlobal); default 4.0 -->
                <samplefrequency>INT</samplefrequency>  <!-- Sampling every INT step; default 50 -->
                <cids>INT,INT,INT,...</cids>            <!-- List of cids to be inserted; starting at 1; default 1 (only first component is inserted) -->
            </chemicalpotential>
* </plugin>
* \endcode
*/
class ExtendedProfileSampling : public PluginBase {

 private:
    friend ExtendedProfileSamplingTest;

    // Control: general
    float _binwidth {1.0f};
    unsigned long _startSampling {0ul};
    unsigned long _writeFrequency {10000ul};
    unsigned long _stopSampling {1000000000ul};
    bool _singleComp {false};

    // Higher moments
    bool _sampleHigherMoms {false};

    // Control: chemical potential
    bool _sampleChemPot {false};
    bool _lattice {true};
    float _factorNumTest {4.0f};
    unsigned long _samplefrequency {50ul};
    // Vector of components to be inserted during chem. pot. sampling; starting at 1
    std::vector<unsigned long> _cidsTest;

    // Auxiliary variables
    unsigned int _numBinsGlobal;
    unsigned long _lenVector;
    std::array<double, 3> _globalBoxLength;
    double _slabVolume;
    CellProcessor* _cellProcessor;
    std::shared_ptr<ParticlePairsHandler> _particlePairsHandler;
    Molecule _mTest;
    unsigned int _numComps;

    // Accumulated quantities over _writeFrequency per bin
    // NOTE: Only the root process knows correct values (except number of molecules and temperature)
    // With only the root process writing data, MPI_Reduce instead of MPI_Allreduce could be used -> only root has correct data
    std::vector<unsigned long> _numMolecules_accum;             // Number of molecules in bin
    std::vector<unsigned long> _doftotal_accum;                 // DOF in bin
    std::vector<double> _mass_accum;                            // Mass
    std::vector<double> _ekin_accum;                            // Kinetic energy including drift
    std::vector<double> _epot_accum;                            // Potential energy
    std::vector<double> _orientation_accum;                     // Orientation order parameter, cf. Eq. 29 in Mecke2001, between z axis of molecule and y plane
    std::vector<double> _virial_accum;                          // Virial
    std::vector<double> _chemPot_accum;                         // Chemical potential as sampled in ms2 (Widom insertion method)
    std::vector<unsigned long> _countNTest_accum;               // Number of inserted test particles for chem. pot. sampling
    std::array<std::vector<double>, 3> _ekinVect_accum;         // Kinetic energy in each direction (drift corrected)
    std::array<std::vector<double>, 3> _velocityVect_accum;     // Drift velocity in each direction
    std::array<std::vector<double>, 3> _virialVect_accum;       // Virial in each direction
    std::array<std::vector<double>, 3> _forceVect_accum;        // Sum of forces on particles in each direction
    std::array<std::vector<double>, 3> _energyfluxVect_accum;   // Energy flux (heat flux plus enthalpy flux) in each direction

    // Higher moments
    std::vector<double> _hmDelta_accum;                         // Higher moment: Delta
    std::array<std::vector<double>, 3> _hmHeatflux_accum;       // Higher moment: Thermal heatflux; x, y, z
    std::array<std::vector<double>, 9> _hmPressure_accum;       // Higher moment: Pressure (traceless); cxcx, cxcy, cxcz, cycx, cycy, cycz, czcx, czcy, czcz
    std::array<std::vector<double>, 9> _hmR_accum;              // Higher moment: R (traceless); cxcx, cxcy, cxcz, cycx, cycy, cycz, czcx, czcy, czcz
    std::array<std::vector<double>, 27> _hmM_accum;             // Higher moment: M (traceless); cicxcx, cicxcy, cicxcz, cicycx, cicycy, cicycz, ciczcx, ciczcy, ciczcz mit i = x,y,z

    std::vector<unsigned long> _countSamples;                   // Number of samples; can vary from bin to bin as some bins could be empty

    // Random number generator
    std::mt19937 _random_generator;  // Mersenne Twister RNG
    std::uniform_real_distribution<> _random_distr;

    void resizeVectors();  // Change size of accumulation vectors
    void resetVectors();   // Set accumulation vectors to zero

 public:
    ExtendedProfileSampling();
    ~ExtendedProfileSampling() override = default;

    void init(ParticleContainer* /* particleContainer */, DomainDecompBase* domainDecomp, Domain* domain) override;

    void readXML(XMLfileUnits& xmlconfig) override;

    void beforeForces(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, unsigned long /* simstep */) override {}

    // Needs to be called when halo cells are still existing (for sampling of chemical potential)
    void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
    
    void endStep(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */, unsigned long /* simstep */) override {}
    
    void finish(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */) override {}

    std::string getPluginName() override { return "ExtendedProfileSampling"; }

    static PluginBase* createInstance() { return new ExtendedProfileSampling(); }

    // Get value of quantity at certain index; Mainly used for unit test
    double getQuantity(DomainDecompBase* domainDecomp, std::string quantityName, unsigned long index, unsigned long simstep);
};


#endif //MARDYN_TRUNK_EXTENDEDPROFILESAMPLING_H
