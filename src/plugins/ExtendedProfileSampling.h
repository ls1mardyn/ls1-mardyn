/*
 * ExtendedProfileSampling.h
 *
 *  Created on: Feb 2022
 *      Author: homes
 */

#ifndef MARDYN_TRUNK_EXTENDEDPROFILESAMPLING_H
#define MARDYN_TRUNK_EXTENDEDPROFILESAMPLING_H

#include "PluginBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "utils/CommVar.h"
#include "parallel/DomainDecompBase.h"

/** @brief
* Plugin: can be enabled via config.xml <br>
*
* Samples residual chemical potential binwise in y-direction for the LJTS fluid (pot. energy corrections not implemented!) using the Widom insertion method
* \code{.xml}
* <plugin name="ExtendedProfileSampling">
*			<binwidth>FLOAT</binwidth>                  <!-- Width of sampling bins; default 1.0 -->
            <start>INT</start>                          <!-- Simstep to start sampling; default 0 -->
            <writefrequency>INT</writefrequency>        <!-- Simstep to write out result file; default 10000 -->
            <stop>INT</stop>                            <!-- Simstep to stop sampling; default 1000000000 -->
            <singlecomponent>BOOL</singlecomponent>     <!-- Ignore component of the particles; default false -->
            <highermoments>BOOL</highermoments>         <!-- Sample higher moments; default false -->
            <chemicalpotential enable=BOOL>             <!-- Enable or disable sampling of chemical potential; default false -->
                <lattice>BOOL</lattice>                 <!-- Choose if lattice or random insertion; Note: The random method does not take local density into account; default true -->
                <factorNumTest>FLOAT</factorNumTest>    <!-- Factor which specifies number of inserted test particles (numTest = factor*numPartsGlobal); default 4.0 -->
                <samplefrequency>INT</samplefrequency>  <!-- Sampling every INT step; default 50 -->
            </chemicalpotential>
* </plugin>
* \endcode
*/
class ExtendedProfileSampling : public PluginBase{

private:
    // Control: general
    float _binwidth;
    unsigned long _startSampling;
    unsigned long _writeFrequency;
    unsigned long _stopSampling;
    bool _singleComp;

    // Higher moments
    bool _sampleHigherMoms;

    // Control: chemical potential
    bool _sampleChemPot;
    bool _lattice;
    float _factorNumTest;
    unsigned long _samplefrequency;

    // Auxiliary variables
    uint16_t _numBinsGlobal;
    unsigned long _lenVector;
    std::array<double, 3> _globalBoxLength;
    double _slabVolume;
    CellProcessor* _cellProcessor;
    ParticlePairsHandler* _particlePairsHandler;
    Molecule _mTest;

    // Accumulated quantities over _writeFrequency per bin
    // NOTE: Only the root process knows correct values (with exceptions)
    std::vector<unsigned long> _numMolecules_accum;             // Number of molecules in bin
    std::vector<double> _density_accum;                         // Local density
    std::vector<double> _temperature_accum;                     // Temperature (drift corrected)
    std::vector<double> _ekin_accum;                            // Kinetic energy
    std::vector<double> _epot_accum;                            // Potential energy
    std::vector<double> _pressure_accum;                        // Pressure
    std::vector<double> _chemPot_accum;                         // Chemical potential as sampled in ms2 (Widom insertion method)
    std::vector<unsigned long> _countNTest_accum;               // Number of inserted test particles for chem. pot. sampling
    std::array<std::vector<double>, 3> _temperatureVect_accum;  // Kinetic temperature in each direction (drift corrected)
    std::array<std::vector<double>, 3> _velocityVect_accum;     // Drift velocity in each direction
    std::array<std::vector<double>, 3> _pressureVect_accum;     // Pressure in each direction
    std::array<std::vector<double>, 3> _energyfluxVect_accum;   // Energy flux (heat flux plus enthalpy flux) in each direction

    // Higher moments
    std::vector<double> _hmDelta_accum;                               // Higher moment: Delta
    std::array<std::vector<double>, 9> _hmPressure_accum;             // Higher moment: Pressure; cxcx, cxcy, cxcz, cycx, cycy, cycz, czcx, czcy, czcz
    std::array<std::vector<double>, 9> _hmR_accum;                    // Higher moment: R; cxcx, cxcy, cxcz, cycx, cycy, cycz, czcx, czcy, czcz
    std::array<std::vector<double>, 27> _hmM_accum;                    // Higher moment: M; cicxcx, cicxcy, cicxcz, cicycx, cicycy, cicycz, ciczcx, ciczcy, ciczcz mit i = x,y,z

    std::vector<unsigned long> _countSamples;                   // Number of samples; can vary from bin to bin as some bins could be empty

    void resizeVectors(); // Change size of accumulation vectors
    void resetVectors();  // Set accumulation vectors to zero

public:
    ExtendedProfileSampling();
	~ExtendedProfileSampling();

    void init(ParticleContainer* /* particleContainer */, DomainDecompBase* domainDecomp, Domain* domain) override;

    void readXML (XMLfileUnits& xmlconfig) override;

    void beforeForces(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, unsigned long /* simstep */) override {}
    void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
    void endStep(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */, unsigned long /* simstep */) override {}
    void finish(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */) override {}

    std::string getPluginName() override { return std::string("ExtendedProfileSampling"); }

    static PluginBase* createInstance() { return new ExtendedProfileSampling(); }

};


#endif //MARDYN_TRUNK_EXTENDEDPROFILESAMPLING_H
