/*
 * ChemPotSampling.h
 *
 *  Created on: Feb 2022
 *      Author: homes
 */

#ifndef MARDYN_TRUNK_CHEMPOTSAMPLING_H
#define MARDYN_TRUNK_CHEMPOTSAMPLING_H

#include "PluginBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "utils/CommVar.h"
#include "parallel/DomainDecompBase.h"

/** @brief
* Plugin: can be enabled via config.xml <br>
*
* Samples residual chemical potential binwise in y-direction for one component systems using the Widom insertion method
* \code{.xml}
* <plugin name="ChemPotSampling">
*			<binwidth>FLOAT</binwidth>              <!-- Width of sampling bins; default 1.0 -->
            <lattice>BOOL</lattice>                 <!-- Choose if lattice or random insertion; Note: The random method does not take local density into account; default true -->
            <factorNumTest>FLOAT</factorNumTest>    <!-- Factor which specifies number of inserted test particles (numTest = factor*numPartsGlobal); default 4.0 -->
            <start>INT</start>                      <!-- Simstep to start sampling; default 0 -->
            <writefrequency>INT</writefrequency>    <!-- Simstep to write out result file; default 10000 -->
            <samplefrequency>INT</samplefrequency>  <!-- Sampling every INT step; default 50 -->
            <stop>INT</stop>                        <!-- Simstep to stop sampling; default 1000000000 -->
* </plugin>
* \endcode
*/
class ChemPotSampling : public PluginBase{

private:
    float _binwidth;
    float _factorNumTest;
    unsigned long _startSampling;
    unsigned long _writeFrequency;
    unsigned long _samplefrequency;
    unsigned long _stopSampling;
    bool _lattice;

    uint16_t _numBinsGlobal;
    std::array<double, 3> _globalBoxLength;
    double _slabVolume;

    // Accumulated quantities over _writeFrequency
    CommVar<std::vector<double>> _chemPotSum;
    std::vector<double> _temperatureSumGlobal;
    std::vector<double> _temperatureWithDriftSumGlobal;
    std::vector<unsigned long> _numMoleculesSumGlobal;
    CommVar<std::vector<unsigned long>> _countNTest;
    std::vector<unsigned long> _countSamples;

    CellProcessor* _cellProcessor;
    ParticlePairsHandler* _particlePairsHandler;
    Molecule _mTest;

    void resetVectors();

public:
    ChemPotSampling();
	~ChemPotSampling();

    void init(ParticleContainer* /* particleContainer */, DomainDecompBase* domainDecomp, Domain* domain) override;

    void readXML (XMLfileUnits& xmlconfig) override;

    void beforeForces(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, unsigned long /* simstep */) override {}
    void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
    void endStep(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */, unsigned long /* simstep */) override {}
    void finish(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */) override {}

    std::string getPluginName() override { return std::string("ChemPotSampling"); }

    static PluginBase* createInstance() { return new ChemPotSampling(); }

};


#endif //MARDYN_TRUNK_CHEMPOTSAMPLING_H
