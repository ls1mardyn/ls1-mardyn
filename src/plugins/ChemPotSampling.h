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
* Samples chemical potential binwise in y-direction ONLY FOR ONE-CENTERED LJ PARTICLES
* \code{.xml}
* <plugin name="ChemPotSampling">
*			<binwidth>1</binwidth> <!-- Width of sampling bins; maximum value: 65534 -->
* </plugin>
* \endcode
*/
class ChemPotSampling : public PluginBase{

private:
    float _binwidth;
    float _factorNumTest;
    unsigned long _startSampling;
    unsigned long _writeFrequency;
    unsigned long _stopSampling;

    uint16_t _numBinsGlobal;
    double _globalBoxLength[3];
    double _slabVolume;

    // Accumulated over _writeFrequency
    CommVar<std::vector<double>> _chemPotSum;
    std::vector<double> _temperatureSumGlobal;
    std::vector<unsigned long> _numMoleculesSumGlobal;

    CellProcessor* _cellProcessor;
    ParticlePairsHandler* _particlePairsHandler;
    Molecule _mTest;

    void resetVectors();

public:
    ChemPotSampling();
	~ChemPotSampling();

    void init(ParticleContainer* /* particleContainer */, DomainDecompBase* domainDecomp, Domain* domain) override;

    void readXML (XMLfileUnits& xmlconfig) override;

    void endStep(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            Domain* /* domain */, unsigned long simstep) override;

    void finish(ParticleContainer* /* particleContainer */,
				DomainDecompBase* /* domainDecomp */, Domain* /* domain */) override {}

    std::string getPluginName() override { return std::string("ChemPotSampling"); }

    static PluginBase* createInstance() { return new ChemPotSampling(); }

};


#endif //MARDYN_TRUNK_CHEMPOTSAMPLING_H
