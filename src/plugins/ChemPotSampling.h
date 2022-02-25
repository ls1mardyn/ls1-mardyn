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
    float _binwidth {1.0f};
    uint16_t _numBinsGlobal;
    double _globalBoxLengthY;
    double _slabVolume;

    unsigned long _startSampling {0ul};
    unsigned long _writeFrequency {10000ul};
    unsigned long _stopSampling {1000000000ul};

    // Accumulated over <writefrequency>
    CommVar<std::vector<double>> _chemPotSum;
    CommVar<std::vector<double>> _temperatureSum;
    CommVar<std::vector<unsigned long>> _numMoleculesSum;

    CellProcessor* _cellProcessor;
    Molecule _mTest;

    void resetVectors();

public:
    ChemPotSampling();
	~ChemPotSampling();

    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;

    void readXML (XMLfileUnits& xmlconfig) override;

    void endStep(
            ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
            Domain *domain, unsigned long simstep) override;

    void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain) override {}

    std::string getPluginName()override {return std::string("ChemPotSampling");}

    static PluginBase* createInstance(){return new ChemPotSampling();}

};


#endif //MARDYN_TRUNK_CHEMPOTSAMPLING_H
