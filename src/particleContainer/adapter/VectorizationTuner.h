/*
 * VectorizationTuner.h
 *
 *  Created on: Jun 12, 2015
 *      Author: tchipevn
 */

#ifndef VECTORIZATIONTUNER_H_
#define VECTORIZATIONTUNER_H_

#include "CellProcessor.h"
#include "FlopCounter.h"
#include <vector>
#include "io/OutputBase.h"
#include "Simulation.h"
#include "ensemble/EnsembleBase.h"

class Component;
class ParticleCell;
//class VectorizedCellProcessor;
//class FlopCounter;

class VectorizationTuner : public OutputBase{
public:

	VectorizationTuner(std::string outputPrefix, double cutoffRadius, double LJCutoffRadius, CellProcessor **cellProcessor);
	VectorizationTuner(double cutoffRadius, double LJCutoffRadius, CellProcessor **cellProcessor);

	~VectorizationTuner();

	//, double& gflopsOwn, double& gflopsPair);
	virtual void initOutput(ParticleContainer* particleContainer,
				DomainDecompBase* domainDecomp, Domain* domain) {
		FlopCounter flopCounter2 = FlopCounter(_cutoffRadius, _LJCutoffRadius);
		tune(*(_simulation.getEnsemble()->components()), **_cellProcessor , flopCounter2);
	}

	virtual void readXML(XMLfileUnits& xmlconfig) {
		_outputPrefix = "mardyn";
		xmlconfig.getNodeValue("outputprefix", _outputPrefix);
		global_log->info() << "Output prefix: " << _outputPrefix << std::endl;
	}

	//do nothing
	virtual void doOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			Domain* domain, unsigned long simstep,
			std::list<ChemicalPotential>* lmu){}

	// do nothing
	virtual void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain){}

	virtual std::string getPluginName() {
		return std::string("VectorizationTuner");
	}

private:
	std::string _outputPrefix;
	CellProcessor** _cellProcessor;
	double _cutoffRadius;
	double _LJCutoffRadius;
	void tune(std::vector<Component> ComponentList, CellProcessor& vcp, FlopCounter& fc);
    void iterate(std::vector<Component> ComponentList,  CellProcessor& vcp, FlopCounter& fc, int numMols, double& gflopsOwn, double& gflopsPair);

	void runOwn(CellProcessor& cp, ParticleCell& cell1, int numRepetitions);
	void runPair(CellProcessor& cp, ParticleCell& cell1, ParticleCell& cell2, int numRepetitions);

	void initSomeMolecules(Component& comp, ParticleCell& cell1, ParticleCell& cell2);
	void initMeshOfMolecules(double boxMin[3], double boxMax[3], Component& comp, ParticleCell& cell1, ParticleCell& cell2);
    void initUniformRandomMolecules(double boxMin[3], double boxMax[3], Component& comp, ParticleCell& cell1, ParticleCell& cell2, int numMols);
    void initNormalRandomMolecules(double boxMin[3], double boxMax[3], Component& comp, ParticleCell& cell1, ParticleCell& cell2, int numMols);
	void clearMolecules(ParticleCell& cell1);
};

#endif /* VECTORIZATIONTUNER_H_ */
