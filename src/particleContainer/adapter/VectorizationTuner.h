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
#include "ensemble/EnsembleBase.h"

class Component;
class ParticleCell;
//class VectorizedCellProcessor;
//class FlopCounter;

/**
 * @brief VectorizationTuner class.
 * This class is used to get detailed information about the performance of the VectorizedCellProcessor.
 * For different scenarios, the performance is evaluated and output.
 * Later this could be used to actually use this class as a tuner, i.e. to use the best possible vectorization method for the actual computation.
 *
 */
class VectorizationTuner : public OutputBase{
public:
	/**
	* @brief Constructor of VectorizationTuner for the old input mode.
	* Here the parameter (outputPrefix) has to be passed explicitly.
	*
	* @param outputPrefix
	* @param cutoffRadius
	* @param LJCutoffRadius
	* @param cellProcessor pointer to the pointer of the cellProcessor. This is needed, since the cell processor is not yet set, when this function is called.
	*/
	VectorizationTuner(std::string outputPrefix, double cutoffRadius, double LJCutoffRadius, CellProcessor **cellProcessor);

	/**
	 * @brief Constructor of VectorizationTuner for the xml input mode.
	 * Here the parameter (outputPrefix) does not have to be passed, it is written from the xml file instead.

	 * @param cutoffRadius
	 * @param LJCutoffRadius
	 * @param cellProcessor pointer to the pointer of the cellProcessor. This is needed, since the cell processor is not yet set, when this function is called.
	 */
	VectorizationTuner(double cutoffRadius, double LJCutoffRadius, CellProcessor **cellProcessor);

	/**
	 * destructor of the VectorizationTuner class.
	 */
	~VectorizationTuner();

	//documentation in OutputBase
	void initOutput(ParticleContainer* particleContainer,
				DomainDecompBase* domainDecomp, Domain* domain);

	//documentation in OutputBase, used to get parameters from xml files.
	void readXML(XMLfileUnits& xmlconfig) {
		_outputPrefix = "mardyn";
		xmlconfig.getNodeValue("outputprefix", _outputPrefix);
		global_log->info() << "Output prefix: " << _outputPrefix << std::endl;
	}

	//documentation in OutputBase, does nothing.
	void doOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			Domain* domain, unsigned long simstep,
			std::list<ChemicalPotential>* lmu){}

	//documentation in OutputBase, does nothing.
	void finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain){}

	//documentation in OutputBase.
	std::string getPluginName() {
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
