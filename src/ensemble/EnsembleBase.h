#pragma once

#include "molecules/Component.h"

#include <string>
#include <vector>
#include <map>

#include "DomainBase.h"
#include "molecules/MoleculeForwardDeclaration.h"

class ParticleContainer;

class MixingRuleBase;

class ChemicalPotential;

class DomainDecompBase;

class CellProcessor;

class Domain;

//! list of updatable values
enum GlobalVariable {
	NUM_PARTICLES = 1 << 0,
	ENERGY = 1 << 1,
	VOLUME = 1 << 2,
	CHEMICAL_POTENTIAL = 1 << 3,
	TEMPERATURE = 1 << 4,
	PRESSURE = 1 << 5
};

enum Type {
	undefined,
	NVT,
	muVT
};

class XMLfileUnits;

//! @brief Base class for ensembles
//! @author Christoph Niethammer <niethammer@hlrs.de>
//! 
//! Each ensemble should provide access to extensive (NVE) and intensive 
//! (\mu p t) variables as well as a function to update global variables.
class Ensemble {
public:
	Ensemble() :
			_domain(nullptr) {
	}

	virtual ~Ensemble();

	virtual void readXML(XMLfileUnits& xmlconfig);

	//! @brief Returns the global number of Molecules of the ensemble.
	virtual unsigned long N() = 0;

	//! @brief Returns the global volume of the ensemble
	virtual double V() = 0;

	//! @brief Returns the global energy of the ensemble
	virtual double E() = 0;

	//! @brief Returns the global chemical potential of the ensemble
	virtual double mu() = 0;

	//! @brief Returns the global presure of the ensemble.
	virtual double p() = 0;

	//! @brief Returns the global Temperature of the ensemble.
	virtual double T() = 0;

	//! @brief Calculate global variables
	//! @param variable Variable to be updated.
	virtual void updateGlobalVariable(ParticleContainer* particleContainer, GlobalVariable variable) = 0;

	DomainBase*& domain() { return _domain; }

	Component* getComponent(int cid) {
		mardyn_assert(cid < static_cast<int>(_components.size()));
		return &_components.at(cid);
	}

	Component* getComponent(std::string name) { return getComponent(_componentnamesToIds[name]); }

	std::vector<Component>* getComponents() { return &_components; }

	void addComponent(Component& component) { _components.push_back(component); }

	//! prepare the _compIDs used by the Vectorized*CellProcessors
	void setComponentLookUpIDs();

	/*! get Ensemble Type (NVT or muVT) */
	int getType() { return _type; }

	/*! Returns _lmu pointer for processing by external plugins */
	virtual std::list<ChemicalPotential>* getLmu() { return nullptr; }

	/*! runs steps only needed in GrandCanonicalEnsemble, does nothing for canonical */
	virtual void initConfigXML(ParticleContainer* moleculeContainer) {};

	/*! runs steps only needed in GrandCanonicalEnsemble, does nothing for canonical */
	virtual void prepare_start() {};

	/*! runs simulate step needed in GrandCanonical, nothing for canonical */
	virtual void
	beforeEventNewTimestep(ParticleContainer* moleculeContainer, DomainDecompBase* domainDecomposition,
						   unsigned long simstep) {};

	/*! runs after forces step for GrandCanonical, nothing for canonical */
	virtual void
	afterForces(ParticleContainer* /* moleculeContainer */, DomainDecompBase* /* domainDecomposition */,
				CellProcessor* /* cellProcessor */,
				unsigned long /* simstep */) {};

	/*! runs before temperature control is applied, but after force calculations */
	virtual void beforeThermostat(unsigned long simstep, unsigned long initStatistics) {};

	/*! Store Sample molecule from old input readers in lmu */
	virtual void storeSample(Molecule* m, uint32_t componentid) {};

protected:


	std::vector<Component> _components;
	std::map<std::string, int> _componentnamesToIds;
	std::vector<MixingRuleBase*> _mixingrules;
	DomainBase* _domain;
	Type _type = undefined;

	/* EnsembleBase has a DomainBase pointer _domain, however this is not compatible
	 * with the needed Domain* for the ChemicalPotential/radial function
	 * and for whatever reason, Domain does not inherit from DomainBase.
	*/
	Domain* _simulationDomain;
};