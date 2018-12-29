#ifndef ENSEMBLE_BASE_H_
#define ENSEMBLE_BASE_H_

#include "molecules/Component.h"

#include <string>
#include <vector>
#include <map>

#include "Domain.h"

class ParticleContainer;
class MixingRuleBase;
class ChemicalPotential;

//! list of updatable values
enum GlobalVariable {
	NUM_PARTICLES      = 1<<0,
	ENERGY             = 1<<1,
	VOLUME             = 1<<2,
	CHEMICAL_POTENTIAL = 1<<3,
	TEMPERATURE        = 1<<4,
	PRESSURE           = 1<<5
};

class XMLfileUnits;
class DomainBase;

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
	virtual void updateGlobalVariable(ParticleContainer *particleContainer, GlobalVariable variable) = 0;

	Domain *& domain() { return _domain; }
	Component* getComponent(int cid) {
		mardyn_assert(cid < static_cast<int>(_components.size()));
		return &_components.at(cid);
	}
	Component* getComponent(std::string name) { return getComponent(_componentnamesToIds[name]); }
	std::vector<Component>* getComponents() { return &_components; }
	void addComponent(Component& component) { _components.push_back(component); }

	//! prepare the _compIDs used by the Vectorized*CellProcessors
	void setComponentLookUpIDs();

	std::string getType(){return _type;}

	// returns nullptrs for canonical ensemble or gets overridden by GrandCanonical
    virtual std::list<ChemicalPotential>* getLmu(){return nullptr;}

    virtual void initConfigXML(ParticleContainer *moleculeContainer, double h) {};

protected:
	std::vector<Component> _components;
	std::map<std::string,int> _componentnamesToIds;
	std::vector<MixingRuleBase*> _mixingrules;
	Domain *_domain;
	std::string _type = "Undefined";
};

#endif /* ENSEMBLE_BASE_H_ */
