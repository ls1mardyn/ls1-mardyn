#include "EnsembleBase.h"

#include "molecules/mixingrules/MixingRuleBase.h"
#include "molecules/mixingrules/LorentzBerthelot.h"
#include "utils/xmlfileUnits.h"
#include "utils/Logger.h"
#include "Simulation.h"
#include "Domain.h"

#include <vector>

using namespace std;
using Log::global_log;

Ensemble::~Ensemble() {
	delete _domain;
	_domain = nullptr;
	for(auto& m : _mixingrules)
		delete (m);
}

void Ensemble::readXML(XMLfileUnits& xmlconfig) {
	long numComponents = 0;
	XMLfile::Query query = xmlconfig.query("components/moleculetype");
	numComponents = query.card();
	global_log->info() << "Number of components: " << numComponents << endl;
	if (numComponents == 0) {
		global_log->fatal() << "No components found. Please verify that you have input them correctly." << std::endl;
		Simulation::exit(96123);
	}
	_components.resize(numComponents);
	XMLfile::Query::const_iterator componentIter;
	string oldpath = xmlconfig.getcurrentnodepath();
	for(componentIter = query.begin(); componentIter; componentIter++) {
		xmlconfig.changecurrentnode(componentIter);
		unsigned int cid = 0;
		xmlconfig.getNodeValue("@id", cid);
		_components[cid - 1].readXML(xmlconfig);
		_componentnamesToIds[_components[cid - 1].getName()] = cid - 1;
		global_log->debug() << _components[cid - 1].getName() << " --> " << cid - 1 << endl;
	}
	xmlconfig.changecurrentnode(oldpath);

	/* mixing rules */
	query = xmlconfig.query("components/mixing/rule");
	XMLfile::Query::const_iterator mixingruletIter;
	uint32_t numMixingrules = 0;
	numMixingrules = query.card();
	global_log->info() << "Found " << numMixingrules << " mixing rules." << endl;
	_mixingrules.resize(numMixingrules);

	// data structure for mixing coefficients of domain class (still in use!!!)
	std::vector<double>& dmixcoeff = global_simulation->getDomain()->getmixcoeff();
	dmixcoeff.clear();

	for(mixingruletIter = query.begin(); mixingruletIter; mixingruletIter++) {
		xmlconfig.changecurrentnode(mixingruletIter);
		MixingRuleBase* mixingrule = nullptr;
		string mixingruletype;

		xmlconfig.getNodeValue("@type", mixingruletype);
		global_log->info() << "Mixing rule type: " << mixingruletype << endl;
		if("LB" == mixingruletype) {
			mixingrule = new LorentzBerthelotMixingRule();

		} else {
			global_log->error() << "Unknown mixing rule " << mixingruletype << endl;
			Simulation::exit(1);
		}
		mixingrule->readXML(xmlconfig);
		_mixingrules.push_back(mixingrule);

		/*
		 * Mixing coefficients
		 *
		 * TODO: information of mixing rules (eta, xi) is stored in Domain class and its actually in use
		 * --> we need to decide where this information should be stored in future, in ensemble class,
		 * in the way it is done above?
		 *
		 */
		double xi, eta;
		xmlconfig.getNodeValue("xi", xi);
		xmlconfig.getNodeValue("eta", eta);
		dmixcoeff.push_back(xi);
		dmixcoeff.push_back(eta);
	}
	xmlconfig.changecurrentnode(oldpath);
	setComponentLookUpIDs();
}

void Ensemble::setComponentLookUpIDs() {
	// Get the maximum Component ID.
	unsigned maxID = 0;
	for(auto c = _components.begin(); c != _components.end(); ++c)
		maxID = std::max(maxID, c->ID());

	// we need a look-up table for the Lennard-Jones centers in the Vectorized*CellProcessors
	// (for no other centers)
	unsigned centers = 0;
	for(auto c = _components.begin(); c != _components.end(); ++c) {
		c->setLookUpId(centers);
		centers += c->numLJcenters();
	}
}
