#include "EnsembleBase.h"

#include "molecules/mixingrules/MixingRuleBase.h"
#include "molecules/mixingrules/LorentzBerthelot.h"
#include "utils/xmlfileUnits.h"
#include "utils/Logger.h"
#include "utils/mardyn_assert.h"
#include "Simulation.h"
#include "Domain.h"

#include <vector>
#include <algorithm>
#include <memory>


Ensemble::~Ensemble() {
	delete _domain;
	_domain = nullptr;
}

void Ensemble::readXML(XMLfileUnits& xmlconfig) {
	long numComponents = 0;
	XMLfile::Query query = xmlconfig.query("components/moleculetype");
	numComponents = query.card();
	Log::global_log->info() << "Number of components: " << numComponents << std::endl;
	if (numComponents == 0) {
		Log::global_log->fatal() << "No components found. Please verify that you have input them correctly." << std::endl;
		mardyn_exit(96123);
	}
	_components.resize(numComponents);
	XMLfile::Query::const_iterator componentIter;
	std::string oldpath = xmlconfig.getcurrentnodepath();
	for(componentIter = query.begin(); componentIter; componentIter++) {
		xmlconfig.changecurrentnode(componentIter);
		unsigned int cid = 0;
		xmlconfig.getNodeValue("@id", cid);
		_components[cid - 1].readXML(xmlconfig);
		_componentnamesToIds[_components[cid - 1].getName()] = cid - 1;
		Log::global_log->debug() << _components[cid - 1].getName() << " --> " << cid - 1 << std::endl;
	}
	xmlconfig.changecurrentnode(oldpath);

	/* mixing rules */
	query = xmlconfig.query("components/mixing/rule");
	XMLfile::Query::const_iterator mixingruletIter;
	uint32_t numMixingrules = 0;
	numMixingrules = query.card();
	Log::global_log->info() << "Found " << numMixingrules << " mixing rules." << std::endl;

	for(mixingruletIter = query.begin(); mixingruletIter; mixingruletIter++) {
		xmlconfig.changecurrentnode(mixingruletIter);
		std::shared_ptr<MixingRuleBase> mixingrule;
		std::string mixingruletype;

		xmlconfig.getNodeValue("@type", mixingruletype);
		Log::global_log->info() << "Mixing rule type: " << mixingruletype << std::endl;
		if ("LB" == mixingruletype) {
			mixingrule = std::make_shared<LorentzBerthelotMixingRule>();

		} else {
			Log::global_log->error() << "Unknown mixing rule " << mixingruletype << std::endl;
			mardyn_exit(1);
		}
		mixingrule->readXML(xmlconfig);

		const int cid1 = mixingrule->getCid1();
		const int cid2 = mixingrule->getCid2();
		// Check if cid is larger than number of components
		// cid starts with 0 and cid2 is always larger than cid1
		if (cid2 >= numComponents) {
			Log::global_log->error() << "Mixing: cid=" << cid2+1 << " is larger than number of components ("
									 << numComponents << ")" << std::endl;
			Simulation::exit(1);
		}
		_mixingrules[cid1][cid2] = mixingrule;
	}
	// Use xi=eta=1.0 as default if no rule was specified
	for (int cidi = 0; cidi < numComponents; ++cidi) {
		for (int cidj = cidi+1; cidj < numComponents; ++cidj) {  // cidj is always larger than cidi
			if (_mixingrules[cidi].count(cidj) == 0) {
				// Only LorentzBerthelot is supported until now
				auto mixingrule = std::make_shared<LorentzBerthelotMixingRule>();
				mixingrule->setCid1(cidi);
				mixingrule->setCid2(cidj);
				_mixingrules[cidi][cidj] = mixingrule;
				Log::global_log->warning() << "Mixing coefficients for components "
										   << mixingrule->getCid1()+1 << " + " << mixingrule->getCid2()+1  // +1 due to internal cid
										   << " set to default (LB with xi=eta=1.0)" << std::endl;
			}
		}
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

void Ensemble::setMixingrule(std::shared_ptr<MixingRuleBase> mixingrule) {
	
	// Symmetry for mixing rules is assumed
	// cid1 must be less than cid2
	const auto [cid1, cid2] = std::minmax(mixingrule->getCid1(), mixingrule->getCid2());

	// Check if cids are valid
	if (cid1 == cid2) {
		Log::global_log->error() << "Mixing setMixingrule: cids must not be the same" << std::endl;
		Simulation::exit(1);
	}
	if (std::min(cid1, cid2) < 0) {
		Log::global_log->error() << "Mixing setMixingrule: cids must not be negative" << std::endl;
		Simulation::exit(1);
	}
	if (std::max(cid1, cid2) >= _components.size()) {
		Log::global_log->error() << "Mixing setMixingrule: cids must not exceed number of components ("
								 << _components.size() << ")" << std::endl;
		Simulation::exit(1);
	}
	
	_mixingrules[cid1][cid2] = mixingrule;
}
