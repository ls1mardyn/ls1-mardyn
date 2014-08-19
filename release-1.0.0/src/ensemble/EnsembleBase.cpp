#include "EnsembleBase.h"

#include "molecules/mixingrules/MixingRuleBase.h"
#include "molecules/mixingrules/LorentzBerthelot.h"
#include "utils/xmlfileUnits.h"
#include "utils/Logger.h"

using namespace std;
using Log::global_log;

void Ensemble::readXML(XMLfileUnits& xmlconfig) {
	long numComponents = 0;
	XMLfile::Query query = xmlconfig.query("components/moleculetype");
	numComponents = query.card();
	global_log->info() << "Number of components: " << numComponents << endl;
	_components.resize(numComponents);
	XMLfile::Query::const_iterator componentIter;
	string oldpath = xmlconfig.getcurrentnodepath();
	for( componentIter = query.begin(); componentIter; componentIter++) {
		xmlconfig.changecurrentnode( componentIter );
		unsigned int cid = 0;
		xmlconfig.getNodeValue( "@id", cid );
		_components[cid - 1].readXML(xmlconfig);
		_componentnamesToIds[_components[cid - 1].getName()] = cid - 1;
		global_log->debug() << _components[cid - 1].getName() << " --> " << cid - 1 << endl;
	}
	xmlconfig.changecurrentnode(oldpath);

	/* mixing rules */
	query = xmlconfig.query("components/mixing/rule");
	XMLfile::Query::const_iterator mixingruletIter;
	long numMixingrules = 0;
	numMixingrules = query.card();
	_mixingrules.resize(numMixingrules);

	for( mixingruletIter = query.begin(); mixingruletIter; mixingruletIter++ ) {
		xmlconfig.changecurrentnode( mixingruletIter );
		MixingRuleBase *mixingrule;
		string mixingruletype;

		xmlconfig.getNodeValue("@type", mixingruletype);
		global_log->info() << "Mixing rule type: " << mixingruletype << endl;
		if( "LB" == mixingruletype ) {
			mixingrule = new LorentzBerthelotMixingRule();

		}
		else {
			global_log->error() << "Unknown mixing rule " << mixingruletype << endl;
			exit(1);
		}
		mixingrule->readXML(xmlconfig);
		_mixingrules.push_back(mixingrule);
	}
	xmlconfig.changecurrentnode(oldpath);
}
