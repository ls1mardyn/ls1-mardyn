#include "molecules/mixingrules/MixingRuleBase.h"

#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"

#include <algorithm>

using namespace std;
using Log::global_log;

void MixingRuleBase::readXML(XMLfileUnits& xmlconfig) {
	unsigned int readFirst, readSecond;

	xmlconfig.getNodeValue("@cid1", readFirst);
	xmlconfig.getNodeValue("@cid2", readSecond);

	global_log->info() << "Component id1: " << readFirst << endl;
	global_log->info() << "Component id2: " << readSecond << endl;

	_cid1 = std::min(readFirst, readSecond);
	_cid2 = std::max(readFirst, readSecond);

	// sanity checked in Ensemble::checkMixingRules()
}

