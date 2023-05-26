#include "molecules/mixingrules/MixingRuleBase.h"

#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"



void MixingRuleBase::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("@cid1", _cid1);
	xmlconfig.getNodeValue("@cid2", _cid2);
	Log::global_log->info() << "Component id1: " << _cid1 << std::endl;
	Log::global_log->info() << "Component id2: " << _cid2 << std::endl;
}

