#include "molecules/mixingrules/LorentzBerthelot.h"

#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"

using Log::global_log;

void LorentzBerthelotMixingRule::readXML(XMLfileUnits& xmlconfig) {
	MixingRuleBase::readXML(xmlconfig);

	xmlconfig.getNodeValue("eta", _eta);
	xmlconfig.getNodeValue("xi", _xi);
	global_log->info() << "Mixing coefficients: (eta, xi) = (" << _eta << ", " << _xi << ")" << std::endl;
}
