#include "molecules/mixingrules/LorentzBerthelot.h"

#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"



void LorentzBerthelotMixingRule::readXML(XMLfileUnits& xmlconfig) {
	MixingRuleBase::readXML(xmlconfig);

	xmlconfig.getNodeValue("eta", _eta);
	xmlconfig.getNodeValue("xi", _xi);
	Log::global_log->info() << "Mixing coefficients: (eta, xi) = (" << _eta << ", " << _xi << ")" << std::endl;
}
