#include "LorentzBerthelot.h"

#include <ostream>

#include "MixingRuleBase.h"
#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"

void LorentzBerthelotMixingRule::readXML(const XMLfileUnits& xmlconfig) {
	MixingRuleBase::readXML(xmlconfig);
	double eta, xi;
	xmlconfig.getNodeValue("eta", eta);
	xmlconfig.getNodeValue("xi", xi);
	_parameters = {eta, xi};
	Log::global_log->info() << "Mixing coefficients for components "
							<< this->getCid1()+1 << " + " << this->getCid2()+1  // +1 due to internal cid
							<< ": (eta, xi) = (" << _parameters.at(0)
							<< ", " <<  _parameters.at(1) << ")" << std::endl;
}
