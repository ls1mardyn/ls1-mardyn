#ifndef LORENTZBERTHELOT_H_
#define LORENTZBERTHELOT_H_

#include "molecules/mixingrules/MixingRuleBase.h"

#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"

class LorentzBerthelotMixingRule : public MixingRuleBase {
public:

	void readXML(XMLfileUnits& xmlconfig) override {
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

	// _parameters of MixingRuleBase.cpp contains (eta, xi) in case of this mixing rule
	double getEta() const { return _parameters.at(0); }
	double getXi() const { return _parameters.at(1); }

	void setEta(double eta) { _parameters.at(0) = eta; }
	void setXi(double xi) { _parameters.at(1) = xi; }

	std::string getType() const override { return "LB"; };

};

#endif /* LORENTZBERTHELOT_H_ */
