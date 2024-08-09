#ifndef LORENTZBERTHELOT_H_
#define LORENTZBERTHELOT_H_

#include <string>

#include "MixingRuleBase.h"

class XMLfileUnits;

class LorentzBerthelotMixingRule : public MixingRuleBase {
 public:
	void readXML(const XMLfileUnits& xmlconfig) override;

	// _parameters of MixingRuleBase.cpp contains (eta, xi) in case of this mixing rule
	double getEta() const { return _parameters.at(0); }
	double getXi() const { return _parameters.at(1); }

	void setEta(double eta) { _parameters.at(0) = eta; }
	void setXi(double xi) { _parameters.at(1) = xi; }

	std::string getType() const override { return "LB"; }
};

#endif /* LORENTZBERTHELOT_H_ */
