#ifndef LORENTZBERTHELOT_H_
#define LORENTZBERTHELOT_H_

#include "molecules/mixingrules/MixingRuleBase.h"

class XMLfileUnits;

class LorentzBerthelotMixingRule : public MixingRuleBase {

public:
	virtual ~LorentzBerthelotMixingRule() {}
	void readXML(XMLfileUnits& xmlconfig);

	double getEta() const {
		return _eta;
	}

	double getXi() const {
		return _xi;
	}

private:
	double _eta;
	double _xi;
};

#endif /* LORENTZBERTHELOT_H_ */
