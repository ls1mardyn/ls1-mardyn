#ifndef LORENTZBERTHELOT_H_
#define LORENTZBERTHELOT_H_

#include "molecules/mixingrules/MixingRuleBase.h"

class XMLfileUnits;

class LorentzBerthelotMixingRule : public MixingRuleBase {

public:
	void readXML(XMLfileUnits& xmlconfig);

private:
	double _eta;
	double _xi;
};

#endif /* LORENTZBERTHELOT_H_ */
