#ifndef LORENTZBERTHELOT_H
#define LORENTZBERTHELOT_H

#include "molecules/mixingrules/MixingRuleBase.h"

class XMLfileUnits;

class LorentzBerthelotMixingRule : public MixingRuleBase {

public:
	void readXML(XMLfileUnits& xmlconfig);

private:
	double _eta;
	double _xi;
};

#endif /* LORENTZBERTHELOT_H */