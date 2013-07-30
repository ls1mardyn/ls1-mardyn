#ifndef MIXINGRULE_BASE_H
#define MIXINGRULE_BASE_H

#include <string>

class XMLfileUnits;

class MixingRuleBase {
public:
	MixingRuleBase(){};
	~MixingRuleBase(){};

	virtual void readXML(XMLfileUnits& xmlconfig);
private:
	std::string _cid1;
	std::string _cid2;
};

#endif /* MixingRuleBase */