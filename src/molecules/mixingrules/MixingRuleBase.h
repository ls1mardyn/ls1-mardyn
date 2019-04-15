#ifndef MIXINGRULE_BASE_H_
#define MIXINGRULE_BASE_H_

#include <string>

class XMLfileUnits;

class MixingRuleBase {
public:
	MixingRuleBase(){};
	virtual ~MixingRuleBase(){};

	/** Automatically sets cid1 <= cid2.
	 * Validity checks are performed elsewhere.
	 */
	virtual void readXML(XMLfileUnits& xmlconfig);

	unsigned int getCid1() const {
		return _cid1;
	}

	unsigned int getCid2() const {
		return _cid2;
	}

private:
	unsigned int _cid1;
	unsigned int _cid2;
};

#endif /* MIXINGRULE_BASE_H_ */
