#ifndef MIXINGRULE_BASE_H_
#define MIXINGRULE_BASE_H_

#include <string>
#include <vector>

class XMLfileUnits;

class MixingRuleBase {
 public:
	virtual ~MixingRuleBase() = default;

	virtual void readXML(const XMLfileUnits& xmlconfig);

	int getCid1() const { return _cid1; }
	int getCid2() const { return _cid2; }

	void setCid1(int cid1) { _cid1 = cid1; }  // Internal cid starting with 0
	void setCid2(int cid2) { _cid2 = cid2; }  // Internal cid starting with 0

	virtual std::string getType() const = 0;

	const std::vector<double> &getParameters() const { return _parameters; }
	void setParameters(std::vector<double> params) { _parameters = params; }

	// Parameters of mixing rule
	// e.g. for LB it contains (eta, xi)
	std::vector<double> _parameters {1.0, 1.0};  // default for eta, xi;

 private:
	// Internal cids starting with 0
	int _cid1 {};
	int _cid2 {};
};

#endif /* MIXINGRULE_BASE_H_ */
