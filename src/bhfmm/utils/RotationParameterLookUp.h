/*
 * RotationParameterLookUp.h
 *
 *  Created on: Dec 15, 2014
 *      Author: uwe
 */

#ifndef ROTATIONPARAMETERLOOKUP_H_
#define ROTATIONPARAMETERLOOKUP_H_

#include <cassert>

namespace bhfmm {

class RotationParameterLookUp {
public:
	RotationParameterLookUp(unsigned ord) : numSlices(ord+1){
		totalNumEntries = numSlices*(numSlices+1)*(2*numSlices+1)/6;
		table = new double[totalNumEntries];
	}

	~RotationParameterLookUp() {
		delete[] table;
	}

	inline double acc_c(unsigned l, int m, int k) const {
		return table[index(l,m,k)];
	}

	inline double& acc(unsigned l, int m, int k) const {
		return table[index(l,m,k)];
	}

	void initFromFile(char* filename);
	void initFromDirEval();
	void initFromRecursion();

	void print(int maxl);

	static RotationParameterLookUp* tab;
private:

	inline unsigned index(unsigned l, int m, int k) const {
		m = abs(m);
		k = abs(k); // Symmetries in m and k!
		assert(l <= numSlices);
		assert(static_cast<unsigned>(m) <= l);
		assert(static_cast<unsigned>(k) <= l);
		return l*(l+1)*(2*l+1)/6 + (l+1)*m + k; // valid for l < 1290
	}

	double* table;
	unsigned totalNumEntries;
	unsigned numSlices;
};

} // namespace bhfmm

#endif /* ROTATIONPARAMETERLOOKUP_H_ */
