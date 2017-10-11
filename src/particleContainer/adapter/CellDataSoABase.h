/*
 * CellDataSoABase.h
 *
 *  Created on: 21 Jan 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_ADAPTER_CELLDATASOABASE_H_
#define SRC_PARTICLECONTAINER_ADAPTER_CELLDATASOABASE_H_

#include <cstddef>

class CellDataSoABase {
public:
	// no virtual destructor - objects are never deleted through a pointer to base
	// they only exist as references FullParticleCell and ParticleCellRMM.
	void setMolNum(size_t molNum) {_molNum = molNum;}
	size_t getMolNum() const {return _molNum;}
	void incrementMolNum() {++_molNum;}
	void decrementMolNum() {--_molNum;}

private:
	size_t _molNum;
};

#endif /* SRC_PARTICLECONTAINER_ADAPTER_CELLDATASOABASE_H_ */
