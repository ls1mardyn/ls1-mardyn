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
	virtual ~CellDataSoABase() {}
	virtual size_t getMolNum() const = 0;
};

#endif /* SRC_PARTICLECONTAINER_ADAPTER_CELLDATASOABASE_H_ */
