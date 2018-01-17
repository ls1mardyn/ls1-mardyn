/*
 * SoA.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_SOA_H_
#define DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_SOA_H_

#include <vector>

namespace autopas {

class SoA {
public:
	std::vector<double> _positions;
	std::vector<double> _forces;
};

} /* namespace autopas */

#endif /* DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_SOA_H_ */
