/*
 * ParticleContainer.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLECONTAINER_H_
#define DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLECONTAINER_H_

#include "iterators/ParticleIterator.h"

namespace autopas {

template<class Particle, class ParticleCell>
class ParticleContainer {
public:
	typedef ParticleIterator<Particle, ParticleCell> iterator;

	void init(int numCells) {
		_data.resize(numCells);
	}

	void addParticle(Particle& p, int cellNumber) {
		_data.at(cellNumber).addParticle(p);
	}

	iterator begin() {return ParticleIterator<Particle, ParticleCell>(&_data);}
private:
	std::vector<ParticleCell> _data;
};

} /* namespace autopas */

#endif /* DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLECONTAINER_H_ */
