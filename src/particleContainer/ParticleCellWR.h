/*
 * ParticleCellWR.h
 *
 *  Created on: 20 Jan 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_PARTICLECELLWR_H_
#define SRC_PARTICLECONTAINER_PARTICLECELLWR_H_

#include "ParticleCellBase.h"
#include "particleContainer/adapter/CellDataSoA_WR.h"

class ParticleCell_WR: public ParticleCellBase {
public:
	ParticleCell_WR();
	~ParticleCell_WR();

	/**
	 * \brief Structure of arrays for VectorizedCellProcessor.
	 * \author Johannes Heckl
	 */
	CellDataSoA_WR _cellDataSoA;
};

#endif /* SRC_PARTICLECONTAINER_PARTICLECELLWR_H_ */
