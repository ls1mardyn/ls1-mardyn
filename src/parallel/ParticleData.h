#pragma once
/**
 * the old class ParticleData is now called ParticleDataFull.
 */

#include "ParticleDataFull.h"
#include "ParticleDataWR.h"

#ifndef MARDYN_WR
	typedef ParticleDataFull ParticleData;
#else
	typedef ParticleDataWR ParticleData;
#endif
