#pragma once
/**
 * the old class ParticleData is now called ParticleDataFull.
 */

#include "ParticleDataFull.h"
#include "ParticleDataWR.h"

#ifndef ENABLE_REDUCED_MEMORY_MODE
	typedef ParticleDataFull ParticleData;
#else
	typedef ParticleDataWR ParticleData;
#endif
