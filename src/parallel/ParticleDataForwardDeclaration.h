/*
 * ParticleDataForwardDeclaration.h
 *
 *  Created on: 21 Apr 2017
 *      Author: seckler
 */

#pragma once

#ifndef ENABLE_REDUCED_MEMORY_MODE
	class ParticleDataFull;
	typedef ParticleDataFull ParticleData;
#else
	class ParticleDataWR;
	typedef ParticleDataWR ParticleData;
#endif


