/*
 * CommunicationBuffer.h
 *
 *  Created on: 16 Oct 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARALLEL_COMMUNICATIONBUFFER_H_
#define SRC_PARALLEL_COMMUNICATIONBUFFER_H_

#include <vector>
#include <stddef.h>
#include "ParticleDataForwardDeclaration.h"


/**
 * This class enables to send only position-data when sending HALO molecules
 * and all data, when sending LEAVING molecules.
 * Sends two integers first, telling how many LEAVING and how many HALO molecules
 * are being sent, and then all the data.
 *
 * Converts all to CHAR internally.
 *
 * TODO: test how this will work when going to Big Endian/Little Endian architectures,
 * due to CHAR conversion.
 */
class CommunicationBuffer {
public:
	CommunicationBuffer();
	~CommunicationBuffer();

	// At the moment, just wrap old functionality.
	size_t getDynamicSize();

	unsigned long getNumElements();

	const ParticleData& at(size_t i) const;

	ParticleData& at(size_t i);

	ParticleData * getData();

	void resize(unsigned long numElements);

	void clear();

private:
	std::vector<ParticleData> _buffer;
};

#endif /* SRC_PARALLEL_COMMUNICATIONBUFFER_H_ */
