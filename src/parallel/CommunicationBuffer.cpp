/*
 * CommunicationBuffer.cpp
 *
 *  Created on: 16 Oct 2017
 *      Author: tchipevn
 */

#include "CommunicationBuffer.h"
#include "ParticleData.h"

CommunicationBuffer::CommunicationBuffer() {
	// TODO Auto-generated constructor stub

}

CommunicationBuffer::~CommunicationBuffer() {
	// TODO Auto-generated destructor stub
}

size_t CommunicationBuffer::getDynamicSize() {
	return _buffer.capacity() * sizeof(ParticleData);
}

unsigned long CommunicationBuffer::getNumElements() {
	return _buffer.size();
}

const ParticleData& CommunicationBuffer::at(size_t i) const {
	return _buffer.at(i);
}

ParticleData& CommunicationBuffer::at(size_t i) {
	return _buffer.at(i);
}

ParticleData* CommunicationBuffer::getData() {
	return _buffer.data();
}

void CommunicationBuffer::resize(unsigned long numElements) {
	_buffer.resize(numElements);
}

void CommunicationBuffer::clear() {
	_buffer.clear();
}
