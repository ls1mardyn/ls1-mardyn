/*
 * CommunicationBuffer.cpp
 *
 *  Created on: 16 Oct 2017
 *      Author: tchipevn
 */

#include "CommunicationBuffer.h"
#include "molecules/Molecule.h"
#include "Simulation.h"
#include "utils/mardyn_assert.h"
#include "ensemble/EnsembleBase.h"

#include <climits> /* UINT64_MAX */

#ifdef ENABLE_REDUCED_MEMORY_MODE
// position, velocity, id
size_t CommunicationBuffer::_numBytesLeaving = 6 * sizeof(vcp_real_calc) + sizeof(unsigned long);
// position, id (for now)
size_t CommunicationBuffer::_numBytesHalo = 3 * sizeof(vcp_real_calc)
	#ifdef LS1_SEND_UNIQUE_ID_FOR_HALO_COPIES
			+ sizeof(unsigned long)
	#endif
		;
size_t CommunicationBuffer::_numBytesForces = sizeof(unsigned long) + 3 * sizeof(vcp_real_calc) + 3 * sizeof(vcp_real_accum);
#else
// position, velocity, orientation, angular momentum, id, cid
size_t CommunicationBuffer::_numBytesLeaving = 13 * sizeof(double) + sizeof(unsigned long) + sizeof(int);
// position, orientation, id, cid
size_t CommunicationBuffer::_numBytesHalo = 7 * sizeof(double) + sizeof(int)
	#ifdef LS1_SEND_UNIQUE_ID_FOR_HALO_COPIES
			+ sizeof(unsigned long)
	#endif
		;
size_t CommunicationBuffer::_numBytesForces = sizeof(unsigned long) + 12 * sizeof(double);
#endif


unsigned char* CommunicationBuffer::getDataForSending() {
	return _buffer.data();
}

size_t CommunicationBuffer::getNumElementsForSending() {
	return _buffer.size();
}

void CommunicationBuffer::clear() {
	_numHalo = 0;
	_numLeaving = 0;
	_numForces = 0;
	_buffer.clear();
}

void CommunicationBuffer::resizeForRawBytes(unsigned long numBytes) {
	_buffer.reserve(numBytes);
	_buffer.resize(numBytes);
}

void CommunicationBuffer::resizeForReceivingMolecules(unsigned long& numLeaving, unsigned long& numHalo) { // adjust for force exchange?
	// message has been received

	// read _numHalo and _numLeaving
	size_t i_runningByte = 0;
	i_runningByte = readValue(i_runningByte, _numLeaving);
	i_runningByte = readValue(i_runningByte, _numHalo);
	numLeaving = _numLeaving;
	numHalo = _numHalo;
}

void CommunicationBuffer::resizeForReceivingMolecules(unsigned long& numForces) {
	// message has been received

	// read _numForces
	size_t i_runningByte = 0;
	//i_runningByte = readValue(i_runningByte, _numForces);
	_numForces = _buffer.size() / _numBytesForces;
	numForces = _numForces;
}

void CommunicationBuffer::resizeForAppendingLeavingMolecules(unsigned long numLeaving) {
	_numLeaving += numLeaving;
	mardyn_assert(_numHalo == 0ul); // assumption: add leaving, add leaving, then add halo, halo, halo, ... but not intertwined.
	size_t numBytes = sizeof(_numHalo) + sizeof(_numLeaving) +
				_numLeaving * _numBytesLeaving +
				_numHalo * _numBytesHalo;
	resizeForRawBytes(numBytes);

	// store _numLeaving
	size_t i_runningByte = 0;
	i_runningByte = emplaceValue(i_runningByte, _numLeaving);
}

void CommunicationBuffer::resizeForAppendingHaloMolecules(unsigned long numHalo) {
	// _numLeaving stays
	_numHalo += numHalo;
	size_t numBytes = sizeof(_numHalo) + sizeof(_numLeaving) +
				_numLeaving * _numBytesLeaving +
				_numHalo * _numBytesHalo;
	resizeForRawBytes(numBytes);

	// store _numHalo
	size_t i_runningByte = sizeof(_numLeaving);
	i_runningByte = emplaceValue(i_runningByte, _numHalo);
}

void CommunicationBuffer::resizeForAppendingForceMolecules(unsigned long numForces) {
	_numForces += numForces;
	// maybe some assert
	size_t numBytes = sizeof(_numForces) + _numForces * _numBytesForces;
	resizeForRawBytes(numBytes);

	// Do NOT write the number of force molecules, here! It is assumed at other places, that they are NOT exchanged!
}

void CommunicationBuffer::addLeavingMolecule(size_t indexOfMolecule, const Molecule& m) {
	mardyn_assert(indexOfMolecule < _numLeaving);

	size_t i_firstByte = getStartPosition(ParticleType_t::LEAVING, indexOfMolecule);
	mardyn_assert(i_firstByte + _numBytesLeaving <= _buffer.capacity());

	size_t i_runningByte = i_firstByte;
#ifdef ENABLE_REDUCED_MEMORY_MODE
	// add id, r, v
	i_runningByte = emplaceValue(i_runningByte, m.getID());
	i_runningByte = emplaceValue(i_runningByte, static_cast<vcp_real_calc>(m.r(0)));
	i_runningByte = emplaceValue(i_runningByte, static_cast<vcp_real_calc>(m.r(1)));
	i_runningByte = emplaceValue(i_runningByte, static_cast<vcp_real_calc>(m.r(2)));
	i_runningByte = emplaceValue(i_runningByte, static_cast<vcp_real_calc>(m.v(0)));
	i_runningByte = emplaceValue(i_runningByte, static_cast<vcp_real_calc>(m.v(1)));
	i_runningByte = emplaceValue(i_runningByte, static_cast<vcp_real_calc>(m.v(2)));
#else
	i_runningByte = emplaceValue(i_runningByte, m.getID());
	i_runningByte = emplaceValue(i_runningByte, m.componentid());
	i_runningByte = emplaceValue(i_runningByte, m.r(0));
	i_runningByte = emplaceValue(i_runningByte, m.r(1));
	i_runningByte = emplaceValue(i_runningByte, m.r(2));
	i_runningByte = emplaceValue(i_runningByte, m.v(0));
	i_runningByte = emplaceValue(i_runningByte, m.v(1));
	i_runningByte = emplaceValue(i_runningByte, m.v(2));
	i_runningByte = emplaceValue(i_runningByte, m.q().qw());
	i_runningByte = emplaceValue(i_runningByte, m.q().qx());
	i_runningByte = emplaceValue(i_runningByte, m.q().qy());
	i_runningByte = emplaceValue(i_runningByte, m.q().qz());
	i_runningByte = emplaceValue(i_runningByte, m.D(0));
	i_runningByte = emplaceValue(i_runningByte, m.D(1));
	i_runningByte = emplaceValue(i_runningByte, m.D(2));
#endif

	mardyn_assert(i_runningByte - i_firstByte == _numBytesLeaving);
}

void CommunicationBuffer::addHaloMolecule(size_t indexOfMolecule, const Molecule& m) {
	mardyn_assert(indexOfMolecule < _numHalo);

	size_t i_firstByte = getStartPosition(ParticleType_t::HALO, indexOfMolecule);
	mardyn_assert(i_firstByte + _numBytesHalo <= _buffer.capacity());

	size_t i_runningByte = i_firstByte;
#ifdef ENABLE_REDUCED_MEMORY_MODE
	#ifdef LS1_SEND_UNIQUE_ID_FOR_HALO_COPIES
		i_runningByte = emplaceValue(i_runningByte, m.getID());
	#endif /*LS1_SEND_UNIQUE_ID_FOR_HALO_COPIES*/
	i_runningByte = emplaceValue(i_runningByte, static_cast<vcp_real_calc>(m.r(0)));
	i_runningByte = emplaceValue(i_runningByte, static_cast<vcp_real_calc>(m.r(1)));
	i_runningByte = emplaceValue(i_runningByte, static_cast<vcp_real_calc>(m.r(2)));
#else
	#ifdef LS1_SEND_UNIQUE_ID_FOR_HALO_COPIES
		i_runningByte = emplaceValue(i_runningByte, m.getID());
	#endif /*LS1_SEND_UNIQUE_ID_FOR_HALO_COPIES*/
	i_runningByte = emplaceValue(i_runningByte, m.componentid());
	i_runningByte = emplaceValue(i_runningByte, m.r(0));
	i_runningByte = emplaceValue(i_runningByte, m.r(1));
	i_runningByte = emplaceValue(i_runningByte, m.r(2));
	i_runningByte = emplaceValue(i_runningByte, m.q().qw());
	i_runningByte = emplaceValue(i_runningByte, m.q().qx());
	i_runningByte = emplaceValue(i_runningByte, m.q().qy());
	i_runningByte = emplaceValue(i_runningByte, m.q().qz());
#endif

	mardyn_assert(i_runningByte - i_firstByte == _numBytesHalo);
}

void CommunicationBuffer::addForceMolecule(size_t indexOfMolecule, const Molecule& m) {
	// some MarDynAssert?
	size_t i_firstByte = getStartPosition(ParticleType_t::FORCE, indexOfMolecule);  // adjust getStartPosition etc.
	// some MarDynAssert?

	size_t i_runningByte = i_firstByte;

	// add force molecule
#ifdef ENABLE_REDUCED_MEMORY_MODE
	i_runningByte = emplaceValue(i_runningByte, m.getID());
	i_runningByte = emplaceValue(i_runningByte, static_cast<vcp_real_calc>(m.r(0)));
	i_runningByte = emplaceValue(i_runningByte, static_cast<vcp_real_calc>(m.r(1)));
	i_runningByte = emplaceValue(i_runningByte, static_cast<vcp_real_calc>(m.r(2)));
	i_runningByte = emplaceValue(i_runningByte, static_cast<vcp_real_accum>(m.F(0)));
	i_runningByte = emplaceValue(i_runningByte, static_cast<vcp_real_accum>(m.F(1)));
	i_runningByte = emplaceValue(i_runningByte, static_cast<vcp_real_accum>(m.F(2)));
#else
	i_runningByte = emplaceValue(i_runningByte, m.getID());
	i_runningByte = emplaceValue(i_runningByte, m.r(0));
	i_runningByte = emplaceValue(i_runningByte, m.r(1));
	i_runningByte = emplaceValue(i_runningByte, m.r(2));
	i_runningByte = emplaceValue(i_runningByte, m.F(0));
	i_runningByte = emplaceValue(i_runningByte, m.F(1));
	i_runningByte = emplaceValue(i_runningByte, m.F(2));
	i_runningByte = emplaceValue(i_runningByte, m.M(0));
	i_runningByte = emplaceValue(i_runningByte, m.M(1));
	i_runningByte = emplaceValue(i_runningByte, m.M(2));
	i_runningByte = emplaceValue(i_runningByte, m.Vi(0));
	i_runningByte = emplaceValue(i_runningByte, m.Vi(1));
	i_runningByte = emplaceValue(i_runningByte, m.Vi(2));
#endif
}

void CommunicationBuffer::readLeavingMolecule(size_t indexOfMolecule, Molecule& m) const {
	mardyn_assert(indexOfMolecule < _numLeaving);

	size_t i_firstByte = getStartPosition(ParticleType_t::LEAVING, indexOfMolecule);
	mardyn_assert(i_firstByte + _numBytesLeaving <= _buffer.capacity());

	size_t i_runningByte = i_firstByte;
#ifdef ENABLE_REDUCED_MEMORY_MODE
	unsigned long idbuf;
	vcp_real_calc rbuf[3], vbuf[3];
	i_runningByte = readValue(i_runningByte, idbuf);
	i_runningByte = readValue(i_runningByte, rbuf[0]);
	i_runningByte = readValue(i_runningByte, rbuf[1]);
	i_runningByte = readValue(i_runningByte, rbuf[2]);
	i_runningByte = readValue(i_runningByte, vbuf[0]);
	i_runningByte = readValue(i_runningByte, vbuf[1]);
	i_runningByte = readValue(i_runningByte, vbuf[2]);
	m.setid(idbuf);
	for (int d = 0; d < 3; ++d) {
		m.setr(d, rbuf[d]);
		m.setv(d, vbuf[d]);
	}
#else
	unsigned long idbuf;
	unsigned int cidbuf;
	double rbuf[3], vbuf[3], qbuf[4], Dbuf[3];
	i_runningByte = readValue(i_runningByte, idbuf);
	i_runningByte = readValue(i_runningByte, cidbuf);
	i_runningByte = readValue(i_runningByte, rbuf[0]);
	i_runningByte = readValue(i_runningByte, rbuf[1]);
	i_runningByte = readValue(i_runningByte, rbuf[2]);
	i_runningByte = readValue(i_runningByte, vbuf[0]);
	i_runningByte = readValue(i_runningByte, vbuf[1]);
	i_runningByte = readValue(i_runningByte, vbuf[2]);
	i_runningByte = readValue(i_runningByte, qbuf[0]);
	i_runningByte = readValue(i_runningByte, qbuf[1]);
	i_runningByte = readValue(i_runningByte, qbuf[2]);
	i_runningByte = readValue(i_runningByte, qbuf[3]);
	i_runningByte = readValue(i_runningByte, Dbuf[0]);
	i_runningByte = readValue(i_runningByte, Dbuf[1]);
	i_runningByte = readValue(i_runningByte, Dbuf[2]);
	m.setid(idbuf);
	Component* component = _simulation.getEnsemble()->getComponent(cidbuf);
	m.setComponent(component);
	m = Molecule(idbuf, component,
		rbuf[0], rbuf[1], rbuf[2],
		vbuf[0], vbuf[1], vbuf[2],
		qbuf[0], qbuf[1], qbuf[2], qbuf[3],
		Dbuf[0], Dbuf[1], Dbuf[2]
	);
#endif

	mardyn_assert(i_runningByte - i_firstByte == _numBytesLeaving);
}

void CommunicationBuffer::readHaloMolecule(size_t indexOfMolecule, Molecule& m) const {
	mardyn_assert(indexOfMolecule < _numHalo);

	size_t i_firstByte = getStartPosition(ParticleType_t::HALO, indexOfMolecule);
	mardyn_assert(i_firstByte + _numBytesHalo <= _buffer.capacity());

	// add id, r, v
	size_t i_runningByte = i_firstByte;
#ifdef ENABLE_REDUCED_MEMORY_MODE
	unsigned long idbuf;
	vcp_real_calc rbuf[3];
	#ifdef LS1_SEND_UNIQUE_ID_FOR_HALO_COPIES
		i_runningByte = readValue(i_runningByte, idbuf);
	#else
		idbuf = UINT64_MAX;
	#endif /*LS1_SEND_UNIQUE_ID_FOR_HALO_COPIES*/
	i_runningByte = readValue(i_runningByte, rbuf[0]);
	i_runningByte = readValue(i_runningByte, rbuf[1]);
	i_runningByte = readValue(i_runningByte, rbuf[2]);
	m.setid(idbuf);
	for (int d = 0; d < 3; ++d) {
		m.setr(d, rbuf[d]);
	}
#else
	unsigned long idbuf;
	unsigned int cidbuf;
	double rbuf[3], vbuf[3] = {0., 0., 0.}, qbuf[4], Dbuf[3] = {0., 0., 0.};
	#ifdef LS1_SEND_UNIQUE_ID_FOR_HALO_COPIES
		i_runningByte = readValue(i_runningByte, idbuf);
	#else
		idbuf = UINT64_MAX;
	#endif /*LS1_SEND_UNIQUE_ID_FOR_HALO_COPIES*/
	i_runningByte = readValue(i_runningByte, cidbuf);
	i_runningByte = readValue(i_runningByte, rbuf[0]);
	i_runningByte = readValue(i_runningByte, rbuf[1]);
	i_runningByte = readValue(i_runningByte, rbuf[2]);
	i_runningByte = readValue(i_runningByte, qbuf[0]);
	i_runningByte = readValue(i_runningByte, qbuf[1]);
	i_runningByte = readValue(i_runningByte, qbuf[2]);
	i_runningByte = readValue(i_runningByte, qbuf[3]);
	Component* component = _simulation.getEnsemble()->getComponent(cidbuf);
	m = Molecule(idbuf, component,
		rbuf[0], rbuf[1], rbuf[2],
		vbuf[0], vbuf[1], vbuf[2],
		qbuf[0], qbuf[1], qbuf[2], qbuf[3],
		Dbuf[0], Dbuf[1], Dbuf[2]
	);
#endif

	mardyn_assert(i_runningByte - i_firstByte == _numBytesHalo);
}

void CommunicationBuffer::readForceMolecule(size_t indexOfMolecule, Molecule& m) const {
	// some mardyn assert
	size_t i_firstByte = getStartPosition(ParticleType_t::FORCE, indexOfMolecule);
	// some mardyn assert

	size_t i_runningByte = i_firstByte;

#ifdef ENABLE_REDUCED_MEMORY_MODE
	vcp_real_calc rbuf[3];
	vcp_real_accum Fbuf[3];
	unsigned long idbuf;

	i_runningByte = readValue(i_runningByte, idbuf);
	i_runningByte = readValue(i_runningByte, rbuf[0]);
	i_runningByte = readValue(i_runningByte, rbuf[1]);
	i_runningByte = readValue(i_runningByte, rbuf[2]);
	i_runningByte = readValue(i_runningByte, Fbuf[0]);
	i_runningByte = readValue(i_runningByte, Fbuf[1]);
	i_runningByte = readValue(i_runningByte, Fbuf[2]);
	m.setid(idbuf);
	for(int d = 0; d < 3; d++) {
		m.setr(d, rbuf[d]);
	}
	for(int d = 0; d < 3; d++) {
		m.setF(d, Fbuf[d]);
	}

#else
	double rbuf[3], Fbuf[3], Mbuf[3], Vibuf[3];
	unsigned long idbuf;

	i_runningByte = readValue(i_runningByte, idbuf);
	i_runningByte = readValue(i_runningByte, rbuf[0]);
	i_runningByte = readValue(i_runningByte, rbuf[1]);
	i_runningByte = readValue(i_runningByte, rbuf[2]);
	i_runningByte = readValue(i_runningByte, Fbuf[0]);
	i_runningByte = readValue(i_runningByte, Fbuf[1]);
	i_runningByte = readValue(i_runningByte, Fbuf[2]);
	i_runningByte = readValue(i_runningByte, Mbuf[0]);
	i_runningByte = readValue(i_runningByte, Mbuf[1]);
	i_runningByte = readValue(i_runningByte, Mbuf[2]);
	i_runningByte = readValue(i_runningByte, Vibuf[0]);
	i_runningByte = readValue(i_runningByte, Vibuf[1]);
	i_runningByte = readValue(i_runningByte, Vibuf[2]);
	m.setid(idbuf);
	for(int d = 0; d < 3; d++) {
		m.setr(d, rbuf[d]);
	}
	m.setF(Fbuf);
	m.setM(Mbuf);
	m.setVi(Vibuf);
#endif
	// some mardyn assert
}

size_t CommunicationBuffer::getStartPosition(ParticleType_t type, size_t indexOfMolecule) const {
	size_t ret = 0;

	// two unsigned longs
	ret += sizeof(_numLeaving) + sizeof(_numHalo);

	if(type == ParticleType_t::LEAVING) {
		ret += indexOfMolecule * _numBytesLeaving;
	} else if(type == ParticleType_t::HALO) {
		ret += _numLeaving * _numBytesLeaving + indexOfMolecule * _numBytesHalo;
	} else if(type == ParticleType_t::FORCE) {
		// the number of forces is NOT stored in the buffer, as they are sent on their own!
		ret = indexOfMolecule * _numBytesForces;
		//additionally no leaving or halo molecules are transmitted, thus _numLeaving, _numHalo are zero!
	}

	return ret;
}

size_t CommunicationBuffer::getDynamicSize() {
	return _buffer.capacity() * sizeof(byte_t);
}


