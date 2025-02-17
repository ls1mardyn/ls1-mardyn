/*
 * CommunicationBuffer.h
 *
 *  Created on: 16 Oct 2017
 *      Author: tchipevn
 */
#ifdef ENABLE_MPI
#ifndef SRC_PARALLEL_COMMUNICATIONBUFFER_H_
#define SRC_PARALLEL_COMMUNICATIONBUFFER_H_

#include "molecules/MoleculeForwardDeclaration.h"
#include "utils/mardyn_assert.h"

#include <vector>
#include <stddef.h>
#include <mpi.h>

// do not uncomment the if, it will break halo copies of the kddecomposition!
//#if (not defined(NDEBUG))
#define LS1_SEND_UNIQUE_ID_FOR_HALO_COPIES
//#pragma message "Compilation info: Unique IDs of Halo-Molecules are always present."
//#endif

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
 *
 * Stores two unsigned long integers, then leaving molecules, then halo molecules.
 */
class CommunicationBuffer {

	friend class CommunicationBufferTest;

public:
	CommunicationBuffer() {
		clear();
	}
	size_t getDynamicSize();

	void clear();

	void resizeForAppendingLeavingMolecules(unsigned long numMols);
	void resizeForAppendingHaloMolecules(unsigned long numMols);
        void resizeForAppendingForceMolecules(unsigned long numMols);

	unsigned char * getDataForSending();
	size_t getNumElementsForSending();
	void resizeForRawBytes(unsigned long numBytes);

	// write
	void addLeavingMolecule(size_t indexOfMolecule, const Molecule& m);
	void addHaloMolecule(size_t indexOfMolecule, const Molecule& m);
        void addForceMolecule(size_t indexOfMolecule, const Molecule& m);

	// read
	void readLeavingMolecule(size_t indexOfMolecule, Molecule& m) const;
	void readHaloMolecule(size_t indexOfMolecule, Molecule& m) const;
	void readForceMolecule(size_t indexOfMolecule, Molecule& m) const;

	void resizeForReceivingMolecules(unsigned long& numLeaving, unsigned long& numHalo);
	void resizeForReceivingMolecules(unsigned long& numForces);

	size_t getNumHalo() const {
		return _numHalo;
	}

	size_t getNumLeaving() const {
		return _numLeaving;
	}

        size_t getNumForces() const {
            return _numForces;
        }

	static MPI_Datatype getMPIDataType() {
		return MPI_CHAR;
	}

private:
	static size_t _numBytesHalo;
	static size_t _numBytesLeaving;
        static size_t _numBytesForces; // where is this set?

	enum class ParticleType_t {HALO=0, LEAVING=1, FORCE=3};
	size_t getStartPosition(ParticleType_t type, size_t indexOfMolecule) const;

	/**
	 * @return the next index for writing
	 */
	template<typename T>
	size_t emplaceValue(size_t indexInBytes, T passByValue);

	template<typename T>
	size_t readValue(size_t indexInBytes, T& passByReference) const;

	typedef unsigned char byte_t;
	std::vector<byte_t> _buffer;
	size_t _numLeaving, _numHalo, _numForces;
};

template<typename T>
inline size_t CommunicationBuffer::emplaceValue(size_t indexInBytes, T passByValue) {
	const size_t numBytesOfT = sizeof(T);
	size_t ret = indexInBytes + numBytesOfT;

	mardyn_assert(_buffer.size() >= ret);

	const byte_t * pointer = reinterpret_cast<byte_t *> (&passByValue);
	for (size_t i = 0; i < numBytesOfT; ++i) {
		_buffer[indexInBytes + i] = pointer[i];
	}

	return ret;
}

template<typename T>
inline size_t CommunicationBuffer::readValue(size_t indexInBytes, T& passByReference) const {
	const size_t numBytesOfT = sizeof(T);
	size_t ret = indexInBytes + numBytesOfT;

	mardyn_assert(_buffer.size() >= ret);

	byte_t * pointer = reinterpret_cast<byte_t *> (&passByReference);
	for (size_t i = 0; i < numBytesOfT; ++i) {
		pointer[i] = _buffer[indexInBytes + i];
	}

	return ret;
}


#endif /* SRC_PARALLEL_COMMUNICATIONBUFFER_H_ */
#endif /* ENABLE_MPI */
