/*
 * CommunicationPartner.h
 *
 *  Created on: Nov 16, 2015
 *      Author: tchipevn
 */

#ifndef COMMUNICATIONPARTNER_H_
#define COMMUNICATIONPARTNER_H_

#include "parallel/ParticleData.h"
#include "mpi.h"
#include <vector>

typedef enum {
	LEAVING_AND_HALO_COPIES = 0, /** send process-leaving particles and halo-copies together in one message */
	HALO_COPIES = 1, /** send halo-copies only */
	LEAVING_ONLY = 2 /** send process-leaving particles only */
} MessageType;

class ParticleContainer;

/**
 * (Bi-Directional) MPI Communication Partner.
 */
class CommunicationPartner {
public:
	CommunicationPartner(const int r, const double hLo[3], const double hHi[3], const double bLo[3],
			const double bHi[3], const double sh[3], const int offset[3]);
	CommunicationPartner(const int r);
	CommunicationPartner(const int r, const double leavingLo[3], const double leavingHi[3]);

	CommunicationPartner(const CommunicationPartner& o);

	// TODO: no operator= implemented!

	~CommunicationPartner();

	void initSend(ParticleContainer* moleculeContainer, const MPI_Comm& comm, const MPI_Datatype& type,
			MessageType msgType, bool removeFromContainer = false);

	bool testSend();

	bool iprobeCount(const MPI_Comm& comm, const MPI_Datatype& type);

	bool testRecv(ParticleContainer* moleculeContainer, bool removeRecvDuplicates);

	void initRecv(int numParticles, const MPI_Comm& comm, const MPI_Datatype& type);

	void deadlockDiagnosticSendRecv();
	void deadlockDiagnosticSend();
	void deadlockDiagnosticRecv();

	int getRank() {
		return _rank;
	}

	const int* getOffset() {
		return _offset;
	}

	//! Specifies, whether the communication to the CommunicationPartner is along a shared face (_offset has only one entry != 0)
	//! @return returns whether they are direct face sharing neighbours
	bool isFaceCommunicator() const {
		return (!!_offset[0] + !!_offset[1] + !!_offset[2]) == 1;
	}
	//! @return returns in which direction the face is shared. If it is not a face communicator, -1 is returned
	int getFaceCommunicationDirection() const {
		if (!isFaceCommunicator())
			return -1;
		return !!_offset[1] * 1 + !!_offset[2] * 2;
	}

	void enlargeInOtherDirections(unsigned int d, double enlargement) {
		for (unsigned int d2 = 0; d2 < 3; d2++) {
			if (d2 == d)
				continue;
			_bothLow[d2] -= enlargement;
			_bothHigh[d2] += enlargement;
			_leavingLow[d2] -= enlargement;
			_leavingHigh[d2] += enlargement;
			_copiesLow[d2] -= enlargement;
			_copiesHigh[d2] += enlargement;
		}
	}
private:
	int _rank;
	double _bothLow[3], _bothHigh[3];
	double _leavingLow[3], _leavingHigh[3];
	double _copiesLow[3], _copiesHigh[3];
	double _shift[3]; //! for periodic boundaries
	int _offset[3];

// technical variables
	MPI_Request *_sendRequest, *_recvRequest;
	MPI_Status *_sendStatus, *_recvStatus;
	std::vector<ParticleData> _sendBuf, _recvBuf;
	bool _msgSent, _countReceived, _msgReceived;

};

#endif /* COMMUNICATIONPARTNER_H_ */
