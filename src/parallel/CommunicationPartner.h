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
	LEAVING_AND_HALO_COPIES = 0, 	/** send process-leaving particles and halo-copies together in one message */
	HALO_COPIES = 1, 				/** send halo-copies only */
	LEAVING_ONLY = 2 				/** send process-leaving particles only */
} MessageType;

class ParticleContainer;

/**
 * (Bi-Directional) MPI Communication Partner.
 */
class CommunicationPartner {
public:
	CommunicationPartner(int r, double hLo[3], double hHi[3], double bLo[3], double bHi[3], double sh, int signDir);

	CommunicationPartner(const CommunicationPartner& o);

	// TODO: no operator= implemented!

	~CommunicationPartner();

	void initCommunication(unsigned short d,
			ParticleContainer* moleculeContainer, const MPI_Comm& comm,
			const MPI_Datatype& type, MessageType msgType);

	bool testSend();

	bool iprobeCount(const MPI_Comm& comm, const MPI_Datatype& type);

	bool testRecv(ParticleContainer* moleculeContainer);

	void deadlockDiagnostic();

private:
	int _rank;
	double _regionLow[3], _regionHigh[3];
	double _haloLow[3], _haloHigh[3];
	double _boundaryLow[3], _boundaryHigh[3];

	double _shift; //! for periodic boundaries
	int _signedDirection;

	// technical variables
	MPI_Request *_sendRequest, *_recvRequest;
	MPI_Status *_sendStatus, *_recvStatus;
	std::vector<ParticleData> _sendBuf, _recvBuf;
	bool _msgSent, _countReceived, _msgReceived;

};

#endif /* COMMUNICATIONPARTNER_H_ */
