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
	CommunicationPartner(int r, double hLo[3], double hHi[3], double bLo[3], double bHi[3], double sh[3]);
	CommunicationPartner(int r);
	CommunicationPartner(int r, double leavingLo[3], double leavingHi[3]);

	CommunicationPartner(const CommunicationPartner& o);

	// TODO: no operator= implemented!

	~CommunicationPartner();

	void initSend(
			ParticleContainer* moleculeContainer, const MPI_Comm& comm,
			const MPI_Datatype& type, MessageType msgType, bool removeFromContainer = false);

	bool testSend();

	bool iprobeCount(const MPI_Comm& comm, const MPI_Datatype& type);

	bool testRecv(ParticleContainer* moleculeContainer, bool removeRecvDuplicates);

	void initRecv(int numParticles, const MPI_Comm& comm, const MPI_Datatype& type);

	void deadlockDiagnosticSendRecv();
	void deadlockDiagnosticSend();
	void deadlockDiagnosticRecv();

	int getRank(){
		return _rank;
	}
private:
	int _rank;
	double _bothLow[3], _bothHigh[3];
	double _leavingLow[3], _leavingHigh[3];
	double _copiesLow[3], _copiesHigh[3];
	double _shift[3]; //! for periodic boundaries

	// technical variables
	MPI_Request *_sendRequest, *_recvRequest;
	MPI_Status *_sendStatus, *_recvStatus;
	std::vector<ParticleData> _sendBuf, _recvBuf;
	bool _msgSent, _countReceived, _msgReceived;

};

#endif /* COMMUNICATIONPARTNER_H_ */
