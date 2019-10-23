/*
 * CommunicationPartner.h
 *
 *  Created on: Nov 16, 2015
 *      Author: tchipevn
 */

#ifndef COMMUNICATIONPARTNER_H_
#define COMMUNICATIONPARTNER_H_

#include <mpi.h>
#include <vector>
#include <stddef.h>
#include "CommunicationBuffer.h"

typedef enum {
	LEAVING_AND_HALO_COPIES = 0, /** send process-leaving particles and halo-copies together in one message */
	HALO_COPIES = 1, /** send halo-copies only */
	LEAVING_ONLY = 2, /** send process-leaving particles only */
	FORCES = 3 /** send forces */
} MessageType;

class ParticleContainer;

struct PositionInfo {
	double _bothLow[3], _bothHigh[3];
	double _leavingLow[3], _leavingHigh[3];
	double _copiesLow[3], _copiesHigh[3];
	double _shift[3]; //! for periodic boundaries
	int _offset[3];
	bool _enlarged[3][2];
};


/**
 * (Bi-Directional) MPI Communication Partner.
 */
class CommunicationPartner {
public:
	CommunicationPartner(const int r, const double hLo[3], const double hHi[3], const double bLo[3],
			const double bHi[3], const double sh[3], const int offset[3], const bool enlarged[3][2]);
	CommunicationPartner(const int r);
	CommunicationPartner(const int r, const double leavingLo[3], const double leavingHi[3]);

	CommunicationPartner(const CommunicationPartner& o);

	CommunicationPartner() = delete;

	CommunicationPartner& operator =(const CommunicationPartner& b);

	~CommunicationPartner();

	void initSend(ParticleContainer* moleculeContainer, const MPI_Comm& comm, const MPI_Datatype& type,
				  MessageType msgType, std::vector<Molecule>& invalidParticles, bool mightUseInvalidParticles,
				  bool doHaloPositionCheck, bool removeFromContainer = false);

	bool testSend();

	void resetReceive();

	bool iprobeCount(const MPI_Comm& comm, const MPI_Datatype& type);

	bool testRecv(ParticleContainer* moleculeContainer, bool removeRecvDuplicates, bool force = false);

	void initRecv(int numParticles, const MPI_Comm& comm, const MPI_Datatype& type);

	void deadlockDiagnosticSendRecv();
	void deadlockDiagnosticSend();
	void deadlockDiagnosticRecv();

	int getRank() const {
		return _rank;
	}

	const int* getOffset() {
		return _haloInfo[0]._offset;
	}

	//! Specifies, whether the communication to the CommunicationPartner is along a shared face (_offset has only one entry != 0)
	//! @return returns whether they are direct face sharing neighbours
	bool isFaceCommunicator() const {
		return (!!_haloInfo[0]._offset[0] + !!_haloInfo[0]._offset[1] + !!_haloInfo[0]._offset[2]) == 1;
	}
	//! @return returns in which direction the face is shared. If it is not a face communicator, -1 is returned
	int getFaceCommunicationDirection() const {
		if (!isFaceCommunicator())
			return -1;
		return !!_haloInfo[0]._offset[1] * 1 + !!_haloInfo[0]._offset[2] * 2;
	}

	void enlargeInOtherDirections(unsigned int d, double enlargement) {
		for (unsigned int p = 0; p < _haloInfo.size(); p++) {
			for (unsigned int d2 = 0; d2 < 3; d2++) {
				if (d2 == d)
					continue;
				if (!_haloInfo[p]._enlarged[d2][0]) {
					_haloInfo[p]._bothLow[d2] -= enlargement;
					_haloInfo[p]._leavingLow[d2] -= enlargement;
					_haloInfo[p]._copiesLow[d2] -= enlargement;
					_haloInfo[p]._enlarged[d2][0] = true;
				}
				if (!_haloInfo[p]._enlarged[d2][1]) {
					_haloInfo[p]._bothHigh[d2] += enlargement;
					_haloInfo[p]._leavingHigh[d2] += enlargement;
					_haloInfo[p]._copiesHigh[d2] += enlargement;
					_haloInfo[p]._enlarged[d2][1] = true;
				}
			}
		}
	}

	//! Combines current CommunicationPartner with the given partner
	//! @param partner which to add to the current CommunicationPartner
	void add(CommunicationPartner partner);

	size_t getDynamicSize();

	void print(std::ostream& stream) const;
private:
	enum HaloOrLeavingCorrection{
		HALO,
		LEAVING,
		NONE,
		FORCES // necessary?
	};
	void collectMoleculesInRegion(ParticleContainer* moleculeContainer, const double lowCorner[3],
			const double highCorner[3], const double shift[3], bool removeFromContainer,
			HaloOrLeavingCorrection haloLeaveCorr, bool doHaloPositionCheck = true);

	int _rank;
	int _countTested;
	std::vector<PositionInfo> _haloInfo;

	// technical variables
	MPI_Request *_sendRequest, *_recvRequest;
	MPI_Status *_sendStatus, *_recvStatus;
	CommunicationBuffer _sendBuf, _recvBuf; // used to be ParticleData and 
	bool _msgSent, _countReceived, _msgReceived;

	void collectLeavingMoleculesFromInvalidParticles(std::vector<Molecule>& invalidParticles, double lowCorner [3], double highCorner [3], double shift [3]);

	friend class NeighborAcquirerTest;
};

#endif /* COMMUNICATIONPARTNER_H_ */
