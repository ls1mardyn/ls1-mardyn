/*
 * CommunicationPartner.h
 *
 *  Created on: Nov 16, 2015
 *      Author: tchipevn
 */

#ifndef COMMUNICATIONPARTNER_H_
#define COMMUNICATIONPARTNER_H_

#include "particleContainer/ParticleContainer.h"

#include "mpi.h"

/**
 * (Bi-Directional) MPI Communication Partner.
 */
class CommunicationPartner {
public:
	CommunicationPartner(int r, double lo[3], double hi[3], double sh) {
		_rank = r;

		for(int d = 0; d < 3; ++d) {
			_regionLow[d] = lo[d];
			_regionHigh[d] = hi[d];
		}

		_shift = sh;

		// some values, to silence the warnings:
		_sendRequest = new MPI_Request;
		_recvRequest = new MPI_Request;
		_sendStatus = new MPI_Status;
		_recvStatus = new MPI_Status;
		_msgSent = _countReceived = _msgReceived = false;
	}

	CommunicationPartner(const CommunicationPartner& o) {
		_rank = o._rank;

		for(int d = 0; d < 3; ++d) {
			_regionLow[d] = o._regionLow[d];
			_regionHigh[d] = o._regionHigh[d];
		}

		_shift = o._shift;

		// some values, to silence the warnings:
		_sendRequest = new MPI_Request;
		_recvRequest = new MPI_Request;
		_sendStatus = new MPI_Status;
		_recvStatus = new MPI_Status;
		_msgSent = _countReceived = _msgReceived = false;
	}

	~CommunicationPartner() {
		delete _sendRequest;
		delete _recvRequest;
		delete _sendStatus;
		delete _recvStatus;
	}

	void initCommunication(unsigned short d,
			ParticleContainer* moleculeContainer, const MPI_Comm& comm,
			const MPI_Datatype& type) {
		using std::vector;
		vector<Molecule*> particles;
		moleculeContainer->getRegion(_regionLow, _regionHigh, particles);

		const int n = particles.size();

		// initialize send buffer
		_sendBuf.resize(n);

		for (int i = 0; i < n; ++i) {
			ParticleData::MoleculeToParticleData(_sendBuf[i], *(particles[i]));
			// add offsets for particles transfered over the periodic boundary
			_sendBuf[i].r[d] += _shift;
		}

		MPI_CHECK(
				MPI_Isend(&(_sendBuf[0]), (int) _sendBuf.size(), type, _rank,
						99, comm, _sendRequest));
		_msgSent = _countReceived = _msgReceived = false;
	}

	bool testSend() {
		if (not _msgSent) {
			int flag = 0;
			MPI_CHECK(MPI_Test(_sendRequest, &flag, _sendStatus));
			if (flag == 1) {
				_msgSent = true;
				_sendBuf.clear();
			}
		}
		return _msgSent;
	}

	bool iprobeCount(const MPI_Comm& comm, const MPI_Datatype& type) {
		if (not _countReceived) {
			int flag = 0;
			MPI_CHECK(MPI_Iprobe(_rank, 99, comm, &flag, _recvStatus));
			if (flag == 1) {
				_countReceived = true;
				int numrecv;
				MPI_CHECK(MPI_Get_count(_recvStatus, type, &numrecv));
				_recvBuf.resize(numrecv);
				MPI_CHECK(
						MPI_Irecv(&(_recvBuf[0]), numrecv, type, _rank, 99,
								comm, _recvRequest));
			}

		}
		return _countReceived;
	}

	bool testRecv(ParticleContainer* moleculeContainer) {
		if (_countReceived and not _msgReceived) {
			int flag = 0;
			MPI_CHECK(MPI_Test(_recvRequest, &flag, _recvStatus));
			if (flag == 1) {
				_msgReceived = true;
				int numrecv = _recvBuf.size();

				for (int i = 0; i < numrecv; i++) {
					Molecule *m;
					ParticleData::ParticleDataToMolecule(_recvBuf[i], &m);
					moleculeContainer->addParticlePointer(m);
				}
				_recvBuf.clear();
			}
		}
		return _msgReceived;
	}

private:
	int _rank;
	double _regionLow[3], _regionHigh[3];
	double _shift; //! for periodic boundaries

	// technical variables
	MPI_Request *_sendRequest, *_recvRequest;
	MPI_Status *_sendStatus, *_recvStatus;
	std::vector<ParticleData> _sendBuf, _recvBuf;
	bool _msgSent, _countReceived, _msgReceived;

};

#endif /* COMMUNICATIONPARTNER_H_ */
