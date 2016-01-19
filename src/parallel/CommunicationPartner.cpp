/*
 * CommunicationPartner.cpp
 *
 *  Created on: Nov 23, 2015
 *      Author: tchipevn
 */

#include "CommunicationPartner.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"
#include <cmath>
#include <sstream>

CommunicationPartner::CommunicationPartner(int r, double hLo[3], double hHi[3],
		double bLo[3], double bHi[3], double sh[3]) {
	_rank = r;

	for (int d = 0; d < 3; ++d) {
		_leavingLow[d] = hLo[d];
		_leavingHigh[d] = hHi[d];
		_copiesLow[d] = bLo[d];
		_copiesHigh[d] = bHi[d];
		_bothLow[d] = fmin(hLo[d], bLo[d]);
		_bothHigh[d] = fmax(hHi[d], bHi[d]);
		_shift[d] = sh[d];
	}

	// some values, to silence the warnings:
	_sendRequest = new MPI_Request;
	_recvRequest = new MPI_Request;
	_sendStatus = new MPI_Status;
	_recvStatus = new MPI_Status;
	_msgSent = _countReceived = _msgReceived = false;
}

CommunicationPartner::CommunicationPartner(int r) {
	_rank = r;

	for (int d = 0; d < 3; ++d) {
		_leavingLow[d] = 0.;
		_leavingHigh[d] = 0.;
		_copiesLow[d] = 0.;
		_copiesHigh[d] = 0.;
		_bothLow[d] = 0.;
		_bothHigh[d] = 0.;
		_shift[d] = 0.;
	}

	// some values, to silence the warnings:
	_sendRequest = new MPI_Request;
	_recvRequest = new MPI_Request;
	_sendStatus = new MPI_Status;
	_recvStatus = new MPI_Status;
	_msgSent = _countReceived = _msgReceived = false;
}

CommunicationPartner::CommunicationPartner(int r, double leavingLo[3], double leavingHigh[3]) {
	_rank = r;

	for (int d = 0; d < 3; ++d) {
		_leavingLow[d] = leavingLo[d];
		_leavingHigh[d] = leavingHigh[d];
		_copiesLow[d] = 0.;
		_copiesHigh[d] = 0.;
		_bothLow[d] = 0.;
		_bothHigh[d] = 0.;
		_shift[d] = 0.;
	}

	// some values, to silence the warnings:
	_sendRequest = new MPI_Request;
	_recvRequest = new MPI_Request;
	_sendStatus = new MPI_Status;
	_recvStatus = new MPI_Status;
	_msgSent = _countReceived = _msgReceived = false;
}

CommunicationPartner::CommunicationPartner(const CommunicationPartner& o) {
	_rank = o._rank;

	for(int d = 0; d < 3; ++d) {
		_bothLow[d] = o._bothLow[d];
		_bothHigh[d] = o._bothHigh[d];
		_leavingLow[d] = o._leavingLow[d];
		_leavingHigh[d] = o._leavingHigh[d];
		_copiesLow[d] = o._copiesLow[d];
		_copiesHigh[d] = o._copiesHigh[d];
		_shift[d] = o._shift[d];
	}

	// some values, to silence the warnings:
	_sendRequest = new MPI_Request;
	_recvRequest = new MPI_Request;
	_sendStatus = new MPI_Status;
	_recvStatus = new MPI_Status;
	_msgSent = _countReceived = _msgReceived = false;
}

CommunicationPartner::~CommunicationPartner() {
	delete _sendRequest;
	delete _recvRequest;
	delete _sendStatus;
	delete _recvStatus;
}

void CommunicationPartner::initSend(
		ParticleContainer* moleculeContainer, const MPI_Comm& comm,
		const MPI_Datatype& type, MessageType msgType) {
	using std::vector;
	using Log::global_log;

	global_log->debug() << _rank << std::endl;

	vector<Molecule*> particles;

	switch(msgType) {
	case LEAVING_AND_HALO_COPIES: {
		moleculeContainer->getRegionSimple(_bothLow, _bothHigh, particles);
		global_log->debug() << "sending halo and boundary particles together" << std::endl;
		break;
	}
	case LEAVING_ONLY: {
		moleculeContainer->getRegionSimple(_leavingLow, _leavingHigh, particles, true);
		global_log->debug() << "sending halo particles only" << std::endl;
		break;
	}
	case HALO_COPIES: {
		moleculeContainer->getRegionSimple(_copiesLow, _copiesHigh, particles);
		global_log->debug() << "sending boundary particles only" << std::endl;
		break;
	}
	}

	const int n = particles.size();

#ifndef NDEBUG
	global_log->debug() << "Buffer contains " << n << " particles with IDs " << std::endl;
	std::ostringstream buf;
	for (int i = 0; i < n; ++i) {
		buf << particles[i]->id() << " ";
	}
	global_log->debug() << buf.str() << std::endl;
#endif

	// initialize send buffer
	_sendBuf.resize(n);

	for (int i = 0; i < n; ++i) {
		ParticleData::MoleculeToParticleData(_sendBuf[i], *(particles[i]));
		// add offsets for particles transfered over the periodic boundary
		for (int d = 0; d < 3; ++d) {
			_sendBuf[i].r[d] += _shift[d];
		}
	}

	MPI_CHECK(
			MPI_Isend(&(_sendBuf[0]), (int) _sendBuf.size(), type, _rank,
					99, comm, _sendRequest));
	_msgSent = _countReceived = _msgReceived = false;
}

bool CommunicationPartner::testSend() {
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

bool CommunicationPartner::iprobeCount(const MPI_Comm& comm, const MPI_Datatype& type) {
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

bool CommunicationPartner::testRecv(ParticleContainer* moleculeContainer, bool removeRecvDuplicates) {
	using Log::global_log;
	if (_countReceived and not _msgReceived) {
		int flag = 0;
		MPI_CHECK(MPI_Test(_recvRequest, &flag, _recvStatus));
		if (flag == 1) {
			_msgReceived = true;
			int numrecv = _recvBuf.size();

#ifndef NDEBUG
			global_log->debug() << "Receiving particles from" << _rank << std::endl;
			global_log->debug() << "Buffer contains " << numrecv << " particles with IDs " << std::endl;
			std::ostringstream buf;
#endif

			for (int i = 0; i < numrecv; i++) {
				Molecule *m;
				ParticleData::ParticleDataToMolecule(_recvBuf[i], &m);
				const bool inBoxCheckedAlready = false;
				moleculeContainer->addParticlePointer(m, inBoxCheckedAlready, removeRecvDuplicates);
#ifndef NDEBUG
				buf << m->id() << " ";
#endif
			}
#ifndef NDEBUG
			global_log->debug() << buf.str() << std::endl;
#endif
			_recvBuf.clear();
		}
	}
	return _msgReceived;
}

void CommunicationPartner::initRecv(int numParticles, const MPI_Comm& comm, const MPI_Datatype& type) {
	_countReceived = true;
	_recvBuf.resize(numParticles);
	MPI_CHECK(
			MPI_Irecv(&(_recvBuf[0]), numParticles, type, _rank, 99,
					comm, _recvRequest));
}

void CommunicationPartner::deadlockDiagnosticSendRecv() {
	using Log::global_log;

	deadlockDiagnosticSend();

	if (not _countReceived) {
		global_log->warning() << "Probe request to " << _rank
				<< " not yet completed" << std::endl;
	}

	deadlockDiagnosticRecv();
}

void CommunicationPartner::deadlockDiagnosticSend() {
	// intentionally using std::cout instead of global_log, we want the messages from all processes
	if (not _msgSent) {
		Log::global_log->warning() << "Send request to " << _rank
				<< " not yet completed" << std::endl;
	}
}

void CommunicationPartner::deadlockDiagnosticRecv() {
	if (not _msgReceived) {
		Log::global_log->warning() << "Recv request to " << _rank
				<< " not yet completed" << std::endl;
	}
}
