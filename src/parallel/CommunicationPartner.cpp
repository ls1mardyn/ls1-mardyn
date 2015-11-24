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
		double bLo[3], double bHi[3], double sh, int signDir) {
	_rank = r;

	for (int d = 0; d < 3; ++d) {
		_haloLow[d] = hLo[d];
		_haloHigh[d] = hHi[d];
		_boundaryLow[d] = bLo[d];
		_boundaryHigh[d] = bHi[d];
		_regionLow[d] = fmin(hLo[d], bLo[d]);
		_regionHigh[d] = fmax(hHi[d], bHi[d]);
	}

	_shift = sh;
	_signedDirection = signDir;

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
		_regionLow[d] = o._regionLow[d];
		_regionHigh[d] = o._regionHigh[d];
		_haloLow[d] = o._haloLow[d];
		_haloHigh[d] = o._haloHigh[d];
		_boundaryLow[d] = o._boundaryLow[d];
		_boundaryHigh[d] = o._boundaryHigh[d];
	}

	_shift = o._shift;
	_signedDirection = o._signedDirection;

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

void CommunicationPartner::initCommunication(unsigned short d,
		ParticleContainer* moleculeContainer, const MPI_Comm& comm,
		const MPI_Datatype& type, MessageType msgType) {
	using std::vector;
	using Log::global_log;

	global_log->debug() << _rank << std::endl;

	vector<Molecule*> particles;
	vector<Molecule*> part2;

	switch(msgType) {
	case LEAVING_AND_HALO_COPIES:
//			moleculeContainer->getRegion(_regionLow, _regionHigh, particles);
		global_log->debug() << "sending halo and boundary particles together" << std::endl;
		moleculeContainer->getHaloParticlesDirection(_signedDirection, particles, false);
		moleculeContainer->getBoundaryParticlesDirection(_signedDirection, particles);
		break;
	case LEAVING_ONLY:
//			moleculeContainer->getRegion(_haloLow, _haloHigh, particles, true);
		global_log->debug() << "sending halo particles only" << std::endl;
		moleculeContainer->getHaloParticlesDirection(_signedDirection, particles, true);
		break;
	case HALO_COPIES:
//			moleculeContainer->getRegion(_boundaryLow, _boundaryHigh, particles);
		global_log->debug() << "sending boundary particles only" << std::endl;
		moleculeContainer->getBoundaryParticlesDirection(_signedDirection, particles);
		break;
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
		_sendBuf[i].r[d] += _shift;
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

bool CommunicationPartner::testRecv(ParticleContainer* moleculeContainer) {
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
				moleculeContainer->addParticlePointer(m);
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

void CommunicationPartner::deadlockDiagnostic() {
	using Log::global_log;
	// intentionally using std::cout instead of global_log, we want the messages from all processes
	if (not _msgSent) {
		global_log->warning() << "Send request to " << _rank
				<< " not yet completed" << std::endl;
	}
	if (not _countReceived) {
		global_log->warning() << "Probe request to " << _rank
				<< " not yet completed" << std::endl;
	}
	if (not _msgReceived) {
		global_log->warning() << "Recv request to " << _rank
				<< " not yet completed" << std::endl;
	}
}
