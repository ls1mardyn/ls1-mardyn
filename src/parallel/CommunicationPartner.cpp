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

CommunicationPartner::CommunicationPartner(const int r, const double hLo[3], const double hHi[3], const double bLo[3],
		const double bHi[3], const double sh[3], const int offset[3]) {
	_rank = r;

	PositionInfo p;
	for (int d = 0; d < 3; ++d) {
		p._leavingLow[d] = hLo[d];
		p._leavingHigh[d] = hHi[d];
		p._copiesLow[d] = bLo[d];
		p._copiesHigh[d] = bHi[d];
		p._bothLow[d] = fmin(hLo[d], bLo[d]);
		p._bothHigh[d] = fmax(hHi[d], bHi[d]);
		p._shift[d] = sh[d];
		p._offset[d] = offset[d];
	}
	_haloInfo.push_back(p);

	// some values, to silence the warnings:
	_sendRequest = new MPI_Request;
	_recvRequest = new MPI_Request;
	_sendStatus = new MPI_Status;
	_recvStatus = new MPI_Status;
	_msgSent = _countReceived = _msgReceived = false;
}

CommunicationPartner::CommunicationPartner(const int r) {
	_rank = r;
	PositionInfo p;
	for (int d = 0; d < 3; ++d) {
		p._leavingLow[d] = 0.;
		p._leavingHigh[d] = 0.;
		p._copiesLow[d] = 0.;
		p._copiesHigh[d] = 0.;
		p._bothLow[d] = 0.;
		p._bothHigh[d] = 0.;
		p._shift[d] = 0.;
		p._offset[d] = 0;
	}
	_haloInfo.push_back(p);

	// some values, to silence the warnings:
	_sendRequest = new MPI_Request;
	_recvRequest = new MPI_Request;
	_sendStatus = new MPI_Status;
	_recvStatus = new MPI_Status;
	_msgSent = _countReceived = _msgReceived = false;
}

CommunicationPartner::CommunicationPartner(const int r, const double leavingLo[3], const double leavingHigh[3]) {
	_rank = r;
	PositionInfo p;
	for (int d = 0; d < 3; ++d) {
		p._leavingLow[d] = leavingLo[d];
		p._leavingHigh[d] = leavingHigh[d];
		p._copiesLow[d] = 0.;
		p._copiesHigh[d] = 0.;
		p._bothLow[d] = 0.;
		p._bothHigh[d] = 0.;
		p._shift[d] = 0.;
		p._offset[d] = 0;
	}
	_haloInfo.push_back(p);
	// some values, to silence the warnings:
	_sendRequest = new MPI_Request;
	_recvRequest = new MPI_Request;
	_sendStatus = new MPI_Status;
	_recvStatus = new MPI_Status;
	_msgSent = _countReceived = _msgReceived = false;
}

CommunicationPartner::CommunicationPartner(const CommunicationPartner& o) {
	_rank = o._rank;

	_haloInfo = o._haloInfo;

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

void CommunicationPartner::initSend(ParticleContainer* moleculeContainer, const MPI_Comm& comm,
		const MPI_Datatype& type, MessageType msgType, bool removeFromContainer) {

	global_log->debug() << _rank << std::endl;

	std::vector<Molecule*> particles;
	std::vector<size_t> endings(_haloInfo.size() + 1, 0);  // stores last positions of the particles for each haloInfo

	switch (msgType) {
	case LEAVING_AND_HALO_COPIES: {
		for (unsigned int p = 0; p < _haloInfo.size(); p++) {
			moleculeContainer->getRegionSimple(_haloInfo[p]._bothLow, _haloInfo[p]._bothHigh, particles);
			endings[p+1] = particles.size();
		}
		global_log->debug() << "sending halo and boundary particles together" << std::endl;
		break;
	}
	case LEAVING_ONLY: {
		for (unsigned int p = 0; p < _haloInfo.size(); p++) {
			moleculeContainer->getRegionSimple(_haloInfo[p]._leavingLow, _haloInfo[p]._leavingHigh, particles,
					removeFromContainer);
			endings[p+1] = particles.size();
		}
		global_log->debug() << "sending halo particles only" << std::endl;
		break;
	}
	case HALO_COPIES: {
		for (unsigned int p = 0; p < _haloInfo.size(); p++) {
			moleculeContainer->getRegionSimple(_haloInfo[p]._copiesLow, _haloInfo[p]._copiesHigh, particles);
			endings[p+1] = particles.size();
		}
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

	for (unsigned int p = 0; p < _haloInfo.size(); p++) {
		for (unsigned int i = endings[p]; i < endings[p+1]; ++i) {
			ParticleData::MoleculeToParticleData(_sendBuf[i], *(particles[i]));
			// add offsets for particles transfered over the periodic boundary
			for (int d = 0; d < 3; ++d) {
				_sendBuf[i].r[d] += _haloInfo[p]._shift[d];
			}
		}
	}
	if (removeFromContainer) {
		for (int i = 0; i < n; ++i) {
			delete particles[i];
		}
	}

	MPI_CHECK(MPI_Isend(&(_sendBuf[0]), (int ) _sendBuf.size(), type, _rank, 99, comm, _sendRequest));
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
			MPI_CHECK(MPI_Irecv(&(_recvBuf[0]), numrecv, type, _rank, 99, comm, _recvRequest));
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
				Molecule m;
				ParticleData::ParticleDataToMolecule(_recvBuf[i], m);
				const bool inBoxCheckedAlready = false;
				moleculeContainer->addParticle(m, inBoxCheckedAlready, removeRecvDuplicates);
#ifndef NDEBUG
				buf << m.id() << " ";
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
	MPI_CHECK(MPI_Irecv(&(_recvBuf[0]), numParticles, type, _rank, 99, comm, _recvRequest));
}

void CommunicationPartner::deadlockDiagnosticSendRecv() {
	using Log::global_log;

	deadlockDiagnosticSend();

	if (not _countReceived) {
		global_log->warning() << "Probe request to " << _rank << " not yet completed" << std::endl;
	}

	deadlockDiagnosticRecv();
}

void CommunicationPartner::deadlockDiagnosticSend() {
	// intentionally using std::cout instead of global_log, we want the messages from all processes
	if (not _msgSent) {
		Log::global_log->warning() << "Send request to " << _rank << " not yet completed" << std::endl;
	}
}

void CommunicationPartner::deadlockDiagnosticRecv() {
	if (not _msgReceived) {
		Log::global_log->warning() << "Recv request to " << _rank << " not yet completed" << std::endl;
	}
}

void CommunicationPartner::add(CommunicationPartner partner) {
	assert(partner._rank == _rank);
	_haloInfo.push_back(partner._haloInfo[0]);
}
