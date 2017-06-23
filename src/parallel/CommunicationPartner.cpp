/*
 * CommunicationPartner.cpp
 *
 *  Created on: Nov 23, 2015
 *      Author: tchipevn
 */

#include "CommunicationPartner.h"
#include <cmath>
#include <sstream>
#include "Simulation.h"

CommunicationPartner::CommunicationPartner(const int r, const double hLo[3], const double hHi[3], const double bLo[3], 
		const double bHi[3], const double sh[3], const int offset[3], const bool enlarged[3][2]) {
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
		p._enlarged[d][0] = enlarged[d][0];
		p._enlarged[d][1] = enlarged[d][1];
	}
	_haloInfo.push_back(p);

	// some values, to silence the warnings:
	_sendRequest = new MPI_Request;
	_recvRequest = new MPI_Request;
	_sendStatus = new MPI_Status;
	_recvStatus = new MPI_Status;
	_msgSent = _countReceived = _msgReceived = false;
	_countTested = 0;
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
		p._enlarged[d][0] = false;
		p._enlarged[d][1] = false;
	}
	_haloInfo.push_back(p);

	// some values, to silence the warnings:
	_sendRequest = new MPI_Request;
	_recvRequest = new MPI_Request;
	_sendStatus = new MPI_Status;
	_recvStatus = new MPI_Status;
	_msgSent = _countReceived = _msgReceived = false;
	_countTested = 0;
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
		p._enlarged[d][0] = false;
		p._enlarged[d][1] = false;
	}
	_haloInfo.push_back(p);
	// some values, to silence the warnings:
	_sendRequest = new MPI_Request;
	_recvRequest = new MPI_Request;
	_sendStatus = new MPI_Status;
	_recvStatus = new MPI_Status;
	_msgSent = _countReceived = _msgReceived = false;
	_countTested = 0;
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
	_countTested = 0;
}

CommunicationPartner& CommunicationPartner::operator =(const CommunicationPartner& o){
// make sure, that the send requests are properly initialized. the delete and new operators are probably unimportant...
	if (this != &o) {
		_rank = o._rank;
		_haloInfo = o._haloInfo;
		delete _sendRequest;
		delete _recvRequest;
		delete _sendStatus;
		delete _recvStatus;
		_sendRequest = new MPI_Request;
		_recvRequest = new MPI_Request;
		_sendStatus = new MPI_Status;
		_recvStatus = new MPI_Status;
		_msgSent = _countReceived = _msgReceived = false;
		_countTested = 0;
	}
	return *this;
}

CommunicationPartner::~CommunicationPartner() {
	delete _sendRequest;
	delete _recvRequest;
	delete _sendStatus;
	delete _recvStatus;
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
	mardyn_assert(partner._rank == _rank);
	_haloInfo.push_back(partner._haloInfo[0]);
}

void CommunicationPartner::timer(const std::string name, const bool start){
	if(start)
		global_simulation->startTimer(name);
	else
		global_simulation->stopTimer(name);
}


template <>
void CommunicationPartner::testRecvHandle<ParticleData>(ParticleContainer* moleculeContainer, bool removeRecvDuplicates, int numrecv) {

	auto& recvBuf = getRecvBuf<ParticleData>();

	timer("COMMUNICATION_PARTNER_TEST_RECV", true);
	static std::vector<Molecule> mols;
	mols.resize(numrecv);
	#if defined(_OPENMP)
	#pragma omp for schedule(static)
	#endif
	for (int i = 0; i < numrecv; i++) {
		Molecule m;
		ParticleData::ParticleDataToMolecule(recvBuf[i], m);
		mols[i] = m;
	}
	timer("COMMUNICATION_PARTNER_TEST_RECV", false);

	#ifndef NDEBUG
	std::ostringstream buf;
		for (int i = 0; i < numrecv; i++) {
			buf << mols[i].id() << " ";
		}
		global_log->debug() << buf.str() << std::endl;
	#endif

	timer("COMMUNICATION_PARTNER_TEST_RECV", true);
	moleculeContainer->addParticles(mols, removeRecvDuplicates);
	mols.clear();
	recvBuf.clear();
	timer("COMMUNICATION_PARTNER_TEST_RECV", false);

}

//! Handle receive for ParticleForceData

//typename std::enable_if<std::is_same<BufferType, ParticleForceData>::value, void>::type
template<>
void CommunicationPartner::testRecvHandle<ParticleForceData>(ParticleContainer* moleculeContainer, bool removeRecvDuplicates, int numrecv) {
	auto& recvBuf = getRecvBuf<ParticleForceData>();

	#if defined(_OPENMP)
	#pragma omp for schedule(static)
	#endif
	for (int i = 0; i < numrecv; i++) {
		ParticleForceData& pData = recvBuf[i];

		Molecule* original;

		if (!moleculeContainer->getMoleculeAtPosition(pData.r, &original)) {
			// This should not happen
			global_log->error()<< "Original molecule not found!" << std::endl;
			mardyn_exit(1);
		}
		mardyn_assert(original->id() == pData.id);

		ParticleForceData::AddParticleForceDataToMolecule(pData, *original);
	}

}
