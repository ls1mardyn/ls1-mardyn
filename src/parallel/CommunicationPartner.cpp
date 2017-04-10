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
#include "WrapOpenMP.h"
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

void CommunicationPartner::initSend(ParticleContainer* moleculeContainer, const MPI_Comm& comm,
		const MPI_Datatype& type, MessageType msgType, bool removeFromContainer) {

	global_log->debug() << _rank << std::endl;

	const unsigned int numHaloInfo = _haloInfo.size();
	switch (msgType){
		case LEAVING_AND_HALO_COPIES: {
			global_log->debug() << "sending halo and boundary particles together" << std::endl;
			for(unsigned int p = 0; p < numHaloInfo; p++){
				collectMoleculesInRegion(moleculeContainer, _haloInfo[p]._bothLow, _haloInfo[p]._bothHigh, _haloInfo[p]._shift);
			}
			break;
		}
		case LEAVING_ONLY: {
			global_log->debug() << "sending leaving particles only" << std::endl;
			for(unsigned int p = 0; p < numHaloInfo; p++){
				collectMoleculesInRegion(moleculeContainer, _haloInfo[p]._leavingLow, _haloInfo[p]._leavingHigh, _haloInfo[p]._shift, removeFromContainer);
			}
			break;
		}
		case HALO_COPIES: {
			global_log->debug() << "sending halo particles only" << std::endl;
			for(unsigned int p = 0; p < numHaloInfo; p++){
				collectMoleculesInRegion(moleculeContainer, _haloInfo[p]._copiesLow, _haloInfo[p]._copiesHigh, _haloInfo[p]._shift);
			}
			break;
		}
	}

	#ifndef NDEBUG
		const int n = _sendBuf.size();
		global_log->debug() << "Buffer contains " << n << " particles with IDs " << std::endl;
		std::ostringstream buf;
		for (int i = 0; i < n; ++i) {
			buf << _sendBuf[i].id << " ";
		}
		global_log->debug() << buf.str() << std::endl;
	#endif

	MPI_CHECK(MPI_Isend(_sendBuf.data(), (int ) _sendBuf.size(), type, _rank, 99, comm, _sendRequest));
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
		if (flag == true) {
			_countReceived = true;
			_countTested = 0;
			int numrecv;
			MPI_CHECK(MPI_Get_count(_recvStatus, type, &numrecv));
                        #ifndef NDEBUG
                                global_log->debug() << "Received particleCount from " << _rank << std::endl;
                                global_log->debug() << "Preparing to receive " << numrecv << " particles." << std::endl;
                        #endif
			_recvBuf.resize(numrecv);
			MPI_CHECK(MPI_Irecv(_recvBuf.data(), numrecv, type, _rank, 99, comm, _recvRequest));
		}
	}
	return _countReceived;
}
bool CommunicationPartner::testRecv(ParticleContainer* moleculeContainer, bool removeRecvDuplicates) {
	using Log::global_log;
	if (_countReceived and not _msgReceived) {
		int flag = 1;
		if (_countTested > 10) {
			// some MPI (Intel, IBM) implementations can produce deadlocks using MPI_Test without any MPI_Wait
			// this fallback just ensures, that messages get received properly.
			MPI_Wait(_recvRequest, _recvStatus);
			_countTested = 0;
			flag = 1;
		} else {
			MPI_CHECK(MPI_Test(_recvRequest, &flag, _recvStatus));
		}
		if (flag == true) {
			_msgReceived = true;
			int numrecv = _recvBuf.size();

			#ifndef NDEBUG
				global_log->debug() << "Receiving particles from " << _rank << std::endl;
				global_log->debug() << "Buffer contains " << numrecv << " particles with IDs " << std::endl;
				std::ostringstream buf;
			#endif

			static std::vector<Molecule> mols;
			mols.resize(numrecv);
			#if defined(_OPENMP)
			#pragma omp for schedule(static)
			#endif
			for (int i = 0; i < numrecv; i++) {
				Molecule m;
				ParticleData::ParticleDataToMolecule(_recvBuf[i], m);
				mols[i] = m;
			}

			#ifndef NDEBUG
				for (int i = 0; i < numrecv; i++) {
					buf << mols[i].id() << " ";
				}
				global_log->debug() << buf.str() << std::endl;
			#endif

			moleculeContainer->addParticles(mols, removeRecvDuplicates);
			mols.clear();
			_recvBuf.clear();
		} else {
			++_countTested;
		}
	}
	return _msgReceived;
}

void CommunicationPartner::initRecv(int numParticles, const MPI_Comm& comm, const MPI_Datatype& type) {
	_countReceived = true;
	_recvBuf.resize(numParticles);
	MPI_CHECK(MPI_Irecv(_recvBuf.data(), numParticles, type, _rank, 99, comm, _recvRequest));
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

void CommunicationPartner::collectMoleculesInRegion(ParticleContainer* moleculeContainer, const double lowCorner[3], const double highCorner[3], const double shift[3], const bool removeFromContainer){
	using std::vector;
	int prevNumMols = _sendBuf.size();
	vector<vector<Molecule>> threadData;
	vector<int> prefixArray;

	#if defined (_OPENMP)
	#pragma omp parallel shared(threadData)
	#endif
	{
		const int numThreads = mardyn_get_num_threads();
		const int threadNum = mardyn_get_thread_num();
		RegionParticleIterator begin = moleculeContainer->iterateRegionBegin(lowCorner, highCorner);
		RegionParticleIterator end = moleculeContainer->iterateRegionEnd();

		#if defined (_OPENMP)
		#pragma omp master
		#endif
		{
			threadData.resize(numThreads);
			prefixArray.resize(numThreads + 1);
		}

		#if defined (_OPENMP)
		#pragma omp barrier
		#endif

		for (RegionParticleIterator i = begin; i != end; ++i) {
			//traverse and gather all molecules in the cells containing part of the box specified as parameter
			//i is a pointer to a Molecule; (*i) is the Molecule
			threadData[threadNum].push_back(*i);
			if (removeFromContainer) {
				i.deleteCurrentParticle();
			}
		}

		prefixArray[threadNum + 1] = threadData[threadNum].size();

		#if defined (_OPENMP)
		#pragma omp barrier
		#endif

		//build the prefix array and resize the send buffer
		#if defined (_OPENMP)
		#pragma omp master
		#endif
		{
			int totalNumMols = 0;
			//build the prefix array
			prefixArray[0] = 0;
			for(int i = 1; i <= numThreads; i++){
				prefixArray[i] += prefixArray[i - 1];
				totalNumMols += threadData[i - 1].size();
			}

			//resize the send buffer
			_sendBuf.resize(prevNumMols + totalNumMols);
		}

		#if defined (_OPENMP)
		#pragma omp barrier
		#endif

		//reduce the molecules in the send buffer and also apply the shift
		int myThreadMolecules = prefixArray[threadNum + 1] - prefixArray[threadNum];
		for(int i = 0; i < myThreadMolecules; i++){
			ParticleData m;
			ParticleData::MoleculeToParticleData(m, threadData[threadNum][i]);
			m.r[0] += shift[0];
			m.r[1] += shift[1];
			m.r[2] += shift[2];
			_sendBuf[prevNumMols + prefixArray[threadNum] + i] = m;
		}
	}
}
