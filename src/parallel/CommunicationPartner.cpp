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
#include "ParticleData.h"
#include "parallel/DomainDecompBase.h"
#include "Domain.h"

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

template<typename BufferType>
void CommunicationPartner::initSend(ParticleContainer* moleculeContainer, const MPI_Comm& comm, const MPI_Datatype& type,
		MessageType msgType, bool removeFromContainer) {

	auto& sendBuf = getSendBuf<BufferType>();

	global_log->debug() << _rank << std::endl;
	_sendBuf.clear();

	const unsigned int numHaloInfo = _haloInfo.size();
	switch (msgType){
		case LEAVING_AND_HALO_COPIES: {
			//"This message type requires a ParticleData buffer."
			mardyn_assert(!isForceData<BufferType>());
			global_log->debug() << "sending halo and boundary particles together" << std::endl;
<<<<<<< .working
			for(unsigned int p = 0; p < numHaloInfo; p++){ // -
				collectMoleculesInRegion<BufferType>(moleculeContainer, _haloInfo[p]._bothLow, _haloInfo[p]._bothHigh, _haloInfo[p]._shift); // -
||||||| .merge-left.r4919
			for(unsigned int p = 0; p < numHaloInfo; p++){ // -
				collectMoleculesInRegion(moleculeContainer, _haloInfo[p]._bothLow, _haloInfo[p]._bothHigh, _haloInfo[p]._shift); // -
=======
			// first leaving particles:
			for (unsigned int p = 0; p < numHaloInfo; p++) { // +
				collectMoleculesInRegion(moleculeContainer, _haloInfo[p]._leavingLow, _haloInfo[p]._leavingHigh, // +
						_haloInfo[p]._shift, removeFromContainer, LEAVING); // +
>>>>>>> .merge-right.r5797
			}

			// then halo particles/copies:
			for (unsigned int p = 0; p < numHaloInfo; p++) {
				collectMoleculesInRegion(moleculeContainer, _haloInfo[p]._copiesLow, _haloInfo[p]._copiesHigh,
						_haloInfo[p]._shift, false, HALO);
			}
			break;
		}
		case LEAVING_ONLY: {
			//"This message type requires a ParticleData buffer."
			mardyn_assert(!isForceData<BufferType>());
			global_log->debug() << "sending leaving particles only" << std::endl;
			for(unsigned int p = 0; p < numHaloInfo; p++){
<<<<<<< .working
				collectMoleculesInRegion<BufferType>(moleculeContainer, _haloInfo[p]._leavingLow, _haloInfo[p]._leavingHigh, _haloInfo[p]._shift, removeFromContainer); // -
||||||| .merge-left.r4919
				collectMoleculesInRegion(moleculeContainer, _haloInfo[p]._leavingLow, _haloInfo[p]._leavingHigh, _haloInfo[p]._shift, removeFromContainer); // -
=======
				collectMoleculesInRegion(moleculeContainer, _haloInfo[p]._leavingLow, _haloInfo[p]._leavingHigh, // +
						_haloInfo[p]._shift, removeFromContainer, LEAVING); // +
>>>>>>> .merge-right.r5797
			}
			break;
		}
		case HALO_COPIES: {
			//"This message type requires a ParticleData buffer."
			mardyn_assert(!isForceData<BufferType>());
			global_log->debug() << "sending halo particles only" << std::endl;
			for(unsigned int p = 0; p < numHaloInfo; p++){
<<<<<<< .working
				collectMoleculesInRegion<BufferType>(moleculeContainer, _haloInfo[p]._copiesLow, _haloInfo[p]._copiesHigh, _haloInfo[p]._shift); // -
||||||| .merge-left.r4919
				collectMoleculesInRegion(moleculeContainer, _haloInfo[p]._copiesLow, _haloInfo[p]._copiesHigh, _haloInfo[p]._shift); // -
=======
				collectMoleculesInRegion(moleculeContainer, _haloInfo[p]._copiesLow, _haloInfo[p]._copiesHigh, // +
						_haloInfo[p]._shift, false, HALO); // +
>>>>>>> .merge-right.r5797
			}
			break;
		}
		case FORCES: {
			// "A force message type requires a ParticleForceData buffer."
			mardyn_assert(isForceData<BufferType>());
			global_log->debug() << "sending forces" << std::endl;
			for(unsigned int p = 0; p < numHaloInfo; p++){
				collectMoleculesInRegion<ParticleForceData>(moleculeContainer, _haloInfo[p]._leavingLow, _haloInfo[p]._leavingHigh, _haloInfo[p]._shift);
			}
			break;
		}
	}

	#ifndef NDEBUG
<<<<<<< .working
		const int n = sendBuf.size();
		global_log->debug() << "Buffer contains " << n << " particles with IDs " << std::endl;
		std::ostringstream buf;
		for (int i = 0; i < n; ++i) {
			buf << sendBuf[i].id << " "; // How do we know, that his is not a CommunicationBuffer?
||||||| .merge-left.r4919
		const int n = _sendBuf.size();
		global_log->debug() << "Buffer contains " << n << " particles with IDs " << std::endl;
		std::ostringstream buf;
		for (int i = 0; i < n; ++i) {
			buf << _sendBuf[i].id << " ";
=======
		const int numLeaving = _sendBuf.getNumLeaving(); // new
		const int numHalo = _sendBuf.getNumHalo(); // new
		global_log->debug() << "Buffer contains " << numLeaving << " leaving particles with IDs " << std::endl; // same
		std::ostringstream buf1; // buf now buf1
		for (int i = 0; i < numLeaving; ++i) {
			Molecule m;
			_sendBuf.readLeavingMolecule(i, m); // writes into m
			buf1 << m.id() << " "; // why is this a molecule?
>>>>>>> .merge-right.r5797
		}
		global_log->debug() << buf1.str() << std::endl;

		global_log->debug() << "and " << numHalo << " halo particles with IDs " << std::endl;
		std::ostringstream buf2;
		for (int i = 0; i < numHalo; ++i) {
			Molecule m;
			_sendBuf.readHaloMolecule(i, m);
			buf2 << m.id() << " ";
		}
		global_log->debug() << buf2.str() << std::endl;


	#endif

<<<<<<< .working
	MPI_CHECK(MPI_Isend(sendBuf.data(), (int ) sendBuf.size(), type, _rank, 99, comm, _sendRequest)); // -
||||||| .merge-left.r4919
	MPI_CHECK(MPI_Isend(_sendBuf.data(), (int ) _sendBuf.size(), type, _rank, 99, comm, _sendRequest)); // -
=======
	MPI_CHECK(MPI_Isend(_sendBuf.getDataForSending(), (int ) _sendBuf.getNumElementsForSending(), _sendBuf.getMPIDataType(), _rank, 99, comm, _sendRequest)); // CommunicationBuffer changed
>>>>>>> .merge-right.r5797
	_msgSent = _countReceived = _msgReceived = false;
}

template<typename BufferType>
bool CommunicationPartner::testSend() {

	auto& sendBuf = getSendBuf<BufferType>();

	if (not _msgSent) {
		int flag = 0;
		MPI_CHECK(MPI_Test(_sendRequest, &flag, _sendStatus));
		if (flag == 1) {
			_msgSent = true;
			sendBuf.clear();
		}
	}
	return _msgSent;
}

<<<<<<< .working
template<typename BufferType>
bool CommunicationPartner::iprobeCount(const MPI_Comm& comm, const MPI_Datatype& type) {
||||||| .merge-left.r4919
bool CommunicationPartner::iprobeCount(const MPI_Comm& comm, const MPI_Datatype& type) {
=======
bool CommunicationPartner::iprobeCount(const MPI_Comm& comm, const MPI_Datatype& /*type*/) {
>>>>>>> .merge-right.r5797

	auto& recvBuf = getRecvBuf<BufferType>();

	if (not _countReceived) {
		int flag = 0;
		MPI_CHECK(MPI_Iprobe(_rank, 99, comm, &flag, _recvStatus));
		if (flag == true) {
			_countReceived = true;
			_countTested = 0;
			int numrecv;
			MPI_CHECK(MPI_Get_count(_recvStatus, _sendBuf.getMPIDataType(), &numrecv));
                        #ifndef NDEBUG
                                global_log->debug() << "Received byteCount from " << _rank << std::endl;
                                global_log->debug() << "Preparing to receive " << numrecv << " bytes." << std::endl;
                        #endif
<<<<<<< .working
			recvBuf.resize(numrecv); // -
			MPI_CHECK(MPI_Irecv(recvBuf.data(), numrecv, type, _rank, 99, comm, _recvRequest)); // -
||||||| .merge-left.r4919
			_recvBuf.resize(numrecv); // -
			MPI_CHECK(MPI_Irecv(_recvBuf.data(), numrecv, type, _rank, 99, comm, _recvRequest)); // -
=======
			_recvBuf.resizeForRawBytes(numrecv); // +
			MPI_CHECK(MPI_Irecv(_recvBuf.getDataForSending(), numrecv, _sendBuf.getMPIDataType(), _rank, 99, comm, _recvRequest)); // +
>>>>>>> .merge-right.r5797
		}
	}
	return _countReceived;
}

template<typename BufferType>
bool CommunicationPartner::testRecv(ParticleContainer* moleculeContainer, bool removeRecvDuplicates) {
	using Log::global_log;

	auto& recvBuf = getRecvBuf<BufferType>(); // is this necessary?

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
<<<<<<< .working
			int numrecv = recvBuf.size(); // again the buffer thing
||||||| .merge-left.r4919
			int numrecv = _recvBuf.size(); // prob. should use this.
=======
>>>>>>> .merge-right.r5797

			unsigned long numHalo, numLeaving;
			_recvBuf.resizeForReceivingMolecules(numLeaving, numHalo);

			#ifndef NDEBUG
				global_log->debug() << "Receiving particles from " << _rank << std::endl;
<<<<<<< .working
				global_log->debug() << "Buffer contains " << numrecv << " particles with IDs " << std::endl; // -
||||||| .merge-left.r4919
				global_log->debug() << "Buffer contains " << numrecv << " particles with IDs " << std::endl; // -
				std::ostringstream buf; // -
=======
				// the new stuff
				global_log->debug() << "Buffer contains " << numLeaving << " leaving particles with IDs " << std::endl;
				std::ostringstream buf1;
				for (unsigned long i = 0; i < numLeaving; ++i) {
					Molecule m;
					_recvBuf.readLeavingMolecule(i, m);
					buf1 << m.id() << " ";
				}
				global_log->debug() << buf1.str() << std::endl;

				global_log->debug() << "and " << numHalo << " halo particles with IDs " << std::endl;
				std::ostringstream buf2;
				for (unsigned long i = 0; i < numHalo; ++i) {
					Molecule m;
					_recvBuf.readHaloMolecule(i, m);
					buf2 << m.id() << " ";
				}
				global_log->debug() << buf2.str() << std::endl;
>>>>>>> .merge-right.r5797
			#endif

<<<<<<< .working
			//Code for different buffer types
			testRecvHandle<BufferType>(moleculeContainer, removeRecvDuplicates, numrecv); // don't know
||||||| .merge-left.r4919
			global_simulation->startTimer("COMMUNICATION_PARTNER_TEST_RECV");
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
			global_simulation->stopTimer("COMMUNICATION_PARTNER_TEST_RECV");
=======
			global_simulation->timers()->start("COMMUNICATION_PARTNER_TEST_RECV");
			unsigned long totalNumMols = numLeaving + numHalo;
			static std::vector<Molecule> mols;
			mols.resize(totalNumMols);

			/*#if defined(_OPENMP) and not defined (ADVANCED_OVERLAPPING)
			#pragma omp parallel for schedule(static)
			#endif*/
			for (unsigned long i = 0; i < totalNumMols; i++) {
				Molecule m;
				if (i < numLeaving) {
					// leaving
					_recvBuf.readLeavingMolecule(i, m);
				} else {
					// halo
					_recvBuf.readHaloMolecule(i - numLeaving, m);
				}
				mols[i] = m;
			}
			global_simulation->timers()->stop("COMMUNICATION_PARTNER_TEST_RECV");
>>>>>>> .merge-right.r5797

<<<<<<< .working
||||||| .merge-left.r4919
			#ifndef NDEBUG
				for (int i = 0; i < numrecv; i++) {
					buf << mols[i].id() << " ";
				}
				global_log->debug() << buf.str() << std::endl;
			#endif

			global_simulation->startTimer("COMMUNICATION_PARTNER_TEST_RECV");
			moleculeContainer->addParticles(mols, removeRecvDuplicates);
			mols.clear();
			_recvBuf.clear();
			global_simulation->stopTimer("COMMUNICATION_PARTNER_TEST_RECV");

======= // why did the sauermann branch delete this? 
			global_simulation->timers()->start("COMMUNICATION_PARTNER_TEST_RECV");
			moleculeContainer->addParticles(mols, removeRecvDuplicates);
			mols.clear();
			_recvBuf.clear();
			global_simulation->timers()->stop("COMMUNICATION_PARTNER_TEST_RECV");

>>>>>>> .merge-right.r5797
		} else {
			++_countTested;
		}
	}
	return _msgReceived;
}


template<typename BufferType>
void CommunicationPartner::initRecv(int numParticles, const MPI_Comm& comm, const MPI_Datatype& type){

	auto& recvBuf = getRecvBuf<BufferType>();

	// one single call from KDDecomposition::migrate particles.
	// So all molecules, which arrive area leaving molecules.
	_countReceived = true;
<<<<<<< .working
	recvBuf.resize(numParticles); // -
	MPI_CHECK(MPI_Irecv(recvBuf.data(), numParticles, type, _rank, 99, comm, _recvRequest)); // -
||||||| .merge-left.r4919
	_recvBuf.resize(numParticles); // -
	MPI_CHECK(MPI_Irecv(_recvBuf.data(), numParticles, type, _rank, 99, comm, _recvRequest)); // -
=======

	// hackaround - resizeForAppendingLeavingMolecules is intended for the send-buffer, not the recv one.
	_recvBuf.resizeForAppendingLeavingMolecules(numParticles); // +

	MPI_CHECK(MPI_Irecv(_recvBuf.getDataForSending(), _recvBuf.getNumElementsForSending(), _sendBuf.getMPIDataType(), _rank, 99, comm, _recvRequest)); // +
>>>>>>> .merge-right.r5797
}


template <>
void CommunicationPartner::testRecvHandle<ParticleData>(ParticleContainer* moleculeContainer, bool removeRecvDuplicates, int numrecv) {

	auto& recvBuf = getRecvBuf<ParticleData>();

	global_simulation->startTimer("COMMUNICATION_PARTNER_TEST_RECV");
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
	global_simulation->stopTimer("COMMUNICATION_PARTNER_TEST_RECV");

	#ifndef NDEBUG
	std::ostringstream buf;
		for (int i = 0; i < numrecv; i++) {
			buf << mols[i].id() << " ";
		}
		global_log->debug() << buf.str() << std::endl;
	#endif

	global_simulation->startTimer("COMMUNICATION_PARTNER_TEST_RECV");
	moleculeContainer->addParticles(mols, removeRecvDuplicates);
	mols.clear();
	recvBuf.clear();
	global_simulation->stopTimer("COMMUNICATION_PARTNER_TEST_RECV");

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


/* -------------------------------------------------------------------------- */
/* - COLLECT MOLECULES IN REGION -------------------------------------------- */
/* -------------------------------------------------------------------------- */

<<<<<<< .working
template<typename BufferType> // -
void CommunicationPartner::collectMoleculesInRegion(ParticleContainer* moleculeContainer, const double lowCorner[3], const double highCorner[3], const double shift[3], // -
		const bool removeFromContainer){ // -
||||||| .merge-left.r4919 
void CommunicationPartner::collectMoleculesInRegion(ParticleContainer* moleculeContainer, const double lowCorner[3], // -
		const double highCorner[3], const double shift[3], const bool removeFromContainer) { // -
=======
// new declaration - use this
void CommunicationPartner::collectMoleculesInRegion(ParticleContainer* moleculeContainer, const double lowCorner[3], // +
		const double highCorner[3], const double shift[3], const bool removeFromContainer, HaloOrLeavingCorrection haloLeaveCorr) { // +
>>>>>>> .merge-right.r5797
	using std::vector;
<<<<<<< .working

	auto& sendBuf = getSendBuf<BufferType>(); // WHAT DOES THIS DO?

	global_simulation->startTimer("COMMUNICATION_PARTNER_INIT_SEND");
	int prevNumMols = sendBuf.size(); // WHY IS THIS GONE? - New methods in CommunicationBuffer
||||||| .merge-left.r4919
	global_simulation->startTimer("COMMUNICATION_PARTNER_INIT_SEND");
	int prevNumMols = _sendBuf.size();
=======
	global_simulation->timers()->start("COMMUNICATION_PARTNER_INIT_SEND"); // new timer start function
>>>>>>> .merge-right.r5797
	vector<vector<Molecule>> threadData;
	vector<int> prefixArray;

	// compute how many molecules are already in of this type:
	unsigned long numMolsAlreadyIn = 0;
	if (haloLeaveCorr == LEAVING) {
		numMolsAlreadyIn = _sendBuf.getNumLeaving();
	} else if (haloLeaveCorr == HALO) {
		numMolsAlreadyIn = _sendBuf.getNumHalo();
	}

	#if defined (_OPENMP)
	#pragma omp parallel shared(threadData, numMolsAlreadyIn)
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
			int totalNumMolsAppended = 0;
			//build the prefix array
			prefixArray[0] = 0;
			for(int i = 1; i <= numThreads; i++){
				prefixArray[i] += prefixArray[i - 1];
				totalNumMolsAppended += threadData[i - 1].size();
			}

			//resize the send buffer
<<<<<<< .working
			sendBuf.resize(prevNumMols + totalNumMols); // -
||||||| .merge-left.r4919
			_sendBuf.resize(prevNumMols + totalNumMols); // -
=======

			// New parameter is used - new methods in CommunicationBuffer - resize does not exist anymore
			if (haloLeaveCorr == LEAVING) { // +
				_sendBuf.resizeForAppendingLeavingMolecules(totalNumMolsAppended); // +
			} else if (haloLeaveCorr == HALO) { // +
				_sendBuf.resizeForAppendingHaloMolecules(totalNumMolsAppended); // +
			}
>>>>>>> .merge-right.r5797
		}

		#if defined (_OPENMP)
		#pragma omp barrier
		#endif


		Domain* domain = global_simulation->getDomain();

		//reduce the molecules in the send buffer and also apply the shift
		int myThreadMolecules = prefixArray[threadNum + 1] - prefixArray[threadNum];
		for(int i = 0; i < myThreadMolecules; i++){
<<<<<<< .working
			BufferType m; // -
			BufferType::MoleculeToParticleData(m, threadData[threadNum][i]); // -
			m.r[0] += shift[0]; // -
			m.r[1] += shift[1]; // -
			m.r[2] += shift[2]; // -

			sendBuf[prevNumMols + prefixArray[threadNum] + i] = m; // -
||||||| .merge-left.r4919
			ParticleData m; // -
			ParticleData::MoleculeToParticleData(m, threadData[threadNum][i]); // -
			m.r[0] += shift[0]; // -
			m.r[1] += shift[1]; // -
			m.r[2] += shift[2]; // -
			_sendBuf[prevNumMols + prefixArray[threadNum] + i] = m; // -
=======
			// accept this block
			Molecule mCopy = threadData[threadNum][i];
			mCopy.move(0, shift[0]);
			mCopy.move(1, shift[1]);
			mCopy.move(2, shift[2]);
			for (int dim = 0; dim < 3; dim++) {
				if (haloLeaveCorr == HALO) {
					// checks if the molecule has been shifted to inside the domain due to rounding errors.
					if (shift[dim] < 0.) { // if the shift was negative, it is now in the lower part of the domain -> min
						if (mCopy.r(dim) >= 0.) {  // in the lower part it was wrongly shifted
							//std::cout << std::endl << "shifting: molecule" << m.id << std::endl;
							vcp_real_calc r = 0;
							mCopy.setr(dim, std::nexttoward(r, r - 1.f)); // ensures that r is smaller than the boundingboxmin
						}
					} else if (shift[dim] > 0.) {  // shift > 0
						if (mCopy.r(dim) < domain->getGlobalLength(dim)) { // in the higher part it was wrongly shifted
							// std::nextafter: returns the next bigger value of _boundingBoxMax
							//std::cout << std::endl << "shifting: molecule" << m.id << std::endl;
							vcp_real_calc r = domain->getGlobalLength(dim);
							mCopy.setr(dim, std::nexttoward(r, r + 1.f));  // ensures that r is bigger than the boundingboxmax
						}
					}
				} else if (haloLeaveCorr == LEAVING) {
					// some additional shifting to ensure that rounding errors do not hinder the correct placement
					if (shift[dim] < 0) {  // if the shift was negative, it is now in the lower part of the domain -> min
						if (mCopy.r(dim) < 0.) { // in the lower part it was wrongly shifted if
							mCopy.setr(dim, 0.); // ensures that r is at least the boundingboxmin
							//std::cout << std::endl << "shifting: molecule" << m.id << std::endl;
						}
					} else  if (shift[dim] > 0.) {  // shift > 0
						if (mCopy.r(dim) >= domain->getGlobalLength(dim)) { // in the lower part it was wrongly shifted if
						// std::nexttoward: returns the next bigger value of _boundingBoxMax
							vcp_real_calc r = domain->getGlobalLength(dim);
							mCopy.setr(dim, std::nexttoward(r, r - 1.f)); // ensures that r is smaller than the boundingboxmax
							//std::cout << std::endl << "shifting: molecule" << m.id << std::endl;
						}
					}
				} /* if-else HALO-LEAVING */
			} /* for-loop dim */
			if (haloLeaveCorr == LEAVING) {
				_sendBuf.addLeavingMolecule(numMolsAlreadyIn + prefixArray[threadNum] + i, mCopy);
			} else if (haloLeaveCorr == HALO) {
				_sendBuf.addHaloMolecule(numMolsAlreadyIn + prefixArray[threadNum] + i, mCopy);
			}
>>>>>>> .merge-right.r5797
		}
	}
	global_simulation->timers()->stop("COMMUNICATION_PARTNER_INIT_SEND");
}
<<<<<<< .working

//--------------------------------------------------------------------------------------------
//------------------------Explicit template instantiations below------------------------------
//--------------------------------------------------------------------------------------------

// what are these for?

template
void CommunicationPartner::initSend<ParticleData>(ParticleContainer* moleculeContainer, const MPI_Comm& comm, const MPI_Datatype& type,
		MessageType msgType, bool removeFromContainer);

template
void CommunicationPartner::initSend<ParticleForceData>(ParticleContainer* moleculeContainer, const MPI_Comm& comm, const MPI_Datatype& type,
		MessageType msgType, bool removeFromContainer);

template
bool CommunicationPartner::testSend<ParticleData>();

template
bool CommunicationPartner::testSend<ParticleForceData>();

template
bool CommunicationPartner::iprobeCount<ParticleData>(const MPI_Comm& comm, const MPI_Datatype& type);

template
bool CommunicationPartner::iprobeCount<ParticleForceData>(const MPI_Comm& comm, const MPI_Datatype& type);

template
bool CommunicationPartner::testRecv<ParticleData>(ParticleContainer* moleculeContainer, bool removeRecvDuplicates);

template
bool CommunicationPartner::testRecv<ParticleForceData>(ParticleContainer* moleculeContainer, bool removeRecvDuplicates);

template
void CommunicationPartner::initRecv<ParticleData>(int numParticles, const MPI_Comm& comm, const MPI_Datatype& type);

template
void CommunicationPartner::initRecv<ParticleForceData>(int numParticles, const MPI_Comm& comm, const MPI_Datatype& type);

//! Handle receive for ParticleData
template
void CommunicationPartner::testRecvHandle<ParticleData>(ParticleContainer* moleculeContainer, bool removeRecvDuplicates, int numrecv) ;

//! Handle receive for ParticleForceData
template
void CommunicationPartner::testRecvHandle<ParticleForceData>(ParticleContainer* moleculeContainer, bool removeRecvDuplicates, int numrecv) ;

template
void CommunicationPartner::collectMoleculesInRegion<ParticleData>(ParticleContainer* moleculeContainer, const double lowCorner[3],
		const double highCorner[3], const double shift[3], const bool removeFromContainer);

template
void CommunicationPartner::collectMoleculesInRegion<ParticleForceData>(ParticleContainer* moleculeContainer, const double lowCorner[3],
		const double highCorner[3], const double shift[3], const bool removeFromContainer);
||||||| .merge-left.r4919
=======

size_t CommunicationPartner::getDynamicSize() {
	return _sendBuf.getDynamicSize() + _recvBuf.getDynamicSize() + _haloInfo.capacity() * sizeof(PositionInfo);
}
>>>>>>> .merge-right.r5797
