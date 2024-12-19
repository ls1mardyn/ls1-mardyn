/*
 * CommunicationPartner.cpp
 *
 *  Created on: Nov 23, 2015
 *      Author: tchipevn
 */

#include "CommunicationPartner.h"
#include <cmath>
#include <sstream>
#include "Domain.h"
#include "ForceHelper.h"
#include "ParticleData.h"
#include "Simulation.h"
#include "utils/mardyn_assert.h"
#include "WrapOpenMP.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"

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
	_isSending = _msgSent = _isReceiving = _countReceived = _msgReceived = false;
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
	_isSending = _msgSent = _isReceiving = _countReceived = _msgReceived = false;
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
	_isSending = _msgSent = _isReceiving = _countReceived = _msgReceived = false;
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
	_isSending = _msgSent = _isReceiving = _countReceived = _msgReceived = false;
	_countTested = 0;
}

CommunicationPartner& CommunicationPartner::operator =(const CommunicationPartner& o) {
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
		_isSending = _msgSent = _isReceiving = _countReceived = _msgReceived = false;
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
									const MPI_Datatype& type, MessageType msgType,
									std::vector<Molecule>& invalidParticles, bool mightUseInvalidParticles,
									bool doHaloPositionCheck, bool removeFromContainer) {
	_sendBuf.clear();

	const unsigned int numHaloInfo = _haloInfo.size();
	switch (msgType){
		case MessageType::LEAVING_AND_HALO_COPIES: {
			Log::global_log->debug() << "sending halo and boundary particles together" << std::endl;
			// first leaving particles:
			for (unsigned int p = 0; p < numHaloInfo; p++) {
				if (moleculeContainer->isInvalidParticleReturner() and mightUseInvalidParticles) {
					collectLeavingMoleculesFromInvalidParticles(invalidParticles, _haloInfo[p]._leavingLow,
																_haloInfo[p]._leavingHigh, _haloInfo[p]._shift);
				} else {
					collectMoleculesInRegion(moleculeContainer, _haloInfo[p]._leavingLow, _haloInfo[p]._leavingHigh,
											 _haloInfo[p]._shift, removeFromContainer, LEAVING);
				}
			}

			// then halo particles/copies:
			for (unsigned int p = 0; p < numHaloInfo; p++) {
				collectMoleculesInRegion(moleculeContainer, _haloInfo[p]._copiesLow, _haloInfo[p]._copiesHigh,
						_haloInfo[p]._shift, false, HALO, doHaloPositionCheck);
			}
			break;
		}
		case MessageType::LEAVING_ONLY: {
			Log::global_log->debug() << "sending leaving particles only" << std::endl;
			for(unsigned int p = 0; p < numHaloInfo; p++){
				if (moleculeContainer->isInvalidParticleReturner() and mightUseInvalidParticles) {
					collectLeavingMoleculesFromInvalidParticles(invalidParticles, _haloInfo[p]._leavingLow,
																_haloInfo[p]._leavingHigh, _haloInfo[p]._shift);
				} else {
					collectMoleculesInRegion(moleculeContainer, _haloInfo[p]._leavingLow, _haloInfo[p]._leavingHigh,
											 _haloInfo[p]._shift, removeFromContainer, LEAVING);
				}
			}
			break;
		}
		case MessageType::HALO_COPIES: {
			Log::global_log->debug() << "sending halo particles only" << std::endl;
			for(unsigned int p = 0; p < numHaloInfo; p++){
				collectMoleculesInRegion(moleculeContainer, _haloInfo[p]._copiesLow, _haloInfo[p]._copiesHigh,
						_haloInfo[p]._shift, false, HALO, doHaloPositionCheck);
			}
			break;
		}
		case MessageType::FORCES: {
			Log::global_log->debug() << "sending forces" << std::endl;
			for(unsigned int p = 0; p < numHaloInfo; p++){
				collectMoleculesInRegion(moleculeContainer, _haloInfo[p]._leavingLow, _haloInfo[p]._leavingHigh,
					_haloInfo[p]._shift, false, FORCES);
			}
			break;
		}
		default:
			std::ostringstream error_message;
			error_message << "[CommunicationPartner] MessageType unknown!" << std::endl;
			MARDYN_EXIT(error_message.str());
	}

	#ifndef NDEBUG
		const int numLeaving = _sendBuf.getNumLeaving();
		const int numHalo = _sendBuf.getNumHalo();
		Log::global_log->debug() << "Buffer contains " << numLeaving << " leaving particles with IDs " << std::endl;
		std::ostringstream buf1;
		for (int i = 0; i < numLeaving; ++i) {
			Molecule m;
			_sendBuf.readLeavingMolecule(i, m);
			buf1 << m.getID() << " ";
		}
		Log::global_log->debug() << buf1.str() << std::endl;

		Log::global_log->debug() << "and " << numHalo << " halo particles with IDs " << std::endl;
		std::ostringstream buf2;
		for (int i = 0; i < numHalo; ++i) {
			Molecule m;
			_sendBuf.readHaloMolecule(i, m);
			buf2 << m.getID() << " ";
		}
		Log::global_log->debug() << buf2.str() << std::endl;


	#endif

	MPI_CHECK(MPI_Isend(_sendBuf.getDataForSending(), (int ) _sendBuf.getNumElementsForSending(), _sendBuf.getMPIDataType(), _rank, 99, comm, _sendRequest));
	_msgSent = false;
	_isSending = true;
}

bool CommunicationPartner::testSend() {
	if (not _msgSent) {
		int flag = 0;
		MPI_CHECK(MPI_Test(_sendRequest, &flag, _sendStatus)); // THIS CAUSES A SEG FAULT IN PUSH_PULL_NEIGHBOURS
		if (flag == 1) {
			_msgSent = true;
			_isSending = false;
			_sendBuf.clear();
		}
	}
	return _msgSent;
}

void CommunicationPartner::resetReceive() {
	_countReceived = _msgReceived = _isReceiving = false;
}

bool CommunicationPartner::iprobeCount(const MPI_Comm& comm, const MPI_Datatype& /*type*/) {
	if (not _countReceived) {
		_isReceiving = true;
		int flag = 0;
		MPI_CHECK(MPI_Iprobe(_rank, 99, comm, &flag, _recvStatus));
		if (flag != 0) {
			_countReceived = true;
			_countTested = 0;
			int numrecv;
			MPI_CHECK(MPI_Get_count(_recvStatus, _sendBuf.getMPIDataType(), &numrecv));
                        #ifndef NDEBUG
                                Log::global_log->debug() << "Received byteCount from " << _rank << std::endl;
                                Log::global_log->debug() << "Preparing to receive " << numrecv << " bytes." << std::endl;
                        #endif
			_recvBuf.resizeForRawBytes(numrecv);
			MPI_CHECK(MPI_Irecv(_recvBuf.getDataForSending(), numrecv, _sendBuf.getMPIDataType(), _rank, 99, comm, _recvRequest));
		}
	}
	return _countReceived;
}
bool CommunicationPartner::testRecv(ParticleContainer* moleculeContainer, bool removeRecvDuplicates, bool force) {
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
		if (flag != 0) {
			_msgReceived = true;
			_isReceiving = false;

			if(!force) { // Buffer is particle data

				unsigned long numHalo, numLeaving;
				_recvBuf.resizeForReceivingMolecules(numLeaving, numHalo);

#ifndef NDEBUG
				Log::global_log->debug() << "Receiving particles from " << _rank << std::endl;
				Log::global_log->debug() << "Buffer contains " << numLeaving << " leaving particles with IDs " << std::endl;
				std::ostringstream buf1;
				for (unsigned long i = 0; i < numLeaving; ++i) {
					Molecule m;
					_recvBuf.readLeavingMolecule(i, m);
					buf1 << m.getID() << " ";
				}
				Log::global_log->debug() << buf1.str() << std::endl;

				Log::global_log->debug() << "and " << numHalo << " halo particles with IDs " << std::endl;
				std::ostringstream buf2;
				for (unsigned long i = 0; i < numHalo; ++i) {
					Molecule m;
					_recvBuf.readHaloMolecule(i, m);
					buf2 << m.getID() << " ";
				}
				Log::global_log->debug() << buf2.str() << std::endl;
#endif

				global_simulation->timers()->start("COMMUNICATION_PARTNER_TEST_RECV");
				unsigned long totalNumMols = numLeaving + numHalo;


				/*#if defined(_OPENMP) and not defined (ADVANCED_OVERLAPPING)
				#pragma omp parallel for schedule(static)
				#endif*/
				for (unsigned long i = 0; i < totalNumMols; i++) {
					Molecule m;
					if (i < numLeaving) {
						// leaving
						_recvBuf.readLeavingMolecule(i, m);
						moleculeContainer->addParticle(m, false, removeRecvDuplicates);
					} else {
						// halo
						_recvBuf.readHaloMolecule(i - numLeaving, m);
						moleculeContainer->addHaloParticle(m, false, removeRecvDuplicates);
					}
				}
			} else { // Buffer is force data

				/*
				#if defined(_OPENMP)
				#pragma omp for schedule(static)
				#endif
				*/

				unsigned long numForces;
				_recvBuf.resizeForReceivingMolecules(numForces);


#ifndef NDEBUG
				Log::global_log->debug() << "Receiving particles from " << _rank << std::endl;
				Log::global_log->debug() << "Buffer contains " << numForces << " force particles with IDs " << std::endl;
				std::ostringstream buf1;

				for(unsigned long i = 0; i < numForces; ++i) {
					Molecule m;
					_recvBuf.readForceMolecule(i, m);
					buf1 << m.getID() << " ";
				}
				Log::global_log->debug() << buf1.str() << std::endl;


#endif
				global_simulation->timers()->start("COMMUNICATION_PARTNER_TEST_RECV");
				//mols.resize(numForces);

				/*#if defined(_OPENMP) and not defined (ADVANCED_OVERLAPPING)
				#pragma omp parallel for schedule(static)
				#endif*/

				double pos[3];
				decltype(moleculeContainer->getMoleculeAtPosition(pos)) originalPreviousIter{};

				for(unsigned i = 0; i < numForces; ++i) {
					Molecule m;
					_recvBuf.readForceMolecule(i, m);
					//mols[i] = m;
					const double position[3] = { m.r(0), m.r(1), m.r(2) };

					originalPreviousIter =
						addValuesAndGetIterator(moleculeContainer, position, originalPreviousIter, m);
				}

				//moleculeContainer->addParticles(mols, removeRecvDuplicates);

			}



			_recvBuf.clear();
			global_simulation->timers()->stop("COMMUNICATION_PARTNER_TEST_RECV");


		} else {
			++_countTested;
		}
	}
	return _msgReceived;
}

void CommunicationPartner::initRecv(int numParticles, const MPI_Comm& comm, const MPI_Datatype& type) {
	// one single call from KDDecomposition::migrate particles.
	// So all molecules, which arrive area leaving molecules.
	_countReceived = true;
	_isReceiving = true;

	// hackaround - resizeForAppendingLeavingMolecules is intended for the send-buffer, not the recv one.
	_recvBuf.resizeForAppendingLeavingMolecules(numParticles);

	MPI_CHECK(MPI_Irecv(_recvBuf.getDataForSending(), _recvBuf.getNumElementsForSending(), _sendBuf.getMPIDataType(), _rank, 99, comm, _recvRequest));
}

void CommunicationPartner::deadlockDiagnosticSendRecv() {

	deadlockDiagnosticSend();

	if (not _countReceived and _isReceiving) {
		Log::global_log->warning() << "Probe request to " << _rank << " not yet completed" << std::endl;
	}

	deadlockDiagnosticRecv();
}

void CommunicationPartner::deadlockDiagnosticSend() {
	// intentionally using std::cout instead of Log::global_log, we want the messages from all processes
	if (not _msgSent and _isSending) {
		Log::global_log->warning() << "Send request to " << _rank << " not yet completed" << std::endl;
	}
}

void CommunicationPartner::deadlockDiagnosticRecv() {
	if (not _msgReceived and _isReceiving) {
		Log::global_log->warning() << "Recv request to " << _rank << " not yet completed" << std::endl;
	}
}

void CommunicationPartner::add(CommunicationPartner partner) {
	mardyn_assert(partner._rank == _rank);
	_haloInfo.push_back(partner._haloInfo[0]);
}

void CommunicationPartner::collectMoleculesInRegion(ParticleContainer* moleculeContainer, const double lowCorner[3],
													const double highCorner[3], const double shift[3],
													const bool removeFromContainer,
													const HaloOrLeavingCorrection haloLeaveCorr,
													bool doHaloPositionCheck) {
	using std::vector;
	global_simulation->timers()->start("COMMUNICATION_PARTNER_INIT_SEND");
	std::vector<std::vector<Molecule>> threadData;
	std::vector<int> prefixArray;

	// compute how many molecules are already in of this type: - adjust for Forces
	unsigned long numMolsAlreadyIn = 0;
	if (haloLeaveCorr == HaloOrLeavingCorrection::LEAVING) {
		numMolsAlreadyIn = _sendBuf.getNumLeaving();
	} else if (haloLeaveCorr == HaloOrLeavingCorrection::HALO) {
		numMolsAlreadyIn = _sendBuf.getNumHalo();
	} else if(haloLeaveCorr == HaloOrLeavingCorrection::FORCES) {
		numMolsAlreadyIn = _sendBuf.getNumForces();
	}

	#if defined (_OPENMP)
	#pragma omp parallel shared(threadData, numMolsAlreadyIn)
	#endif
	{
		// in the case of an autopas container, we only want to iterate over inner particle cells if we are sending
		// halos (this assumes the direct communication scheme!)
		ParticleIterator::Type iteratorType =
			haloLeaveCorr == HaloOrLeavingCorrection::HALO and moleculeContainer->isInvalidParticleReturner()
				? ParticleIterator::Type::ONLY_INNER_AND_BOUNDARY
				: ParticleIterator::Type::ALL_CELLS;
		const int numThreads = mardyn_get_num_threads();
		const int threadNum = mardyn_get_thread_num();
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

		for (auto i = moleculeContainer->regionIterator(lowCorner, highCorner, iteratorType); i.isValid(); ++i) {
			//traverse and gather all molecules in the cells containing part of the box specified as parameter
			//i is a pointer to a Molecule; (*i) is the Molecule
			threadData[threadNum].push_back(*i);
			mardyn_assert(i->inBox(lowCorner, highCorner));
			if (removeFromContainer) {
				moleculeContainer->deleteMolecule(i, false);
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
			if (haloLeaveCorr == HaloOrLeavingCorrection::LEAVING) {
				_sendBuf.resizeForAppendingLeavingMolecules(totalNumMolsAppended);
			} else if (haloLeaveCorr == HaloOrLeavingCorrection::HALO) {
				_sendBuf.resizeForAppendingHaloMolecules(totalNumMolsAppended);
			} else if (haloLeaveCorr == HaloOrLeavingCorrection::FORCES) {
				_sendBuf.resizeForAppendingForceMolecules(totalNumMolsAppended);
			}
		}

		#if defined (_OPENMP)
		#pragma omp barrier
		#endif


		Domain* domain = global_simulation->getDomain();

		//reduce the molecules in the send buffer and also apply the shift
		int myThreadMolecules = prefixArray[threadNum + 1] - prefixArray[threadNum];
		for(int i = 0; i < myThreadMolecules; i++){
			Molecule mCopy = threadData[threadNum][i];
			mCopy.move(0, shift[0]);
			mCopy.move(1, shift[1]);
			mCopy.move(2, shift[2]);
			for (int dim = 0; dim < 3; dim++) {
				if (haloLeaveCorr == HaloOrLeavingCorrection::HALO and doHaloPositionCheck) {
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
				} else if (haloLeaveCorr == HaloOrLeavingCorrection::LEAVING || haloLeaveCorr == HaloOrLeavingCorrection::FORCES) {
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
				}
				/* else if(haloLeaveCorr == HaloOrLeavingCorrection::FORCES) { // using Leaving correction for now.
					// THIS IS STILL MISSING!
				} */
			} /* for-loop dim */
			if (haloLeaveCorr == HaloOrLeavingCorrection::LEAVING) {
				_sendBuf.addLeavingMolecule(numMolsAlreadyIn + prefixArray[threadNum] + i, mCopy);
			} else if (haloLeaveCorr == HaloOrLeavingCorrection::HALO) {
				_sendBuf.addHaloMolecule(numMolsAlreadyIn + prefixArray[threadNum] + i, mCopy);
			} else if (haloLeaveCorr == HaloOrLeavingCorrection::FORCES) {
				_sendBuf.addForceMolecule(numMolsAlreadyIn + prefixArray[threadNum] + i, mCopy);
			}
		}
	}
	global_simulation->timers()->stop("COMMUNICATION_PARTNER_INIT_SEND");
}

void CommunicationPartner::collectLeavingMoleculesFromInvalidParticles(std::vector<Molecule>& invalidParticles, double* lowCorner,
                                                                       double* highCorner, double* shift) {

	global_simulation->timers()->start("COMMUNICATION_PARTNER_INIT_SEND");

	// compute how many molecules are already in of this type: - adjust for Forces

	const auto removeBegin = std::partition(invalidParticles.begin(), invalidParticles.end(), [=](const Molecule& m) {
	  // if this returns true, the particle will be put in the first part of the partition, else in the second.
	  return not m.inBox(lowCorner, highCorner);
	});

	auto numMolsAlreadyIn = _sendBuf.getNumLeaving();
	const auto totalNumMolsAppended = invalidParticles.end() - removeBegin;
	// resize the send-buffer
	_sendBuf.resizeForAppendingLeavingMolecules(totalNumMolsAppended);

	Domain* domain = global_simulation->getDomain();

	// This lambda has the following functionality:
	// it will add the given molecule to _sendBuf with the necessary shift.
	auto shiftAndAdd = [domain, lowCorner, highCorner, shift, this, &numMolsAlreadyIn](Molecule& m) {
		if (not m.inBox(lowCorner, highCorner)) {
			std::ostringstream error_message;
			error_message << "trying to remove a particle that is not in the halo region" << std::endl;
			MARDYN_EXIT(error_message.str());
		}
		for (int dim = 0; dim < 3; dim++) {
			if (shift[dim] != 0) {
				m.setr(dim, m.r(dim) + shift[dim]);
				// some additional shifting to ensure that rounding errors do not hinder the correct
				// placement
				if (shift[dim] < 0) {    // if the shift was negative, it is now in the lower part of the domain -> min
					if (m.r(dim) < 0) {  // in the lower part it was wrongly shifted if
						m.setr(dim, 0);  // ensures that r is at least the boundingBoxMin
					}
				} else {                                             // shift > 0
					if (m.r(dim) >= domain->getGlobalLength(dim)) {  // in the lower part it was wrongly shifted if
						// std::nexttoward: returns the next bigger value of _boundingBoxMax
						vcp_real_calc r = domain->getGlobalLength(dim);
						m.setr(dim, std::nexttoward(r, r - 1.f));  // ensures that r is smaller than the boundingBoxMax
					}
				}
			}
		}
		_sendBuf.addLeavingMolecule(numMolsAlreadyIn, m);
		++numMolsAlreadyIn;
	};

	// now insert all particles that are in the second partition to _sendBuf.
	std::for_each(removeBegin, invalidParticles.end(), shiftAndAdd);

	// remove the now already processed particles from the invalidParticles vector.
	invalidParticles.erase(removeBegin, invalidParticles.end());

	global_simulation->timers()->stop("COMMUNICATION_PARTNER_INIT_SEND");

}

size_t CommunicationPartner::getDynamicSize() {
	return _sendBuf.getDynamicSize() + _recvBuf.getDynamicSize() + _haloInfo.capacity() * sizeof(PositionInfo);
}

void CommunicationPartner::print(std::ostream& stream) const {
	stream << "Partner rank: "<< _rank << std::endl;
	stream << "Halo regions: " << std::endl;
	for(auto& region : _haloInfo){
		stream << " Region:" << std::endl;
		stream << "  both: "
				<< " [" << region._bothLow[0] << ", " << region._bothHigh[0] << ") x "
				<< " [" << region._bothLow[1] << ", " << region._bothHigh[1] << ") x "
				<< " ["	<< region._bothLow[2] << ", " << region._bothHigh[2] << ")" << std::endl;
		stream << "  leaving: "
				<< " [" << region._leavingLow[0] << ", " << region._leavingHigh[0] << ") x "
				<< " [" << region._leavingLow[1] << ", " << region._leavingHigh[1] << ") x "
				<< " ["	<< region._leavingLow[2] << ", " << region._leavingHigh[2] << ")" << std::endl;
		stream << "  copies: "
				<< " [" << region._copiesLow[0] << ", " << region._copiesHigh[0] << ") x "
				<< " [" << region._copiesLow[1] << ", " << region._copiesHigh[1] << ") x "
				<< " ["	<< region._copiesLow[2] << ", " << region._copiesHigh[2] << ")" << std::endl;
		stream << "  offset: (" << region._offset[0] << ", " << region._offset[1] << ", " << region._offset[2] << ")" << std::endl;
		stream << "  shift:	(" << region._shift[0] << ", " << region._shift[1]	<< ", "<< region._shift[2] << ")" << std::endl;
	}
}

