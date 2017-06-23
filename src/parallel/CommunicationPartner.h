/*
 * CommunicationPartner.h
 *
 *  Created on: Nov 16, 2015
 *      Author: tchipevn
 */

#ifndef COMMUNICATIONPARTNER_H_
#define COMMUNICATIONPARTNER_H_

#include "mpi.h"
#include <vector>
#include <type_traits>
#include "ParticleDataForwardDeclaration.h"
#include "utils/Logger.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "ParticleData.h"
#include "ParticleForceData.h"
#include "WrapOpenMP.h"
#include <string>

using Log::global_log;

typedef enum {
	LEAVING_AND_HALO_COPIES = 0, /** send process-leaving particles and halo-copies together in one message */
	HALO_COPIES = 1, /** send halo-copies only */
	LEAVING_ONLY = 2, /** send process-leaving particles only */
	FORCES = 3 //< send forces
} MessageType;

class ParticleContainer;
class ParticleForceData;

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

	template<typename BufferType = ParticleData>
	void initSend(ParticleContainer* moleculeContainer, const MPI_Comm& comm, const MPI_Datatype& type,
			MessageType msgType, bool removeFromContainer = false) {

		auto& sendBuf = getSendBuf<BufferType>();

		global_log->debug() << _rank << std::endl;

		const unsigned int numHaloInfo = _haloInfo.size();
		switch (msgType){
			case LEAVING_AND_HALO_COPIES: {
				//"This message type requires a ParticleData buffer."
				mardyn_assert(!isForceData<BufferType>());
				global_log->debug() << "sending halo and boundary particles together" << std::endl;
				for(unsigned int p = 0; p < numHaloInfo; p++){
					collectMoleculesInRegion<BufferType>(moleculeContainer, _haloInfo[p]._bothLow, _haloInfo[p]._bothHigh, _haloInfo[p]._shift);
				}
				break;
			}
			case LEAVING_ONLY: {
				//"This message type requires a ParticleData buffer."
				mardyn_assert(!isForceData<BufferType>());
				global_log->debug() << "sending leaving particles only" << std::endl;
				for(unsigned int p = 0; p < numHaloInfo; p++){
					collectMoleculesInRegion<BufferType>(moleculeContainer, _haloInfo[p]._leavingLow, _haloInfo[p]._leavingHigh, _haloInfo[p]._shift, removeFromContainer);
				}
				break;
			}
			case HALO_COPIES: {
				//"This message type requires a ParticleData buffer."
				mardyn_assert(!isForceData<BufferType>());
				global_log->debug() << "sending halo particles only" << std::endl;
				for(unsigned int p = 0; p < numHaloInfo; p++){
					collectMoleculesInRegion<BufferType>(moleculeContainer, _haloInfo[p]._copiesLow, _haloInfo[p]._copiesHigh, _haloInfo[p]._shift);
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
			const int n = sendBuf.size();
			global_log->debug() << "Buffer contains " << n << " particles with IDs " << std::endl;
			std::ostringstream buf;
			for (int i = 0; i < n; ++i) {
				buf << sendBuf[i].id << " ";
			}
			global_log->debug() << buf.str() << std::endl;
		#endif

		MPI_CHECK(MPI_Isend(sendBuf.data(), (int ) sendBuf.size(), type, _rank, 99, comm, _sendRequest));
		_msgSent = _countReceived = _msgReceived = false;
	}

	template<typename BufferType = ParticleData>
	bool testSend() {

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

	template<typename BufferType = ParticleData>
	bool iprobeCount(const MPI_Comm& comm, const MPI_Datatype& type) {

		auto& recvBuf = getRecvBuf<BufferType>();

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
				recvBuf.resize(numrecv);
				MPI_CHECK(MPI_Irecv(recvBuf.data(), numrecv, type, _rank, 99, comm, _recvRequest));
			}
		}
		return _countReceived;
	}

	template<typename BufferType = ParticleData >
	bool testRecv(ParticleContainer* moleculeContainer, bool removeRecvDuplicates) {
		using Log::global_log;

		auto& recvBuf = getRecvBuf<BufferType>();

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
				int numrecv = recvBuf.size();

				#ifndef NDEBUG
					global_log->debug() << "Receiving particles from " << _rank << std::endl;
					global_log->debug() << "Buffer contains " << numrecv << " particles with IDs " << std::endl;
				#endif

				//Code for different buffer types
				testRecvHandle<BufferType>(moleculeContainer, removeRecvDuplicates, numrecv);

			} else {
				++_countTested;
			}
		}
		return _msgReceived;
	}

	//! Handle receive for ParticleData
	template <class BufferType>
	void testRecvHandle(ParticleContainer* moleculeContainer, bool removeRecvDuplicates, int numrecv){

	}

	template<typename BufferType = ParticleData >
	void initRecv(int numParticles, const MPI_Comm& comm, const MPI_Datatype& type){

		auto& recvBuf = getRecvBuf<BufferType>();

		_countReceived = true;
		recvBuf.resize(numParticles);
		MPI_CHECK(MPI_Irecv(recvBuf.data(), numParticles, type, _rank, 99, comm, _recvRequest));
	}

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

private:
	template<typename BufferType = ParticleData>
	void collectMoleculesInRegion(ParticleContainer* moleculeContainer, const double lowCorner[3], const double highCorner[3], const double shift[3],
			const bool removeFromContainer = false){
		using std::vector;

		auto& sendBuf = getSendBuf<BufferType>();

		timer("COMMUNICATION_PARTNER_INIT_SEND", true);
		int prevNumMols = sendBuf.size();
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
				sendBuf.resize(prevNumMols + totalNumMols);
			}

			#if defined (_OPENMP)
			#pragma omp barrier
			#endif

			//reduce the molecules in the send buffer and also apply the shift
			int myThreadMolecules = prefixArray[threadNum + 1] - prefixArray[threadNum];
			for(int i = 0; i < myThreadMolecules; i++){
				BufferType m;
				BufferType::MoleculeToParticleData(m, threadData[threadNum][i]);
				m.r[0] += shift[0];
				m.r[1] += shift[1];
				m.r[2] += shift[2];

				sendBuf[prevNumMols + prefixArray[threadNum] + i] = m;
			}
		}
		timer("COMMUNICATION_PARTNER_INIT_SEND", false);
	}

	//! Decide if T is ParticleForceData or ParticleData. Fail if T is something else.
	template<typename BufferType>
	static constexpr inline typename std::enable_if<std::is_same<BufferType, ParticleForceData>::value, bool>::type isForceData(){
		return true;
	}

	template<typename BufferType>
	static constexpr inline typename std::enable_if<std::is_same<BufferType, ParticleData>::value, bool>::type isForceData() {
		return false;
	}

	//! Returns the receive buffer for a specific buffer type
	template<typename BufferType>
	inline typename std::enable_if<std::is_same<BufferType, ParticleForceData>::value, std::vector<BufferType>&>::type getRecvBuf(){
		return _recvBufForces;
	}

	template<typename BufferType>
	inline typename std::enable_if<std::is_same<BufferType, ParticleData>::value, std::vector<BufferType>&>::type getRecvBuf(){
		return _recvBuf;
	}

	//! Returns the send buffer for a specific buffer type
	template<typename BufferType>
	inline typename std::enable_if<std::is_same<BufferType, ParticleForceData>::value, std::vector<BufferType>&>::type  getSendBuf(){
		return _sendBufForces;
	}

	template<typename BufferType>
	inline typename std::enable_if<std::is_same<BufferType, ParticleData>::value, std::vector<BufferType>&>::type  getSendBuf(){
		return _sendBuf;
	}

	//! As global_simulation is only available in the cpp file, start and stop timer there.
	//! @param name Name of the timer
	//! @param start Start (true) or stop (false) the timer
	void timer(const std::string name, const bool start);

	int _rank;
    int _countTested;
	std::vector<PositionInfo> _haloInfo;

	// technical variables
	MPI_Request *_sendRequest, *_recvRequest;
	MPI_Status *_sendStatus, *_recvStatus;
	std::vector<ParticleData> _sendBuf, _recvBuf;
	std::vector<ParticleForceData> _sendBufForces, _recvBufForces;
	bool _msgSent, _countReceived, _msgReceived;

};

//typename std::enable_if<std::is_same<BufferType, ParticleData>::value, void>::type
template <>
void CommunicationPartner::testRecvHandle<ParticleData>(ParticleContainer* moleculeContainer, bool removeRecvDuplicates, int numrecv) ;

//! Handle receive for ParticleForceData

//typename std::enable_if<std::is_same<BufferType, ParticleForceData>::value, void>::type
template<>
void CommunicationPartner::testRecvHandle<ParticleForceData>(ParticleContainer* moleculeContainer, bool removeRecvDuplicates, int numrecv) ;

#endif /* COMMUNICATIONPARTNER_H_ */
