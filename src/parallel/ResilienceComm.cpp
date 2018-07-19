/*
 * CommunicationBuffer.cpp
 *
 *  Created on: 12 July 2018
 *      Author: fernanor
 */

#include "ResilienceComm.h"

#include <climits> /* UINT64_MAX */

#include <sstream>
#include "utils/Logger.h"
using Log::global_log;

//pretty much default
ResilienceComm::ResilienceComm(int numProcs, int rank)
		: _numProcs(numProcs)
		, _rank(rank) {
}

//copy constructor
ResilienceComm::ResilienceComm(ResilienceComm const& rc)
		: _numProcs(rc._numProcs)
		, _rank(rc._rank) {
}

ResilienceComm::~ResilienceComm() {
	/*do nothing*/
}

int ResilienceComm::scatterBackupInfo(std::vector<int>& backupInfo, 
	                      std::vector<int>& backing, 
	                      std::vector<int>& backedBy, 
						  int const numberOfBackups) {
	size_t totalBytesSend = backupInfo.size()*sizeof(int);
	size_t totalBytesRecv = 2*numberOfBackups*sizeof(int);
	std::unique_ptr<char[]> sendArray(new char[totalBytesSend]);
	std::unique_ptr<char[]> recvArray(new char[totalBytesRecv]);
	if (_rank == 0) {
		mardyn_assert(static_cast<uint>(2*numberOfBackups*_numProcs) == backupInfo.size());
		memcpy(reinterpret_cast<void*>(sendArray.get()), reinterpret_cast<void*>(backupInfo.data()), totalBytesSend);
	}
	else {
		memset(sendArray.get(), 0, totalBytesSend);
	}
	memset(recvArray.get(), 0, totalBytesRecv);
	int mpi_error =	MPI_Scatter(sendArray.get(),
	                            totalBytesRecv,        //the call expects the number of bytes to send PER RANK
	 			                MPI_CHAR,
				                recvArray.get(),
				                totalBytesRecv,
				                MPI_CHAR,
				                0,
				                MPI_COMM_WORLD);
	mardyn_assert(mpi_error == MPI_SUCCESS);
	backing.resize(numberOfBackups);
	backedBy.resize(numberOfBackups);
	memcpy(backing.data(), reinterpret_cast<void*>(recvArray.get()), totalBytesRecv/2);
	memcpy(backedBy.data(), reinterpret_cast<void*>(recvArray.get()+totalBytesRecv/2), totalBytesRecv/2);
	for (int i=0; i<numberOfBackups; ++i) {
		mardyn_assert(backing[i]<_numProcs);
		mardyn_assert(backedBy[i]<_numProcs);
	}
	// (re-)init the request buffers
	_exchangeSizesRequests.clear();
	_exchangeSnapshotsRequests.clear();
	_exchangeSizesRequests.resize(backing.size()+backedBy.size());
	_exchangeSnapshotsRequests.resize(backing.size()+backedBy.size());
	return 0;
}

int ResilienceComm::exchangeSnapshotSizes(
		std::vector<int>& backing,
		std::vector<int>& backedBy,
		size_t const snapshotSize,
		std::vector<int>& backupDataSizes) {
	// send the size of this snapshot to all ranks backing it
	int src = -1;
	int dest = -1;
	int tag = -1;
	int status = MPI_ERR_UNKNOWN;
	for (size_t ib=0; ib<backedBy.size(); ++ib) {
		MPI_Request request;
		dest = backedBy[ib];
		tag = _rank+dest; // maybe a better tag 
		status = MPI_Isend(&snapshotSize, sizeof(snapshotSize), MPI_CHAR, dest, tag, MPI_COMM_WORLD, &request);
		mardyn_assert(status == MPI_SUCCESS);
		mardyn_assert(ib < _exchangeSizesRequests.size());
		_exchangeSizesRequests[ib] = request;
	}
	// setup the receiving buffers too for all ranks the current one is backing
	for (size_t ib=0; ib<backing.size(); ++ib) {
		MPI_Request request;
		src = backing[ib];
		tag = _rank+src;
		status = MPI_Irecv(&backupDataSizes.data()[ib], sizeof(snapshotSize), MPI_CHAR, src, tag, MPI_COMM_WORLD, &request);
		mardyn_assert(status == MPI_SUCCESS);
		mardyn_assert(ib+backedBy.size() < _exchangeSizesRequests.size());
		_exchangeSizesRequests[ib+backedBy.size()] = request;
	}
	status = MPI_Waitall(_exchangeSizesRequests.size(), _exchangeSizesRequests.data(), MPI_STATUSES_IGNORE);
	mardyn_assert(status == MPI_SUCCESS);
	return 0;
}

int ResilienceComm::exchangeSnapshots(
		std::vector<int>& backing,
		std::vector<int>& backedBy,
		std::vector<int>& backupDataSizes,
		std::vector<char>& sendData,
		std::vector<char>& recvData) {
	// send the size of this snapshot to all ranks backing it
	int src = -1;
	int dest = -1;
	int tag = -1;
	int status = MPI_ERR_UNKNOWN;
	//prepare memory, create prefix sum
	std::vector<int> recvIndices(backupDataSizes.size());
	recvIndices[0] = 0; //first index always 0
	auto dstIt = recvIndices.begin()+1;
	auto srcIt = backupDataSizes.begin();
	while (dstIt != recvIndices.end()) {
		*dstIt = *(dstIt-1)+*srcIt;
		++dstIt; ++srcIt;
	}
	size_t const totalRecvSize = recvIndices.back()+*srcIt;
	recvData.resize(totalRecvSize);
	for (size_t ib=0; ib<backedBy.size(); ++ib) {
		MPI_Request request;
		dest = backedBy[ib];
		tag = _rank+dest; // maybe a better tag 
		status = MPI_Isend(sendData.data(), sendData.size(), MPI_CHAR, dest, tag, MPI_COMM_WORLD, &request);
		mardyn_assert(status == MPI_SUCCESS);
		mardyn_assert(ib < _exchangeSnapshotsRequests.size());
		_exchangeSnapshotsRequests[ib] = request;
	}
	// setup the receiving buffers too for all ranks the current one is backing
	for (size_t ib=0; ib<backing.size(); ++ib) {
		MPI_Request request;
		src = backing[ib];
		tag = _rank+src;
		size_t const recvIndex = recvIndices[ib];
		status = MPI_Irecv(&recvData.data()[recvIndex], backupDataSizes[ib], MPI_CHAR, src, tag, MPI_COMM_WORLD, &request);
		mardyn_assert(status == MPI_SUCCESS);
		mardyn_assert(ib+backedBy.size() < _exchangeSnapshotsRequests.size());
		_exchangeSnapshotsRequests[ib+backedBy.size()] = request;
	}
	status = MPI_Waitall(_exchangeSnapshotsRequests.size(), _exchangeSnapshotsRequests.data(), MPI_STATUSES_IGNORE);
	mardyn_assert(status == MPI_SUCCESS);
	return 0;
}
