/*
 * CommunicationBuffer.cpp
 *
 *  Created on: 12 July 2018
 *      Author: fernanor
 */

#ifdef ENABLE_MPI
#include "ResilienceComm.h"

#include <climits> /* UINT64_MAX */

#include <sstream>
#include "utils/Logger.h"

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
						  int const numberOfBackups,
						  size_t const sizePerRank,
	                      std::vector<int>& backing, 
	                      std::vector<int>& backedBy, 
	                      std::vector<int>& backingTags, 
	                      std::vector<int>& backedByTags) {
	size_t totalBytesRecv = sizePerRank*numberOfBackups*sizeof(int);
	std::vector<char> recvArray(totalBytesRecv);
	if (_rank == 0) {
		std::stringstream bkinf;
		// bkinf << "    RR: backupInfo: " << totalBytesRecv << "\n";
		// for (auto const& rnk : backupInfo) {
		// 	bkinf << rnk << ", ";
		// }
		// global_log->info() << bkinf.str() << std::endl;
		mardyn_assert(static_cast<uint>(sizePerRank*numberOfBackups*_numProcs) == backupInfo.size());
	}
	else {
		mardyn_assert(backupInfo.empty());
	}
	constexpr int const scatteringRank = 0;
	int mpi_error =	MPI_Scatter(reinterpret_cast<char*>(backupInfo.data()),
	                            totalBytesRecv,
	 			                MPI_CHAR,
				                recvArray.data(),
				                totalBytesRecv,
				                MPI_CHAR,
				                scatteringRank,
				                MPI_COMM_WORLD);
	mardyn_assert(mpi_error == MPI_SUCCESS);
	backing.resize(numberOfBackups);
	backedBy.resize(numberOfBackups);
	backingTags.resize(numberOfBackups);
	backedByTags.resize(numberOfBackups);
	auto backingAsChar = reinterpret_cast<char*>(backing.data());
	auto backedByAsChar = reinterpret_cast<char*>(backedBy.data());
	auto backingTagsAsChar = reinterpret_cast<char*>(backingTags.data());
	auto backedByTagsAsChar = reinterpret_cast<char*>(backedByTags.data());
	std::copy(recvArray.begin()+0*totalBytesRecv/sizePerRank, recvArray.begin()+1*totalBytesRecv/sizePerRank, backingAsChar);
	std::copy(recvArray.begin()+1*totalBytesRecv/sizePerRank, recvArray.begin()+2*totalBytesRecv/sizePerRank, backedByAsChar);
	std::copy(recvArray.begin()+2*totalBytesRecv/sizePerRank, recvArray.begin()+3*totalBytesRecv/sizePerRank, backingTagsAsChar);
	std::copy(recvArray.begin()+3*totalBytesRecv/sizePerRank, recvArray.begin()+4*totalBytesRecv/sizePerRank, backedByTagsAsChar);
	// global_log->info() << "    RR: Dumping scattered backup info: " << std::endl;
	// global_log->set_mpi_output_all();
	// std::stringstream bckd, bckBy, bckdTags, bckByTags;
	// for (int i=0; i<numberOfBackups; ++i) {
	// 	mardyn_assert(backing[i]<_numProcs);
	// 	mardyn_assert(backedBy[i]<_numProcs);
	// 	bckd << backing[i] << ", ";
	// 	bckBy << backedBy[i] << ", ";
	// 	bckdTags << backingTags[i] << ", ";
	// 	bckByTags << backedByTags[i] << ", ";
	// }
	// global_log->info() << "        Backed: " << bckd.str() << " Backed by: " << bckBy.str() << std::endl;
	// global_log->info() << "        Backed tags: " << bckdTags.str() << " Backed by tags: " << bckByTags.str() << std::endl;
	return 0;
}

int ResilienceComm::exchangeSnapshotSizes(
		std::vector<int>& backing,
		std::vector<int>& backedBy,
		std::vector<int>& backingTags,
		std::vector<int>& backedByTags,
		size_t const snapshotSize,
		std::vector<int>& backupDataSizes) {
	// send the size of this snapshot to all ranks backing it
	int src = -1;
	int dest = -1;
	int tag = -1;
	int status = MPI_ERR_UNKNOWN;
	for (size_t ib=0; ib<backedBy.size(); ++ib) {
		MPI_Request request=0;
		dest = backedBy[ib];
		tag = backedByTags[ib];
		// status = MPI_Isend(&snapshotSize, sizeof(snapshotSize), MPI_CHAR, dest, tag, MPI_COMM_WORLD, &request);
		status = MPI_Bsend(&snapshotSize, sizeof(snapshotSize), MPI_CHAR, dest, tag, MPI_COMM_WORLD);
		mardyn_assert(status == MPI_SUCCESS);
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	// setup the receiving buffers too for all ranks the current one is backing
	for (size_t ib=0; ib<backing.size(); ++ib) {
		MPI_Status recvStatus; 
		src = backing[ib];
		tag = backingTags[ib];
		void* target = &(backupDataSizes.data()[ib]);
		// status = MPI_Irecv(target, sizeof(snapshotSize), MPI_CHAR, src, tag, MPI_COMM_WORLD, &request);
		status = MPI_Recv(&(backupDataSizes.data()[ib]), sizeof(snapshotSize), MPI_CHAR, src, tag, MPI_COMM_WORLD, &recvStatus);
		mardyn_assert(status == MPI_SUCCESS);
	}
	mardyn_assert(status == MPI_SUCCESS);
	return 0;
}

int ResilienceComm::exchangeSnapshots(
		std::vector<int>& backing,
		std::vector<int>& backedBy,
		std::vector<int>& backingTags,
		std::vector<int>& backedByTags,
		std::vector<int>& backupDataSizes,
		std::vector<char>& sendData,
		std::vector<char>& recvData) {
	// send the snapshot to all ranks backing it
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
		dest = backedBy[ib];
		tag = backedByTags[ib];
		// global_log->info() << "    RR: Sending " << sendData.size()
		// 		<< " bytes to: " << dest << " using tag: " << tag << std::endl;
		status = MPI_Bsend(sendData.data(), sendData.size(), MPI_CHAR, dest, tag, MPI_COMM_WORLD);
		mardyn_assert(status == MPI_SUCCESS);
	}
	// setup the receiving buffers too for all ranks the current one is backing
	for (size_t ib=0; ib<backing.size(); ++ib) {
		MPI_Status recvStatus;
		src = backing[ib];
		tag = backingTags[ib];
		size_t const recvIndex = recvIndices[ib];
		// global_log->info() << "    RR: Receiving " 
		// 		<< backupDataSizes[ib] << " bytes from " 
		// 		<< src << " at " 
		// 		<< recvIndices[ib] << " using tag: " << tag << std::endl;
		status = MPI_Recv(&recvData.data()[recvIndex], backupDataSizes[ib], MPI_CHAR, src, tag, MPI_COMM_WORLD, &recvStatus);
		mardyn_assert(status == MPI_SUCCESS);
		//verify a bunch of stuff
		int count;
		MPI_Get_count(&recvStatus, MPI_CHAR, &count);
		mardyn_assert(count == backupDataSizes[ib]);
		mardyn_assert(recvStatus.MPI_SOURCE == src);
		mardyn_assert(recvStatus.MPI_TAG == tag);
		mardyn_assert(status == MPI_SUCCESS);
	}
	return 0;
}
#endif /* ENABLE_MPI */