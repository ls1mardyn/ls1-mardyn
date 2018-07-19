/*
 * CommunicationBuffer.h
 *
 *  Created on: 11 July 2018
 *      Author: Oliver Fernandes
 */

#ifndef SRC_PARALLEL_RESILIENCECOMM_H_
#define SRC_PARALLEL_RESILIENCECOMM_H_

#include "utils/mardyn_assert.h" 

#include <memory>
#include <vector>
#include <stddef.h>
#include <mpi.h>

// do not uncomment the if, it will break halo copies of the kddecomposition!
//#if (not defined(NDEBUG))
#define LS1_SEND_UNIQUE_ID_FOR_HALO_COPIES
//#pragma message "Compilation info: Unique IDs of Halo-Molecules are always present."
//#endif

/**
 * This class implements all communication needed for setting up the redundancy resilience scheme.
 * 
 * Converts all to CHAR internally.
 *
 * TODO: test how this will work when going to Big Endian/Little Endian architectures,
 * due to CHAR conversion.
 *
 */
class ResilienceComm {

public:
	explicit ResilienceComm(int numProcs, int rank);     // constructor
	ResilienceComm() = delete;                           // delete default constructor
	ResilienceComm(ResilienceComm const& rc);            // copy constructor
	~ResilienceComm();                                   // destructor
	/**
	 * Scatters which rank is backing which and the inverse, i.e. which rank is being backed by which.
	 * Scatters info about which rank should backup which. This data is only produced on the root node, and needs
	 * to be scattered. The method assumes that: #ranks being backed by current rank = #ranks backing current rank = numberOfBackups
	 * The method modifies the request buffers (RedundancyComm::_exchangeSizesRequests and RedundancyComm::_exchangeSnapshotRequests)
	 * @param[in] backupInfo Contains the rank configuration. Set up in RedundancyResilience::_determineBackups, check there for layout.
	 * @param[out] backing A list of ranks being backed by the current rank
	 * @param[out] backedBy A list of ranks being backed by the current rank
	 * @param[out] numberOfBackups how many backups are made per rank
	 */
	int scatterBackupInfo(
			std::vector<int>& backupInfo, 
	        std::vector<int>& backing, 
	        std::vector<int>& backedBy, 
			int const numberOfBackups
	);
	/**
	 * Send snapshot size from caller to ranks backing it.
	 * Called by each process after setting up the new snapshot.
	 * Needs to be called before calling RedundancyResilience::exchangeSnapshots. If this was done already
	 * should probably be checked in that method. But it isn't.
	 * Updates request buffer _exchangeSizesRequests.
	 * @param[in] backing A list of ranks being backed by the current rank
	 * @param[in] backedBy A list of ranks the current rank is backing up
	 * @param[in] snapshotSize Size of the local snapshot data in bytes
	 * @param[out] backupDataSizes The individual sizes of the snapshots acquired 
	 */
	int exchangeSnapshotSizes(
			std::vector<int>& backing,
			std::vector<int>& backedBy,
			size_t const snapshotSize, 
			std::vector<int>& backupDataSizes
	);
	/**
	 * Send snapshot meta data from caller to ranks backing it.
	 * Called by each process after trading the backup snapshot sizes.
	 * This will perform the actual sending of snapshot data for each of the backup ranks.
	 * Updates request buffer _exchangeSnapshotRequests.
	 * @param[in] backing A list of ranks being backed by the current rank
	 * @param[in] backedBy A list of ranks the current rank is backing up
	 * @param[in] backupDataSizes The individual sizes of the snapshots acquired 
	 * @param[in] sendData The local snapshot data as char vector
	 * @param[out] recvData The snapshot data of the ranks being backed
	 */
	int exchangeSnapshots(
			std::vector<int>& backing,
			std::vector<int>& backedBy,
			std::vector<int>& backupDataSizes,
			std::vector<char>& sendData,
			std::vector<char>& recvData
	);
private:
	std::vector<MPI_Request> _exchangeSizesRequests;
	std::vector<MPI_Request> _exchangeSnapshotsRequests;
	int const _numProcs;
	int const _rank;
};

#endif /* SRC_PARALLEL_RESILIENCECOMM_H_ */
