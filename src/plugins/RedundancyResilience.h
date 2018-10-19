/*
 * RedundancyResilience.h
 *
 *  Created on: 11 Jul 2018
 *      Author: tchipevn / Oliver Fernandes
 */

///
/// \file RedundancyResilience.h
/// Redundancy Resilience Plugin Header. See the RedundancyResilience class description for a manual on how to use the plugin
///

#ifndef SRC_PLUGINS_REDUNDANCYRESILIENCE_H_
#define SRC_PLUGINS_REDUNDANCYRESILIENCE_H_

#include "PluginBase.h"
#include "molecules/MoleculeForwardDeclaration.h"
#include "utils/mardyn_assert.h" 

#ifdef ENABLE_MPI
#include <set>
#include <vector>
#include <memory>
#include <map>
#include <climits> /* UINT64_MAX */

#include <stddef.h>
#include <mpi.h>

class Snapshot;
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
#define RR_INTS_PER_RANK 4
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
	 * The tags are used to identify the messages in send/recv calls.
	 * @param[in] backupInfo Contains the rank configuration. Set up in RedundancyResilience::_determineBackups, check there for layout.
	 * @param[in] numberOfBackups how many backups are made per rank
	 * @param[in] sizePerRank Number of information elements per rank (numberOfBackups*number of arrays, atm 4)
	 * @param[out] backing A list of ranks being backed by the current rank
	 * @param[out] backedBy A list of ranks backing the current rank
	 * @param[out] backingTags The tags associated with the communication of the ranks in backing (used in recv)
	 * @param[out] backedByTags The tags associated with the communication of the ranks in backedBy (used in send)
	 */
	int scatterBackupInfo(
			std::vector<int>& backupInfo, 
			int const numberOfBackups,
			size_t const sizePerRank,
	        std::vector<int>& backing, 
	        std::vector<int>& backedBy, 
	        std::vector<int>& backingTags, 
	        std::vector<int>& backedByTags
	);
	/**
	 * Send snapshot size from caller to ranks backing it.
	 * Called by each process after setting up the new snapshot.
	 * Needs to be called before calling RedundancyResilience::exchangeSnapshots. If this was done already
	 * should probably be checked in that method. But it isn't.
	 * @param[in] backing A list of ranks being backed by the current rank
	 * @param[in] backedBy A list of ranks the current rank is backing up
	 * @param[in] backingTags The tags associated with the communication of the ranks in backing (used in recv)
	 * @param[in] backedByTags The tags associated with the communication of the ranks in backedBy (used in send)
	 * @param[in] snapshotSize Size of the local snapshot data in bytes
	 * @param[out] backupDataSizes The individual sizes of the snapshots acquired 
	 */
	int exchangeSnapshotSizes(
			std::vector<int>& backing,
			std::vector<int>& backedBy,
			std::vector<int>& backingTags,
			std::vector<int>& backedByTags,
			size_t const snapshotSize, 
			std::vector<int>& backupDataSizes
	);
	/**
	 * Send snapshot meta data from caller to ranks backing it.
	 * Called by each process after trading the backup snapshot sizes.
	 * This will perform the actual sending of snapshot data for each of the backup ranks.
	 * @param[in] backing A list of ranks being backed by the current rank
	 * @param[in] backedBy A list of ranks the current rank is backing up
	 * @param[in] backingTags The tags associated with the communication of the ranks in backing (used in recv)
	 * @param[in] backedByTags The tags associated with the communication of the ranks in backedBy (used in send)
	 * @param[in] backupDataSizes The individual sizes of the snapshots acquired 
	 * @param[in] sendData The local snapshot data as char vector
	 * @param[out] recvData The complete snapshot datas of the ranks being backed
	 */
	int exchangeSnapshots(
			std::vector<int>& backing,
			std::vector<int>& backedBy,
			std::vector<int>& backingTags,
			std::vector<int>& backedByTags,
			std::vector<int>& backupDataSizes,
			std::vector<char>& sendData,
			std::vector<char>& recvData
	);
private:
	int const _numProcs;
	int const _rank;
};

class RedundancyResilience: public PluginBase {
public:
	RedundancyResilience() {}; 
	virtual ~RedundancyResilience() {}
	/**
	 * @brief Method to initialize the resilience setup.
	 * 
	 * This callback does everything necessary to be able to use the backup functionality.
	 * In particular, it ensures all ranks know where they should take data to backup, and 
	 * where they should send data to be backed up, and creates RedundancyResilience::_comm 
	 * for abstracting MPI calls.
	 * It also does some simple sanity checking of configuration values.
	 * @param[in] particleContainer Particle data
	 * @param[in] domainDecomp Data about the domain decomposition
	 * @param[in] domain Data about the domain
	 */
	void init(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain);

    void readXML(XMLfileUnits& xmlconfig);

    /**
     * @brief restarting takes place here
     */
	void beforeEventNewTimestep(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	);

    void beforeForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    ) {}

    void afterForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    ) {}
	/**
	 * @brief Check if new snapshot is requested here.
	 * 
	 * Use this plugin callback to check if we need to do a backup in the current time step.
	 * The decision is based on RedundancyResilience::_backupInterval, set in the XML 
	 * configuration of this plugin.
	 * @param[in] particleContainer Contains the rank's particle data
	 * @param[in] domainDecomp The domain decomposition
	 * @param[in] domain Data for the domain
	 * @param[in] simstep Current simulation time step that is being processed
	*/
    void endStep(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            Domain* domain, unsigned long simstep);

    void finish(ParticleContainer* particleContainer,
                              DomainDecompBase* domainDecomp, Domain* domain){}

    std::string getPluginName() {
    	return std::string("RedundancyResilience");
    }

	static PluginBase* createInstance() { return new RedundancyResilience(); }

	class Snapshot {
	public:
		void addMolecule(const Molecule& m) {
			_molecules.push_back(m);
		}

		double getCurrentTime() const {
			return _currentTime;
		}

		void setCurrentTime(double currentTime) {
			_currentTime = currentTime;
		}

		const std::array<double, 3>& getBoxDims() const {
			return _boxDims;
		}

		void setBoxDims(const std::array<double, 3>& boxDims) {
			_boxDims = boxDims;
		}

		unsigned long getGlobalNumberOfMolecules() const {
			return _globalNumberOfMolecules;
		}

		void setGlobalNumberOfMolecules(unsigned long globalNumberOfMolecules) {
			_globalNumberOfMolecules = globalNumberOfMolecules;
		}

		double getTemperature() const {
			return _temperature;
		}

		void setTemperature(double temperature) {
			_temperature = temperature;
		}

		int getRank() const {
			return _rank;
		}

		void setRank(int rank) {
			_rank = rank;
		}

		const std::vector<Molecule>& getMolecules() const {
			return _molecules;
		}

		void clearMolecules() {
			_molecules.clear();
		}

	private:
		// snapshot data
		std::vector<Molecule> _molecules;        ///< the molecule data should be backed up in here
		double _currentTime;                     ///< the time step this snap shot was made should be stored here
		int _rank;                     
		// the following fields are maybe unnecessary, but leaving them here now for consistency to written headers in file-checkpoints
		unsigned long _globalNumberOfMolecules;
		double _temperature;                     // maybe not necessary; for consistency to currently written headers
		std::array<double, 3> _boxDims;          // maybe not necessary; for consistency to currently written headers
	};

protected:
	/**
	 * @brief Method to determine backup pairings.
	 * 
	 * This function fills backupInfo with 2*RedundancyResilience::_numberOfBackups(n)*number of ranks(p) values defining the
	 * backup topology, i.e. which rank backs up which. It is made up of concatenated subarrays p_i, for each subarray you have 2*n
	 * entries. The first n in each p_i are the ranks that are backed up by rank i. The last n are the ranks that are backing up i.
	 * Consider following example with p=4 ranks and n=1 redundant backups requested:
	 * \code{.unparsed}
	 * layout :  Rank 0                        Rank 1                    Rank 2                    Rank 3
	 *           Backing:[1]  BackedBy:[3]     Backing:[2] BackedBy:[0]  Backing:[3] BackedBy:[1]  Backing:[0] BackedBy:[2]
	 * \endcode
	 * The backupInfo array therefore looks like this : [1 3 2 0 3 1 0 2]
	 * This is executed only on the root node, to allow for arbitrary but coordinated backup assignments. The resulting array is
	 * communicated using RedundancyResilience::_comm.
	 * @param[in] domainDecomp Data about the domainDecomposition
	 * @param[out] backupInfo Empty vector used to hold resulting array. Will be resized accordingly.
	 */
	std::vector<int> _determineBackups(DomainDecompBase const* domainDecomp);
private:
	/**
	 * @brief Copy the local simulation data to a snapshot
	 * This function copies all data necessary to restart the process to *this object's _snapshot instance.
	 * This is also the data that is communicated via MPI calls.
	 * @param[in] particleContainer Contains the rank's particle data
	 * @param[in] domainDecomp The domain decomposition
	 * @param[in] domain Data for the domain
	 * @param[in] simstep Current simulation time step that is being processed
	 */
	void _saveLocalSnapshot(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            Domain* domain, unsigned long simstep);
	/**
	 * @brief Converts the snapshot to bytes for sending
	 * This function transforms the local snapshot data to an array of char, to be used for MPI calls.
	 * @return A shared pointer to the byte array containing the converted data
	 */
	std::vector<char> _serializeSnapshot(void) const;
	/**
	 * @brief Stores received data to snapshots
	 * This method transforms the data received as bytes from the ResilienceComm::_exchangeSnapshots() call.
	 * It partitions the data according to the sizes in backupDataSizes, and stores them in RedundancyResilience::_backupSnapshots,
	 * using the RedundancyResilience::_deserializeSnapshot() call.
	 */
	void _storeSnapshots(std::vector<int>& backupDataSizes, std::vector<char>& backupData);
	/**
	 * @brief Converts a byte array into a snapshot
	 * This method transforms the data received as bytes from the ResilienceComm::_exchangeSnapshots() call.
	 * The actual data being worked on is decided by passing the correct iterators. The method assumes a fixed layout of the
	 * serialized data. (atm: sizeof(int) bytes -> int rank, sizeof(double) bytes -> currentTime, n*sizeof(Molecule) -> std::vector<Molecule>)
	 * @param[in] snapshotStart Iterator to the start of the the snapshot's byte data
	 * @param[in] snapshotEnd Iterator to the end of the the snapshot's byte data. snapshotEnd itself is not part of the snapshot
	 * @param[out] snapshot Target for the incoming data
	 */
	std::vector<char>::iterator _deserializeSnapshot(std::vector<char>::iterator const snapshotStart,
                                                     std::vector<char>::iterator const snapshotEnd,
													 Snapshot& snapshot);
	/**
	 * @brief Validate the fake data to be correct for debugging
	 * The fake data for a snapshot is just an array of ints going from 0 to the snapshots rank
	 * In actual scenarios, this would be particle data.
	 */
	bool _validateFakeData(int const rank, std::vector<char>& fakeData);
	
	std::unique_ptr<ResilienceComm> _comm;       ///< store the communication handling object
	std::vector<int> _backing;                   ///< contains the rank ids this rank does redundancy backups for
	std::vector<int> _backingTags;               ///< contains the associated comm tags for the ranks in _backing
	std::vector<int> _backedBy;                  ///< contains the ranks of processes which create backups for this rank
	std::vector<int> _backedByTags;              ///< contains the associated comm tags for the ranks in _backedBy
	std::vector<Snapshot> _backupSnapshots;      ///< contains the data of the backed up ranks

	Snapshot _snapshot; // make an std::vector eventually
	unsigned long _backupInterval;
	int _numberOfBackups;                        // number of ranks to backup
	int _sizePerRank;
};

#else /* ENABLE_MPI */

class Snapshot;

class RedundancyResilience: public PluginBase {
public:
	RedundancyResilience() {};
 	virtual ~RedundancyResilience() {}
 	void init(ParticleContainer* particleContainer,
 			DomainDecompBase* domainDecomp, Domain* domain) {}

    void readXML(XMLfileUnits& xmlconfig) {}
	void beforeEventNewTimestep(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) {}
	MPI_Buffer_attach()

    void beforeForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    ) {}

    void afterForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    ) {}
    void endStep(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            Domain* domain, unsigned long simstep) {};

    void finish(ParticleContainer* particleContainer,
                              DomainDecompBase* domainDecomp, Domain* domain
	) {}

    std::string getPluginName() {
    	return std::string("RedundancyResilience");
    }

	static PluginBase* createInstance() { 
		return new RedundancyResilience(); 
	}
};

#endif /* ENABLE_MPI */

#endif /* SRC_PLUGINS_REDUNDANCYRESILIENCE_H_ */
