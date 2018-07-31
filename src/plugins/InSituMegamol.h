/*
 * InSituMegamol.h
 *
 *  Created on: 11 Jul 2018
 *      Author: tchipevn / Oliver Fernandes
 */

///
/// \file InSituMegamol.h
/// Insitu Megamol Plugin Header. See the InSituMegamol class description for a manual on how to use the plugin
///

#ifndef SRC_PLUGINS_INSITUMEGAMOL_H_
#define SRC_PLUGINS_INSITUMEGAMOL_H_

#include "PluginBase.h"
#include "molecules/MoleculeForwardDeclaration.h"

#include <chrono>
#include <set>
#include <sstream>
#include <iomanip>
#include <vector>
#include <memory>
#include <errno.h>

class Snapshot;

#ifdef INSITU
#include "zmq.h"

//typedef void to easily identify the ZMQ components (which are all essentially pointers to some stuff)
typedef void ZmqContext;
typedef void ZmqRequest;
typedef void ZmqPublish;

class InSituMegamol: public PluginBase {
	/**
	 * @brief Sync status return 
	 * 
	 * 
	 */
	enum ISM_SYNC_STATUS {
		ISM_SYNC_SUCCESS = 0,           // 0
		ISM_SYNC_SYNCHRONIZING,         // 1
		ISM_SYNC_REPLY_BUFFER_OVERFLOW, // 2
		ISM_SYNC_TIMEOUT,               // 3
		ISM_SYNC_SEND_FAILED,           // 4
		ISM_SYNC_UNKNOWN_ERROR          // 5
	};
	/**
	 * @brief Manages the ZMQ calls to Megamol
	 * 
	 * Used to abstract the messages sent to Megamol. Also manages the resources
	 * required for the communication.
	 */
	class ZmqManager {
	public:
		ZmqManager();
		~ZmqManager();
		ZmqManager(ZmqManager const& rhs); //not a real copy: new resources get allocated. should probably share context...

		/**
		 * @brief First 'handshake' communication with Megamol. If this fails (returns false), deactivate the plugin
		 */
		bool performHandshake(void);
		/**
		 * @brief This is called after the data in the shared memory file has been updated.
		 */
		bool triggerUpdate(std::string fname);
		/**
		 * @brief Return a string with the ZMQ version number
		 */
		std::string getZmqVersion(void) const;
		/**
		 * @brief Connect the requester object to the connection name contained in XML input
		 */
		void setConnection(std::string);
		/**
		 * @brief Generate name internal module names strings tied to rank, so Megamol can discriminate
		 */
		void setModuleNames(int rank);
		/**
		 * @brief Sets the maximum size of the buffer for Megamol's reply in zmq_recv call.
		 */
		void setReplyBufferSize(int replyBufferSize) {
			_replyBufferSize = replyBufferSize;
			_replyBuffer.clear();
			_replyBuffer.resize(_replyBufferSize, 0);
		}
		/**
		 * @brief Pass on XML input
		 */
		void setSyncTimeout(int syncTimeout) {
			_syncTimeout = syncTimeout;
		}
	private:
		ZmqManager& operator=(ZmqManager const& rhs); //definitely disallow this guy
		/**
		 * @brief Evaluates replies from Megamol during handshake
		 */
		InSituMegamol::ISM_SYNC_STATUS _synchronizeMegamol(void);
		/**
		 * @brief Waits a sec, returns either ISM_SYNC_TIMEOUT or ISM_SYNC_SYNCHRONIZING
		 */
		InSituMegamol::ISM_SYNC_STATUS _timeoutCheck(void) const;
		/**
		 * @brief Wraps the zmq_send. Wrapper needed to handle the internal send/recv matching.
		 */
		int _send(std::string msg, int blockPolicy);
		/**
		 * @brief Wraps the zmq_recv. Wrapper needed to handle the internal send/recv matching.
		 */
		int _recv(int blockPolicy);

		int _replyBufferSize;
		int _syncTimeout;
		int _sendCount;
		std::vector<char> _replyBuffer;
		std::stringstream _datTag; 
		std::stringstream _geoTag;
	 	std::unique_ptr<ZmqContext, decltype(&zmq_ctx_destroy)> _context;
		std::unique_ptr<ZmqRequest, decltype(&zmq_close)> _requester;
		std::unique_ptr<ZmqPublish, decltype(&zmq_close)> _publisher;
	};
public:
	InSituMegamol(); 
	virtual ~InSituMegamol() {};

	void init(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain
	);

    void readXML(XMLfileUnits& xmlconfig);

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

    void endStep(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            Domain* domain, unsigned long simstep
	);
    
	void finish(ParticleContainer* particleContainer,
            DomainDecompBase* domainDecomp, Domain* domain
	) {}

    std::string getPluginName() {
    	return std::string("InSituMegamol");
    }

	static PluginBase* createInstance() { return new InSituMegamol(); }

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
	Snapshot _snapshot; // make an std::vector eventually
private:
	// XML settings
	int _snapshotInterval;
	int _replyBufferSize;
	int _syncTimeout;
	std::string _connectionName;

	//gather data
	std::vector<char> _generateMmpldSeekTable(std::vector< std::vector<char> >& dataLists);
	std::vector<char> _buildMmpldDataList(ParticleContainer* particleContainer);

	// serialize all
	void _resetMmpldBuffer(void);
	void _addMmpldHeader(float* bbox, float simTime);
	void _addMmpldSeekTable(std::vector<char> seekTable);
	void _addMmpldFrame(std::vector< std::vector<char> > dataLists);
	std::string _writeMmpldBuffer(int rank, unsigned long simstep);

	std::vector<char> _mmpldBuffer;
	std::vector<char>::iterator _mmpldSize;
	ZmqManager _zmqManager;
	bool _isEnabled;
};
#else
class InSituMegamol: public PluginBase {
public:
	InSituMegamol(); 
	virtual ~InSituMegamol() {};
	void init(
			ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain
	) {}

    void readXML(XMLfileUnits& xmlconfig) {};

	void beforeEventNewTimestep(
			ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp,
			unsigned long simstep
	) {}

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
            Domain* domain, unsigned long simstep
	) {}
    
	void finish(
			ParticleContainer* particleContainer,
            DomainDecompBase* domainDecomp, Domain* domain
	) {}

    std::string getPluginName() {
    	return std::string("InSituMegamol");
    }

	static PluginBase* createInstance() { return new InSituMegamol(); }
};
#endif // INSITU
#endif /* SRC_PLUGINS_REDUNDANCYRESILIENCE_H_ */
