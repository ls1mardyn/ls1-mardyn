#ifndef DOMAINDECOMPOSITION_H_
#define DOMAINDECOMPOSITION_H_

#include <mpi.h>

#include "parallel/CollectiveCommunication.h"
#include "parallel/DomainDecompBase.h"
#include "parallel/ParticleData.h"

#define DIM 3

#define LOWER  0
#define HIGHER 1

/** @brief Basic domain decomposition based parallelisation, dividing the
 * domain into \#procs equal sized cuboids
 *
 * In a domain decomposition, each process gets part of the spacial domain.
 * In this implementation, the whole domain has to be a cuboid which is
 * decomposed into several, equally sized smaller cuboids.
 * At the boundary, each process needs molecules from neighbouring domains to
 * be able to calculate the forces on the own molecules.
 * Molecules are moving across the boundaries of local domains. So methods are
 * implemented to transfer those molecules.
 *
 * @cite Griebel-2007
 */
class DomainDecomposition : public DomainDecompBase {
public:
	//! @brief The constructor has to determine the own rank and the number of neighbours and
	//!        sets up the topology
	DomainDecomposition();

	// documentation see father class (DomainDecompBase.h)
	virtual ~DomainDecomposition();

	/** @brief Read in XML configuration for KDDecomposition and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <parallelisation type="DomainDecomposition">
	     <!-- no parameters -->
	   </parallelisation>
	   \endcode
	 */
	virtual void readXML(XMLfileUnits& xmlconfig);

	//! @brief exchange molecules between processes
	//!
	//! molecules which aren't in the domain of their process any
	//! more are transferred to their neighbours. Additionally, the
	//! molecules for the halo-region are transferred. To reduce the number
	//! of neighbours a single process has to communicate with, particles
	//! that i.e. have to be moved to the lower right neighbour are
	//! moved to the right neighbour first and then from the right neighbour
	//! to the lower neighbour.
	//! @param moleculeContainer needed to get those molecules which have to be exchanged
	//! @param domain is e.g. needed to get the size of the local domain
	virtual void exchangeMolecules(ParticleContainer* moleculeContainer, Domain* domain);

	//! @brief this decompositin does no balancing, it just exchanges the particles
	//!
	//! This domain decomposition devides the domain into equally sized smaller cuboids, therefore
	//! a balancing of the load is not possible. The method only has to ensure that the particles
	//! between the processes are exchanged, therefore exchangeMolecules is called.
	//! @param balance has no influence in this implementation
	//! @param moleculeContainer needed for calculating load and to get the particles
	//! @param domain is e.g. needed to get the size of the local domain
	virtual void balanceAndExchange(bool balance, ParticleContainer* moleculeContainer, Domain* domain);

	// documentation see father class (DomainDecompBase.h)
	bool procOwnsPos(double x, double y, double z, Domain* domain);

	// documentation see father class (DomainDecompBase.h)
	virtual double getBoundingBoxMin(int dimension, Domain* domain);

	// documentation see father class (DomainDecompBase.h)
	virtual double getBoundingBoxMax(int dimension, Domain* domain);

	//! @brief writes information about the current decomposition into the given file
	//!
	//! This decomposition first writes a small header with the following lines:
	//! - "size", followed by three double values for the size of the domain
	//! - "cells", followed by three int values for the number of cells in each direction.
	//!            Each cell here stands for the domain of one process, so the product of
	//!            of those three int values is the number of processes
	//! - "procs", follwed by one int value (number of processes)
	//! - "data DomainDecomp": just this text to state that now the header is finished and data follows
	//! - one line per proc, containing the "cellid" using a x-y-z-ordering and the rank of the proc
	//!
	//! An example file for 8 procs could look like this:
	//!
	//! |size 62.0 62.0 62.0 \n
	//!  cells 2 2 2 \n
	//!  procs 8 \n
	//!  data DomainDecomp \n
	//!  1 1 \n
	//!  2 2 \n
	//!  3 3 \n
	//!  4 4 \n
	//!  5 5 \n
	//!  6 6 \n
	//!  7 7 \n
	//!  8 8
	//! @param filename name of the file into which the data will be written
	//! @param domain e.g. needed to get the bounding boxes
	virtual void printDecomp(std::string filename, Domain* domain);

	// documentation see father class (DomainDecompBase.h)
	void barrier() { MPI_CHECK( MPI_Barrier(_comm) ); }

	//! @brief returns total number of molecules
	unsigned Ndistribution(unsigned localN, float* minrnd, float* maxrnd);

	//! @brief checks identity of random number generators
	void assertIntIdentity(int IX);
	void assertDisjunctivity(TMoleculeContainer* mm);

	//##################################################################
	// The following methods with prefix "collComm" are all used
	// in the context of collective communication. Each of the methods
	// basically has to call the corresponding method from the class
	// CollectiveCommunication (or CollectiveCommDummy in the sequential
	// case). To get information about how to use this methods, read
	// the documentation of the class CollectiveCommunication and of the
	// father class of this class (DomainDecompBase.h)
	//##################################################################
	void collCommInit(int numValues) {
		_collCommunication.init(_comm, numValues);
	}

	void collCommFinalize() {
		_collCommunication.finalize();
	}

	void collCommAppendInt(int intValue) {
		_collCommunication.appendInt(intValue);
	}

	void collCommAppendUnsLong(unsigned long unsLongValue) {
		_collCommunication.appendUnsLong(unsLongValue);
	}

	void collCommAppendFloat(float floatValue) {
		_collCommunication.appendFloat(floatValue);
	}

	void collCommAppendDouble(double doubleValue) {
		_collCommunication.appendDouble(doubleValue);
	}

	void collCommAppendLongDouble(long double longDoubleValue) {
		_collCommunication.appendLongDouble(longDoubleValue);
	}

	int collCommGetInt() {
		return _collCommunication.getInt();
	}

	unsigned long collCommGetUnsLong() {
		return _collCommunication.getUnsLong();
	}

	float collCommGetFloat() {
		return _collCommunication.getFloat();
	}

	double collCommGetDouble() {
		return _collCommunication.getDouble();
	}

	long double collCommGetLongDouble() {
		return _collCommunication.getLongDouble();
	}

	void collCommAllreduceSum() {
		_collCommunication.allreduceSum();
	}

	void collCommBroadcast(int root = 0) {
		_collCommunication.broadcast(root);
	}

protected:
	MPI_Datatype _mpi_Particle_data;

	//! new topology after initializing the torus
	MPI_Comm _comm;

	//! flag, which tells whether a processor covers the whole domain along a dimension
	//! if true, we will use the methods provided by the base class for handling the
	//! respective dimension, instead of packing and unpacking messages to self
	bool _coversWholeDomain[DIM];

	struct _CommunicationPartner {
		int _rank;
		double _regionLow[3], _regionHigh[3];
		MPI_Request _sendRequest, _recvRequest;
		MPI_Status _sendStatus, _recvStatus;
		std::vector<ParticleData> _sendBuf, _recvBuf;
		double _shift; //! for periodic boundaries
	};

private:
	void initCommunicationPartners(ParticleContainer * moleculeContainer, Domain * domain);

	//! Number of processes in each dimension (i.e. 2 for 8 processes)
	int _gridSize[DIM];

	//! Grid coordinates of process
	int _coords[DIM];

	//! Array of neighbour ranks.
	//! The first array index specifies the coordinate index,
	//! the second one the direction. For the later use the predefined LOWER and HIGHER macros.
	_CommunicationPartner _partners[DIM][2];

	//! variable used for different kinds of collective operations
	CollectiveCommunication _collCommunication;
};

#endif /* DOMAINDECOMPOSITION_H_ */
