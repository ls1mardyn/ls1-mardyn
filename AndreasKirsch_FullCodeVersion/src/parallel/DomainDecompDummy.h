#ifndef DOMAINDECOMPDUMMY_H_
#define DOMAINDECOMPDUMMY_H_

#include "parallel/DomainDecompBase.h"
#include "parallel/CollectiveCommDummy.h"

#include <iostream>

//! @brief implement the %domain decomposition for a single processor
//! @author Martin Buchholz
//!
//! As explained in the comment to DomainDecompBase, there is only one
//! code for the sequential and the parallel version. This class
//! accomplished the decomposition for the sequential version (1 process)
//!
//! While this class doesn't has to implement communication between processes,
//! there are still some other things to do, especially the handling of the
//! boundary region.
class DomainDecompDummy : public DomainDecompBase {
public:
	//! The constructor has nothing to do
	DomainDecompDummy();

	//! The destructor has nothing to do
	~DomainDecompDummy();

	//! molecules which aren't in the domain of their process any
	//! are moved to the opposite side of the domain (periodic boundary).
	//! Additionally, the molecules from the boundary region are duplicated
	//! and copied into the corresponding halo-regions.
	void exchangeMolecules(ParticleContainer* moleculeContainer, const std::vector<Component>& components, Domain* domain);

	//! @brief in the sequential version, no balancing is necessary --> calls exchangeMolecules
	void balanceAndExchange(bool balance, ParticleContainer* moleculeContainer, const std::vector<Component>& components, Domain* domain);

	//! @brief returns true
	//!
	//! It is assumed that in the sequential version, the only process possesses all particles
	//! Therefore, even for a position outside of the domain, true is returned
	bool procOwnsPos(double x, double y, double z, Domain* domain) {
		return true;
	}

	//! @brief returns 0.0
	//!
	//! In the sequential version, the only process covers the whole domain, so
	//! there is not distance between any point and the region of the process
	double guaranteedDistance(double x, double y, double z, Domain* domain) {
		return 0.0;
	}

	// documentation see father class (DomainDecompBase.h)
	unsigned long countMolecules(ParticleContainer* moleculeContainer, std::vector<unsigned long> &compCount);

	// documentation see father class (DomainDecompBase.h)
	double getBoundingBoxMin(int dimension, Domain* domain);

	// documentation see father class (DomainDecompBase.h)
	double getBoundingBoxMax(int dimension, Domain* domain);

	//! @brief In the sequential mode, there is no decomposition and therefore nothing to be printed
	void printDecomp(std::string filename, Domain* domain) {
		std::cerr << "printDecomp useless in serial mode" << std::endl;
	}

	//! @brief opends the file(append), loops over all molecules and writes each into the file
	void writeMoleculesToFile(std::string filename, ParticleContainer* moleculeContainer);

	//! @brief There is only one process, so this method always returns 0
	int getRank() {
		return 0;
	}

	//! @brief There is only one process, so this method always returns 1
	int getNumProcs() {
		return 1;
	}

	//! @brief for the sequential version, the processor name is not that important
	//! @todo implement this
	const char* getProcessorName() const;

	//! @brief one process doesn't need synchronisation, so nothing is done here
	void barrier() {
	}

	// documentation see father class (DomainDecompBase.h)
	double getTime();

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
	// the documentation of the class CollectiveCommDummy and of the
	// father class of this class (DomainDecompBase.h)
	//##################################################################
	void collCommInit(int numValues) {
		_collComm.init(numValues);
	}

	void collCommFinalize() {
		_collComm.finalize();
	}

	void collCommAppendInt(int intValue) {
		_collComm.appendInt(intValue);
	}

	void collCommAppendUnsLong(unsigned long unsLongValue) {
		_collComm.appendUnsLong(unsLongValue);
	}

	void collCommAppendFloat(float floatValue) {
		_collComm.appendFloat(floatValue);
	}

	void collCommAppendDouble(double doubleValue) {
		_collComm.appendDouble(doubleValue);
	}

	void collCommAppendLongDouble(long double longDoubleValue) {
		_collComm.appendLongDouble(longDoubleValue);
	}

	int collCommGetInt() {
		return _collComm.getInt();
	}

	unsigned long collCommGetUnsLong() {
		return _collComm.getUnsLong();
	}

	float collCommGetFloat() {
		return _collComm.getInt();
	}

	double collCommGetDouble() {
		return _collComm.getDouble();
	}

	long double collCommGetLongDouble() {
		return _collComm.getLongDouble();
	}

	void collCommAllreduceSum() {
	}

	void collCommBroadcast() {
	}

private:
	//! Dummy variable for sequential "collective" communication, basically only
	//! needed to store values and read them again.
	CollectiveCommDummy _collComm;

};

#endif /*DOMAINDECOMPDUMMY_H_*/
