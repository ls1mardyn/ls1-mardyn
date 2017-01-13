#include "DomainDecomposition.h"

#include "Domain.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/ParticleData.h"
#include "utils/xmlfileUnits.h"
#include "utils/Logger.h"
#include "parallel/NeighbourCommunicationScheme.h"
#include "parallel/HaloRegion.h"

using Log::global_log;
using namespace std;

DomainDecomposition::DomainDecomposition() :
		DomainDecompMPIBase() {
	int period[DIMgeom]; // 1(true) when using periodic boundary conditions in the corresponding dimension
	int reorder; // 1(true) if the ranking may be reordered by MPI_Cart_create

	// We create a torus topology, so all boundary conditions are periodic
	for (int d = 0; d < DIMgeom; d++)
		period[d] = 1;

	// Allow reordering of process ranks
	reorder = 1;

	for (int i = 0; i < DIMgeom; i++) {
		_gridSize[i] = 0;
	}
	MPI_CHECK(MPI_Dims_create( _numProcs, DIMgeom, (int *) &_gridSize ));

	// Create the communicator
	MPI_CHECK(MPI_Cart_create(MPI_COMM_WORLD, DIMgeom, _gridSize, period, reorder, &_comm));
	global_log->info() << "MPI grid dimensions: " << _gridSize[0] << ", " << _gridSize[1] << ", " << _gridSize[2]
			<< endl;

	// introduce coordinates
	MPI_CHECK(MPI_Comm_rank(_comm, &_rank));
	MPI_CHECK(MPI_Cart_coords(_comm, _rank, DIMgeom, _coords));
	global_log->info() << "MPI coordinate of current process: " << _coords[0] << ", " << _coords[1] << ", "
			<< _coords[2] << endl;
}

DomainDecomposition::~DomainDecomposition() {
	MPI_Comm_free(&_comm);
}

void DomainDecomposition::initCommunicationPartners(double cutoffRadius, Domain * domain) {
	for (int d = 0; d < DIMgeom; ++d) {
		_neighbourCommunicationScheme->setCoverWholeDomain(d, _gridSize[d] == 1);
	}
	_neighbourCommunicationScheme->initCommunicationPartners(cutoffRadius, domain, this);
}

void DomainDecomposition::prepareNonBlockingStage(bool /*forceRebalancing*/, ParticleContainer* moleculeContainer,
		Domain* domain, unsigned int stageNumber) {
	DomainDecompMPIBase::prepareNonBlockingStageImpl(moleculeContainer, domain, stageNumber, LEAVING_AND_HALO_COPIES);
}

void DomainDecomposition::finishNonBlockingStage(bool /*forceRebalancing*/, ParticleContainer* moleculeContainer,
		Domain* domain, unsigned int stageNumber) {
	DomainDecompMPIBase::finishNonBlockingStageImpl(moleculeContainer, domain, stageNumber, LEAVING_AND_HALO_COPIES);
}

bool DomainDecomposition::queryBalanceAndExchangeNonBlocking(bool /*forceRebalancing*/,
		ParticleContainer* /*moleculeContainer*/, Domain* /*domain*/) {
	return true;
}

void DomainDecomposition::balanceAndExchange(bool /*forceRebalancing*/, ParticleContainer* moleculeContainer,
		Domain* domain) {
	DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_AND_HALO_COPIES);
}

void DomainDecomposition::readXML(XMLfileUnits& xmlconfig) {
	/* TODO: Maybe add decomposition dimensions, default auto. */
	DomainDecompMPIBase::readXML(xmlconfig);
}

bool DomainDecomposition::procOwnsPos(double x, double y, double z, Domain* domain) {
	if (x < getBoundingBoxMin(0, domain) || x >= getBoundingBoxMax(0, domain))
		return false;
	else if (y < getBoundingBoxMin(1, domain) || y >= getBoundingBoxMax(1, domain))
		return false;
	else if (z < getBoundingBoxMin(2, domain) || z >= getBoundingBoxMax(2, domain))
		return false;
	else
		return true;
}

double DomainDecomposition::getBoundingBoxMin(int dimension, Domain* domain) {
	return _coords[dimension] * domain->getGlobalLength(dimension) / _gridSize[dimension];
}

double DomainDecomposition::getBoundingBoxMax(int dimension, Domain* domain) {
	return (_coords[dimension] + 1) * domain->getGlobalLength(dimension) / _gridSize[dimension];
}

void DomainDecomposition::printDecomp(string filename, Domain* domain) {

	if (_rank == 0) {
		ofstream povcfgstrm(filename.c_str());
		povcfgstrm << "size " << domain->getGlobalLength(0) << " " << domain->getGlobalLength(1) << " "
				<< domain->getGlobalLength(2) << endl;
		povcfgstrm << "cells " << _gridSize[0] << " " << _gridSize[1] << " " << _gridSize[2] << endl;
		povcfgstrm << "procs " << _numProcs << endl;
		povcfgstrm << "data DomainDecomp" << endl;
		povcfgstrm.close();
	}

	for (int process = 0; process < _numProcs; process++) {
		if (_rank == process) {
			ofstream povcfgstrm(filename.c_str(), ios::app);
			povcfgstrm << _coords[2] * _gridSize[0] * _gridSize[1] + _coords[1] * _gridSize[0] + _coords[0] << " "
					<< _rank << endl;
			povcfgstrm.close();
		}
		barrier();
	}
}

std::vector<int> DomainDecomposition::getNeighbourRanks() {
#if defined(ENABLE_MPI)
	std::vector<int> neighbours;
	if (_numProcs == 1) {
		for (int i = 0; i < 6; i++)
			neighbours.push_back(_rank);
	} else {
		neighbours = _neighbourCommunicationScheme->get3StageNeighbourRanks();
	}
	return neighbours;
#else
	return std::vector<int>(0);
#endif
}

/**
 * The key of this function is that opposite sites are always neighbouring each other in the array (i.e. leftAreaIndex = 0, righAreaIndex = 1, ...)
 *
 **/
std::vector<int> DomainDecomposition::getNeighbourRanksFullShell() {
#if defined(ENABLE_MPI)

	std::vector<int> neighbours(26, -1);
	if (_numProcs == 1) {
		for (int i = 0; i < 26; i++)
			neighbours[i] = _rank;
	} else {
		neighbours = _neighbourCommunicationScheme->getFullShellNeighbourRanks();
	}
	return neighbours;
#else
	return std::vector<int>(0);
#endif
}

std::vector<CommunicationPartner> DomainDecomposition::getNeighboursFromHaloRegion(Domain* domain, const HaloRegion& haloRegion,
		double cutoff) {
//TODO: change this method for support of midpoint rule, half shell, eighth shell, Neutral Territory
// currently only one process per region is possible.
	int rank;
	int regionCoords[DIMgeom];
	for (unsigned int d = 0; d < DIMgeom; d++) {
		regionCoords[d] = _coords[d] + haloRegion.offset[d];
	}
	//TODO: only full shell! (otherwise more neighbours possible)
	MPI_CHECK(MPI_Cart_rank(getCommunicator(), regionCoords, &rank)); //does automatic shift for periodic boundaries
	double haloLow[3];
	double haloHigh[3];
	double boundaryLow[3];
	double boundaryHigh[3];
	double shift[3];

	for (unsigned int d = 0; d < DIMgeom; d++) {
		haloLow[d] = haloRegion.rmin[d];
		haloHigh[d] = haloRegion.rmax[d];
		//TODO: ONLY FULL SHELL!!!
		boundaryLow[d] = haloRegion.rmin[d] - haloRegion.offset[d] * cutoff; //rmin[d] if offset[d]==0
		boundaryHigh[d] = haloRegion.rmax[d] - haloRegion.offset[d] * cutoff; //if offset[d]!=0 : shift by cutoff in negative offset direction
		if (_coords[d] == 0 and haloRegion.offset[d] == -1) {
			shift[d] = domain->getGlobalLength(d);
		} else if (_coords[d] == _gridSize[d] - 1 and haloRegion.offset[d] == 1) {
			shift[d] = -domain->getGlobalLength(d);
		} else{
			shift[d] = 0.;
		}

	}
	// initialize using initializer list - here a vector with one element is created
	return std::vector<CommunicationPartner> {CommunicationPartner(rank, haloLow, haloHigh, boundaryLow, boundaryHigh, shift, haloRegion.offset)};
}
