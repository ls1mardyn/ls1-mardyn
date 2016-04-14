#include "DomainDecomposition.h"

#include "Domain.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/ParticleData.h"
#include "utils/xmlfileUnits.h"
#include "utils/Logger.h"

using Log::global_log;
using namespace std;

DomainDecomposition::DomainDecomposition() : DomainDecompMPIBase() {
	int period[DIM]; // 1(true) when using periodic boundary conditions in the corresponding dimension
	int reorder; // 1(true) if the ranking may be reordered by MPI_Cart_create

	// We create a torus topology, so all boundary conditions are periodic
	for (int d = 0; d < DIM; d++)
		period[d] = 1;

	// Allow reordering of process ranks
	reorder = 1;

	for (int i = 0; i < DIM; i++) {
		_gridSize[i] = 0;
	}
	MPI_CHECK( MPI_Dims_create( _numProcs, DIM, (int *) &_gridSize ) );

	// Create the communicator
	MPI_CHECK( MPI_Cart_create(MPI_COMM_WORLD, DIM, _gridSize, period, reorder, &_comm) );
	global_log->info() << "MPI grid dimensions: " << _gridSize[0]<<", "<<_gridSize[1]<<", "<<_gridSize[2] << endl;

	// introduce coordinates
	MPI_CHECK( MPI_Comm_rank(_comm, &_rank) );
	MPI_CHECK( MPI_Cart_coords(_comm, _rank, DIM, _coords) );
	global_log->info() << "MPI coordinate of current process: " << _coords[0]<<", "<<_coords[1]<<", "<<_coords[2] << endl;

	for (int d = 0; d < DIM; ++d) {
		if(_gridSize[d] == 1) {
			_coversWholeDomain[d] = true;
		} else {
			_coversWholeDomain[d] = false;
		}
	}

	_neighboursInitialized = false;
}

DomainDecomposition::~DomainDecomposition() {
	MPI_Comm_free(&_comm);
}

void DomainDecomposition::initCommunicationPartners(double cutoffRadius, Domain * domain) {

	if(_neighboursInitialized) {
		return;
	}
	_neighboursInitialized = true;

	// corners of the process-specific domain
	double rmin[DIM]; // lower corner
	double rmax[DIM]; // higher corner
	double halo_L[DIM]; // width of the halo strip

	for (int d = 0; d < DIM; d++) {
		rmin[d] = getBoundingBoxMin(d, domain);
		rmax[d] = getBoundingBoxMax(d, domain);

		// TODO: this should be safe, as long as molecules don't start flying around
		// at the speed of one cutoffRadius per timestep
		halo_L[d] = cutoffRadius;

		_neighbours[d].clear();
	}

	int direction;

	for (unsigned short d = 0; d < DIM; d++) {
		if(_coversWholeDomain[d]) {
			// nothing to do;
			continue;
		}

		// set the ranks
		int ranks[2];

		MPI_CHECK( MPI_Cart_shift(_comm, d, 1, &ranks[LOWER], &ranks[HIGHER]) );

		// when moving a particle across a periodic boundary, the molecule position has to change
		// these offset specify for each dimension (x, y and z) and each direction ("left"/lower
		// neighbour and "right"/higher neighbour, how the paritcle coordinates have to be changed.
		// e.g. for dimension x (d=0) and a process on the left boundary of the domain, particles
		// moving to the left get the length of the whole domain added to their x-value
		double offsetLower[DIM];
		double offsetHigher[DIM];
		offsetLower[d] = 0.0;
		offsetHigher[d] = 0.0;

		// process on the left boundary
		if (_coords[d] == 0)
			offsetLower[d] = domain->getGlobalLength(d);
		// process on the right boundary
		if (_coords[d] == _gridSize[d] - 1)
			offsetHigher[d] = -domain->getGlobalLength(d);

		for (direction = LOWER; direction <= HIGHER; direction++) {
			double regToSendLow[DIM];
			double regToSendHigh[DIM];

			// set the regions
			for (int i = 0; i < DIM; i++) {
				regToSendLow[i] = rmin[i] - halo_L[i];
				regToSendHigh[i] = rmax[i] + halo_L[i];
			}

			double haloLow[3];
			double haloHigh[3];
			double boundaryLow[3];
			double boundaryHigh[3];

			switch (direction) {
			case LOWER:
				regToSendHigh[d] = rmin[d] + halo_L[d];
				for (int i = 0; i < DIM; ++i) {
					haloLow[i] = regToSendLow[i];
					if (i == d) {
						haloHigh[i] = boundaryLow[i] = rmin[i];
					} else {
						haloHigh[i] = regToSendHigh[i];
						boundaryLow[i] = regToSendLow[i];
					}
					boundaryHigh[i] = regToSendHigh[i];
				}
				break;
			case HIGHER:
				regToSendLow[d] = rmax[d] - halo_L[d];
				for (int i = 0; i < DIM; ++i) {
					boundaryLow[i] = regToSendLow[i];
					if (i == d) {
						boundaryHigh[i] = haloLow[i] = rmax[i];
					} else {
						boundaryHigh[i] = regToSendHigh[i];
						haloLow[i] = regToSendLow[i];
					}
					haloHigh[i] = regToSendHigh[i];
				}
				break;
			}

			// set the shift
			double shift[3] = {0., 0., 0.};
			if (direction == LOWER)
				shift[d] = offsetLower[d];
			if (direction == HIGHER)
				shift[d] = offsetHigher[d];

			_neighbours[d].push_back(
					CommunicationPartner(ranks[direction], haloLow, haloHigh, boundaryLow, boundaryHigh, shift));
		}
	}
}

void DomainDecomposition::balanceAndExchange(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain) {
	DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_AND_HALO_COPIES);
}

void DomainDecomposition::readXML(XMLfileUnits& xmlconfig) {
	/* no parameters */
	/* TODO: Maybe add decomposition dimensions, default auto. */
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
		povcfgstrm << "size " << domain->getGlobalLength(0) << " " << domain->getGlobalLength(1) << " " << domain->getGlobalLength(2) << endl;
		povcfgstrm << "cells " << _gridSize[0] << " " << _gridSize[1] << " " << _gridSize[2] << endl;
		povcfgstrm << "procs " << _numProcs << endl;
		povcfgstrm << "data DomainDecomp" << endl;
		povcfgstrm.close();
	}

	for (int process = 0; process < _numProcs; process++) {
		if (_rank == process) {
			ofstream povcfgstrm(filename.c_str(), ios::app);
			povcfgstrm << _coords[2] * _gridSize[0] * _gridSize[1] + _coords[1] * _gridSize[0] + _coords[0] << " " << _rank << endl;
			povcfgstrm.close();
		}
		barrier();
	}
}

virtual std::vector<int> DomainDecomposition::getNeighbourRanks(){
	std::vector<int> neighbours;
	for(int d = 0; d < DIM;d++){
		for(int n = 0; n < 2; n++){
			neighbours.push_back(_neighbours[d][n].getRank());
		}
	}
	return neighbours;
}


