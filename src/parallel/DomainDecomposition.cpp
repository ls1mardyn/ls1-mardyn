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
	double halo_width[DIM]; // width of the halo strip

	for (int d = 0; d < DIM; d++) {
		rmin[d] = getBoundingBoxMin(d, domain);
		rmax[d] = getBoundingBoxMax(d, domain);

		// TODO: this should be safe, as long as molecules don't start flying around
		// at the speed of one cutoffRadius per time step
		halo_width[d] = cutoffRadius;

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

		// When moving a particle across a periodic boundary, the molecule position has to change.
		// These offsets specify for each dimension (x, y and z) and each direction ("left"/lower
		// neighbor and "right"/higher neighbor, how the particle coordinates have to be changed.
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
				regToSendLow[i] = rmin[i] - halo_width[i];
				regToSendHigh[i] = rmax[i] + halo_width[i];
			}

			double haloLow[3];
			double haloHigh[3];
			double boundaryLow[3];
			double boundaryHigh[3];

			switch (direction) {
			case LOWER:
				regToSendHigh[d] = rmin[d] + halo_width[d];
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
				regToSendLow[d] = rmax[d] - halo_width[d];
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

int DomainDecomposition::getNonBlockingStageCount(){
	return 3;  // three stages: first x, then y, then z
}

void DomainDecomposition::prepareNonBlockingStage(bool /*forceRebalancing*/,
		ParticleContainer* moleculeContainer, Domain* domain,
		unsigned int stageNumber) {
	DomainDecompMPIBase::prepareNonBlockingStageImpl(moleculeContainer, domain, stageNumber, LEAVING_AND_HALO_COPIES);
}

void DomainDecomposition::finishNonBlockingStage(bool /*forceRebalancing*/,
		ParticleContainer* moleculeContainer, Domain* domain,
		unsigned int stageNumber) {
	DomainDecompMPIBase::finishNonBlockingStageImpl(moleculeContainer, domain, stageNumber, LEAVING_AND_HALO_COPIES);
}

bool DomainDecomposition::queryBalanceAndExchangeNonBlocking(bool /*forceRebalancing*/,
		ParticleContainer* /*moleculeContainer*/, Domain* /*domain*/) {
	return true;
}

void DomainDecomposition::balanceAndExchange(bool /*forceRebalancing*/, ParticleContainer* moleculeContainer, Domain* domain) {
	DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_AND_HALO_COPIES);
}

void DomainDecomposition::readXML(XMLfileUnits& /*xmlconfig*/) {
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

std::vector<int> DomainDecomposition::getNeighbourRanks(){
#if defined(ENABLE_MPI)
	int numProcs;
	MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
	std::vector<int> neighbours;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if(numProcs == 1){

		for(int i = 0; i<6; i++)
			neighbours.push_back(rank);
	}
	else{
		for(int d = 0; d < DIM;d++){
			for(int n = 0; n < 2; n++){
				if(_coversWholeDomain[d]){
					neighbours.push_back(rank);
				}
				else{
					neighbours.push_back(_neighbours[d][n].getRank());
				}
			}
		}
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
std::vector<int> DomainDecomposition::getNeighbourRanksFullShell(){
#if defined(ENABLE_MPI)
	int numProcs;
	MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
	int myRank;
	MPI_Comm_rank(MPI_COMM_WORLD,&myRank);

	std::vector<int> neighbours(26,-1);
	if(numProcs == 1){

		for(int i = 0; i<26; i++)
			neighbours[i] = myRank;
	}
	else{
		for(int i = 0; i < 2*DIM;i++){
			if(_coversWholeDomain[i/2]){
				neighbours[i] = myRank;
			}
			else{
				neighbours[i] = (_neighbours[i/2][i%2].getRank());
			}
//			std::cout << neighbours[i] << "\n";
		}
		std::vector<int> offsets(6,0);
		//calculate the rank offsets in every dimension in plus and minus direction
		for(int i = 0; i < DIM * 2; i++){
			offsets[i] = neighbours[i]-myRank;
//			std::cout << offsets[i] << " ";
		}
//		std::cout << "\n";

		//calculate remaining 20 neighbours through offsets

		//edges
		//low x direction
		for(int i = 0; i < 4; i++){ //all edges that are adjacent to lower x area (left)
			neighbours[2*i+6] = neighbours[0] + offsets[i+2];
		}
		//higher x direction
		for(int i = 0; i < 4; i++){ //all edges that are adjacent to higher x area (right)
			//always get opposite edges next to each other
			int indexShift = (i%2 == 0)? +1 : -1;
			neighbours[2*i+7] = neighbours[1] + offsets[i+2+indexShift];
		}

		//low y direction
		for(int i = 0; i < 2; i++){ //all edges that are adjacent to lower y  area (bottom) not adjacent to lower x (left)
			neighbours[2*i+14] = neighbours[2] + offsets[i+4];
		}
		//higher y direction
		for(int i = 0; i < 2; i++){ //all edges that are adjacent to higher y area (top) and not adjacent to lower x (right)
			//always get opposite edges next to each other
			int indexShift = (i%2 == 0)? +1 : -1;
			neighbours[2*i+15] = neighbours[3] + offsets[i+4+indexShift];
		}

		//corners
		//lower x direction
		int index = 18;
		for(int i = 0; i < 2; i++){//y offset
			for(int j = 0; j < 2; j++){ // z offset
				neighbours[index] = neighbours[0] + offsets[2+i] + offsets[4+j];
				index = index + 2;
			}
		}
		index = 19;
		for(int i = 1; i >= 0; i--){//y offset
			for(int j = 1; j >= 0; j--){ // z offset
				neighbours[index] = neighbours[1] + offsets[2+i] + offsets[4+j];
				index = index + 2;
			}
		}
	}
	return neighbours;
#else
	return std::vector<int>(0);
#endif
}
