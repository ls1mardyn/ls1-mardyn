#include "DomainDecomposition.h"

#include "Domain.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/xmlfileUnits.h"
#include "utils/Logger.h"
#include "parallel/NeighbourCommunicationScheme.h"
#include "parallel/HaloRegion.h"
#include "ParticleData.h"

using Log::global_log;
using namespace std;

DomainDecomposition::DomainDecomposition() : DomainDecomposition(MPI_COMM_WORLD) {}

DomainDecomposition::DomainDecomposition(MPI_Comm comm) : DomainDecompMPIBase(comm), _gridSize{0}, _coords{0} {
	initMPIGridDims();
}

void DomainDecomposition::initMPIGridDims() {
	mardyn_assert(DIMgeom == 3);
	int period[DIMgeom] = {1, 1, 1}; // 1(true) when using periodic boundary conditions in the corresponding dimension
	int reorder = 1; // 1(true) if the ranking may be reordered by MPI_Cart_create
	{
		auto numProcsGridSize = _gridSize[0] * _gridSize[1] * _gridSize[2];
		if (numProcsGridSize != _numProcs and numProcsGridSize != 0) {
			global_log->error() << "DomainDecomposition: Wrong grid size given!" << std::endl;
			global_log->error() << "\tnumProcs is " << _numProcs << "," << std::endl;
			global_log->error() << "\tbut grid is " << _gridSize[0] << " x " << _gridSize[1] << " x " << _gridSize[2] << std::endl;
			global_log->error() << "\tresulting in " << numProcsGridSize << " subdomains!" << std::endl;
			global_log->error() << "\tplease check your input file!" << std::endl;
			Simulation::exit(2134);
		}
	}

	MPI_CHECK(MPI_Dims_create( _numProcs, DIMgeom, (int *) &_gridSize ));
	MPI_CHECK(MPI_Cart_create(_comm, DIMgeom, _gridSize, period, reorder, &_comm));
	global_log->info() << "MPI grid dimensions: " << _gridSize[0] << ", " << _gridSize[1] << ", " << _gridSize[2] << endl;
	MPI_CHECK(MPI_Comm_rank(_comm, &_rank));
	MPI_CHECK(MPI_Cart_coords(_comm, _rank, DIMgeom, _coords));
	global_log->info() << "MPI coordinate of current process: " << _coords[0] << ", " << _coords[1] << ", " << _coords[2] << endl;
}

DomainDecomposition::~DomainDecomposition() {
	MPI_Comm_free(&_comm);
}

void DomainDecomposition::initCommunicationPartners(double cutoffRadius, Domain * domain, ParticleContainer* moleculeContainer) {
	for (int d = 0; d < DIMgeom; ++d) {
		_neighbourCommunicationScheme->setCoverWholeDomain(d, _gridSize[d] == 1);
	}
	_neighbourCommunicationScheme->initCommunicationPartners(cutoffRadius, domain, this, moleculeContainer);
}

void DomainDecomposition::prepareNonBlockingStage(bool /*forceRebalancing*/, ParticleContainer* moleculeContainer,
		Domain* domain, unsigned int stageNumber) {
	if (sendLeavingWithCopies()) {
		DomainDecompMPIBase::prepareNonBlockingStageImpl(moleculeContainer, domain, stageNumber,
														 LEAVING_AND_HALO_COPIES);
	} else {
		// Would first need to send leaving, then halo -> not good for overlapping!
		global_log->error() << "nonblocking P2P using separate messages for leaving and halo is currently not "
							   "supported. Please use the indirect neighbor communication scheme!"
							<< std::endl;
		Simulation::exit(235861);
	}
}

void DomainDecomposition::finishNonBlockingStage(bool /*forceRebalancing*/, ParticleContainer* moleculeContainer,
		Domain* domain, unsigned int stageNumber) {
	if (sendLeavingWithCopies()) {
		DomainDecompMPIBase::finishNonBlockingStageImpl(moleculeContainer, domain, stageNumber,
														LEAVING_AND_HALO_COPIES);
	} else {
		// Would first need to send leaving, then halo -> not good for overlapping!
		global_log->error()
			<< "nonblocking P2P using separate messages for leaving and halo is currently not supported." << std::endl;
		Simulation::exit(235861);
	}
}

bool DomainDecomposition::queryBalanceAndExchangeNonBlocking(bool /*forceRebalancing*/,
		ParticleContainer* /*moleculeContainer*/, Domain* /*domain*/, double /*etime*/) {
	return true;
}

void DomainDecomposition::balanceAndExchange(double /*lastTraversalTime*/, bool /*forceRebalancing*/, ParticleContainer* moleculeContainer,
		Domain* domain) {
	if (sendLeavingWithCopies()) {
		global_log->debug() << "DD: Sending Leaving and Halos." << std::endl;
		DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_AND_HALO_COPIES);
	} else {
		global_log->debug() << "DD: Sending Leaving." << std::endl;
		DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_ONLY);
#ifndef MARDYN_AUTOPAS
		moleculeContainer->deleteOuterParticles();
#endif
		global_log->debug() << "DD: Sending Halos." << std::endl;
		DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);
	}
}

void DomainDecomposition::readXML(XMLfileUnits& xmlconfig) {
	DomainDecompMPIBase::readXML(xmlconfig);

	if(xmlconfig.changecurrentnode("MPIGridDims")) {
		_gridSize[0] = xmlconfig.getNodeValue_int("x", 0);
		_gridSize[1] = xmlconfig.getNodeValue_int("y", 0);
		_gridSize[2] = xmlconfig.getNodeValue_int("z", 0);
		xmlconfig.changecurrentnode("..");
		initMPIGridDims();
	}
}

double DomainDecomposition::getBoundingBoxMin(int dimension, Domain* domain) {
	return _coords[dimension] * domain->getGlobalLength(dimension) / _gridSize[dimension];
}

double DomainDecomposition::getBoundingBoxMax(int dimension, Domain* domain) {
	if(_coords[dimension] + 1 == _gridSize[dimension]){
		return domain->getGlobalLength(dimension);
	}
	return (_coords[dimension] + 1) * domain->getGlobalLength(dimension) / _gridSize[dimension];
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
	//order of ranks is important in current version!!!
#if defined(ENABLE_MPI) //evil hack to not destroy the necessary order
    int myRank;
	MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
	int numProcs;
	MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
	std::vector<std::vector<std::vector<int>>> ranks = getAllRanks();
	int myCoords[3];
	MPI_Cart_coords(_comm, myRank, 3, myCoords);
	std::vector<int> neighbours(26,-1);
	if(numProcs == 1){
		for(int i = 0; i<26; i++)
			neighbours[i] = myRank;
	}
	else{
		for(int i = 0; i<26; i++){
			int x,y,z;
			switch(i)
			{
			case 0: //faces
				x=-1;y=0;z=0;break;
			case 1:
				x=1;y=0;z=0;break;
			case 2:
				x=0;y=-1;z=0;break;
			case 3:
				x=0;y=1;z=0;break;
			case 4:
				x=0;y=0;z=-1;break;
			case 5:
				x=0;y=0;z=1;break;
			case 6: //edges
				x=-1;y=-1;z=0;break;
			case 7:
				x=1;y=1;z=0;break;
			case 8:
				x=-1;y=1;z=0;break;
			case 9:
				x=1;y=-1;z=0;break;
			case 10:
				x=-1;y=0;z=-1;break;
			case 11:
				x=1;y=0;z=1; break;
			case 12:
				x=-1;y=0;z=1;break;
			case 13:
				x=1;y=0;z=-1;break;
			case 14:
				x=0;y=-1;z=-1;break;
			case 15:
				x=0;y=1;z=1;break;
			case 16:
				x=0;y=-1;z=1;break;
			case 17:
				x=0;y=1;z=-1;break;
			case 18:
				x=-1;y=-1;z=-1;break;
			case 19: //corners
				x=1;y=1;z=1;break;
			case 20:
				x=-1;y=-1;z=1;break;
			case 21:
				x=1;y=1;z=-1;break;
			case 22:
				x=-1;y=1;z=-1;break;
			case 23:
				x=1;y=-1;z=1;break;
			case 24:
				x=-1;y=1;z=1;break;
			case 25:
				x=1;y=-1;z=-1;break;
			}
			int coordsTemp[3];
			coordsTemp[0] = (myCoords[0] + x + ranks.size()) % ranks.size();
			coordsTemp[1] = (myCoords[1] + y + ranks[0].size()) % ranks[0].size();
			coordsTemp[2] = (myCoords[2] + z + ranks[0][0].size()) % ranks[0][0].size();
			int rank;
			MPI_Cart_rank(_comm, coordsTemp, &rank);
			neighbours[i] = rank;
		}
	}
	return neighbours;
#else
	return std::vector<int>(0);
#endif
}


std::vector<std::vector<std::vector<int>>> DomainDecomposition::getAllRanks(){
#ifdef ENABLE_MPI
	std::vector<std::vector<std::vector<int>>> ranks;
	int myRank;
	MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
	int numProcessors;
	MPI_Comm_size(MPI_COMM_WORLD,&numProcessors);

	ranks.resize(_gridSize[0]);
	for(int i = 0; i < _gridSize[0]; i++){
		ranks[i].resize(_gridSize[1]);
		for(int j = 0; j < _gridSize[1]; j++){
			ranks[i][j].resize(_gridSize[2]);
		}
	}
	int coords[3];
	for(int i = 0; i < numProcessors; i++){
		MPI_Cart_coords(_comm, i, 3, coords);
		ranks[coords[0]][coords[1]][coords[2]] = i;
	}
	return ranks;
#else
	return std::vector<std::vector<std::vector<int>>>(0);
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
	bool enlarged[3][2];

	for (unsigned int d = 0; d < DIMgeom; d++) {
		haloLow[d] = haloRegion.rmin[d];
		haloHigh[d] = haloRegion.rmax[d];
		//TODO: ONLY FULL SHELL!!!
		boundaryLow[d] = haloRegion.rmin[d] - haloRegion.offset[d] * haloRegion.width; //rmin[d] if offset[d]==0
		boundaryHigh[d] = haloRegion.rmax[d] - haloRegion.offset[d] * haloRegion.width; //if offset[d]!=0 : shift by cutoff in negative offset direction
		if (_coords[d] == 0 and haloRegion.offset[d] == -1) {
			shift[d] = domain->getGlobalLength(d);
		} else if (_coords[d] == _gridSize[d] - 1 and haloRegion.offset[d] == 1) {
			shift[d] = -domain->getGlobalLength(d);
		} else{
			shift[d] = 0.;
		}
		enlarged[d][0] = false;
		enlarged[d][1] = false;
	}
	// initialize using initializer list - here a vector with one element is created
	std::vector<CommunicationPartner> temp;
	temp.push_back(CommunicationPartner(rank, haloLow, haloHigh, boundaryLow, boundaryHigh, shift, haloRegion.offset, enlarged));
	return temp;
}

