#include "KDDecomposition2.h"

#include <cfloat>
#include <sstream>
#include <fstream>
#include <climits>

#include "Domain.h"
#include "KDNode.h"
#include "molecules/Molecule.h"
#include "Simulation.h"
#include "particleContainer/ParticleContainer.h"
#include "ParticleData.h"
#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"

#include <cmath>

using namespace std;
using Log::global_log;

#define DEBUG_DECOMP

KDDecomposition2::KDDecomposition2(double cutoffRadius, Domain* domain, int updateFrequency, int fullSearchThreshold)
		: _steps(0), _frequency(updateFrequency), _fullSearchThreshold(fullSearchThreshold) {

	MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &_ownRank) );
	MPI_CHECK( MPI_Comm_size(MPI_COMM_WORLD, &_numProcs) );

	int lowCorner[KDDIM] = {0};
	int highCorner[KDDIM] = {0};
	bool coversWholeDomain[KDDIM];
	_globalNumCells = 1;
	
	for (int dim = 0; dim < KDDIM; dim++) {
		_globalCellsPerDim[dim] = (int) floor(domain->getGlobalLength(dim) / cutoffRadius);
		_globalNumCells *= _globalCellsPerDim[dim];
		highCorner[dim] = _globalCellsPerDim[dim] - 1;
		_cellSize[dim] = domain->getGlobalLength(dim) / ((double) _globalCellsPerDim[dim]);
		coversWholeDomain[dim] = true;
	}
	
	_numParticlesPerCell = new unsigned int[_globalNumCells];
	for (int i = 0; i < _globalNumCells; i++)
		_numParticlesPerCell[i] = 0;

	// create initial decomposition
	// ensure that enough cells for the number of procs are avaialble
	_decompTree = new KDNode(_numProcs, lowCorner, highCorner, 0, 0, coversWholeDomain, 0);
	if (!_decompTree->isResolvable()) {
		global_log->error() << "KDDecompsition not possible. Each process needs at least 8 cells." << endl;
		global_log->error() << "The number of Cells is only sufficient for " << _decompTree->getNumMaxProcs() << " Procs!" << endl;
		barrier(); // the messages above are only promoted to std::out if we have the barrier somehow...
		global_simulation->exit(-1);
	}
	_decompTree->buildKDTree();
	_ownArea = _decompTree->findAreaForProcess(_ownRank);

	// initialize the mpi data type for particles once in the beginning
	ParticleData::setMPIType(_mpi_Particle_data);
	KDNode::initMPIDataType();

	global_log->info() << "Created KDDecomposition2 with updateFrequency=" << _frequency << ", fullSearchThreshold=" << _fullSearchThreshold << endl;

#ifdef DEBUG_DECOMP
	global_log->info() << "Initial Decomposition: " << endl;
	if (_ownRank == 0) {
		_decompTree->printTree("");
	}
#endif
}

KDDecomposition2::~KDDecomposition2() {
	delete[] _numParticlesPerCell;
//	_decompTree->serialize(string("kddecomp.dat"));
	if (_ownRank == 0) {
		_decompTree->plotNode("kddecomp.vtu");
	}
	KDNode::shutdownMPIDataType();
}

void KDDecomposition2::readXML(XMLfileUnits& xmlconfig) {
	/* TODO: Maybe add decomposition dimensions, default auto. */
	xmlconfig.getNodeValue("updateFrequency", _frequency);
	global_log->info() << "KDDecomposition2 update frequency: " << _frequency << endl;
	xmlconfig.getNodeValue("fullSearchThreshold", _fullSearchThreshold);
	global_log->info() << "KDDecomposition2 full search threshold: " << _fullSearchThreshold << endl;
}

void KDDecomposition2::exchangeMolecules(ParticleContainer* moleculeContainer, Domain* domain) {
	balanceAndExchange(false, moleculeContainer, domain);
}

void KDDecomposition2::balance() {
	
}


void KDDecomposition2::balanceAndExchange(bool balance, ParticleContainer* moleculeContainer, Domain* domain) {
	_moleculeContainer = moleculeContainer;
	KDNode* newDecompTree = NULL;
	KDNode* newOwnArea = NULL;

	if (_steps % _frequency == 0 || _steps <= 1) {
		global_log->info() << "KDDecomposition2: rebalancing..." << endl;
		getNumParticles(moleculeContainer);
		newDecompTree = new KDNode(_numProcs, &(_decompTree->_lowCorner[0]), &(_decompTree->_highCorner[0]), 0, 0, _decompTree->_coversWholeDomain, 0);

		if (decompose(newDecompTree, newOwnArea, MPI_COMM_WORLD)) {
			global_log->warning() << "Domain too small to achieve a perfect load balancing" << endl;
		}

		completeTreeInfo(newDecompTree, newOwnArea);		global_log->info() << "KDDecomposition2: rebalancing finished" << endl;

#ifdef DEBUG_DECOMP
		if (_ownRank == 0) {
			newDecompTree->printTree("");
		}
#endif
	}

	vector<int> procsToSendTo; // all processes to which this process has to send data
	vector<int> procsToRecvFrom; // all processes from which this process has to recv data
	vector<list<Molecule*> > particlePtrsToSend; // pointer to particles to be send
	vector<ParticleData*> particlesRecvBufs; // buffer used by my recv call
	vector<int> numMolsToSend; // number of particles to be send to other procs
	vector<int> numMolsToRecv; // number of particles to be recieved from other procs
	// collect particles to be send and find out number of particles to be recieved
	if (_steps % _frequency == 0 || _steps <= 1) {
		int haloCellIdxMin[3]; // Assuming a global 3D Cell index, haloCellIdxMin[3] gives the position
		// of the low local domain corner within this global 3D cell index
		int haloCellIdxMax[3]; // same as heloCellIdxMax, only high instead of low Corner
		for (int dim = 0; dim < 3; dim++) {
			haloCellIdxMin[dim] = newOwnArea->_lowCorner[dim] - 1;
			haloCellIdxMax[dim] = newOwnArea->_highCorner[dim] + 1;
		}
		vector<int> neighbHaloAreas; // The areas (unit: cells, including halo) of the neighbouring procs
		// For each proc, 6 int values are reserved (xlow, ylow, zlow, xhigh,...)
		// These values are not used in this context.
		getOwningProcs(haloCellIdxMin, haloCellIdxMax, _decompTree, _decompTree, &procsToRecvFrom, &neighbHaloAreas);
		getPartsToSend(_ownArea, newDecompTree, moleculeContainer, domain, procsToSendTo, numMolsToSend, particlePtrsToSend);
	}
	else {
		getPartsToSend(_ownArea, _decompTree, moleculeContainer, domain, procsToSendTo, numMolsToSend, particlePtrsToSend);
		procsToRecvFrom = procsToSendTo;
	}

	sendReceiveParticleData(procsToSendTo, procsToRecvFrom, numMolsToSend, numMolsToRecv, particlePtrsToSend, particlesRecvBufs);

	if (_steps % _frequency == 0 || _steps <= 1) {
		// find out new bounding boxes (of newOwnArea)
		double bBoxMin[3];
		double bBoxMax[3];
		for (int dim = 0; dim < 3; dim++) {
			bBoxMin[dim] = (newOwnArea->_lowCorner[dim]) * _cellSize[dim];
			bBoxMax[dim] = (newOwnArea->_highCorner[dim] + 1) * _cellSize[dim];
		}
		// shift the region of the moleculeContainer and delete all particles which are no
		// longer in the new region
		// TODO Rebuild problem with particles leaving the domain
		adjustOuterParticles(newOwnArea, moleculeContainer, domain);
		moleculeContainer->rebuild(bBoxMin, bBoxMax);

		_ownArea = newOwnArea;
		delete _decompTree;
		_decompTree = newDecompTree;
	}

	double lowLimit[3];
	double highLimit[3];
	for (int dim = 0; dim < 3; dim++) {
		lowLimit[dim] = moleculeContainer->getBoundingBoxMin(dim) - moleculeContainer->get_halo_L(dim);
		highLimit[dim] = moleculeContainer->getBoundingBoxMax(dim) + moleculeContainer->get_halo_L(dim);
	}

	// store recieved molecules in the molecule container
	// TODO move this to sendReceiveParticleData?
	ParticleData newMol;
	for (int neighbCount = 0; neighbCount < (int) procsToRecvFrom.size(); neighbCount++) {
		if (procsToRecvFrom[neighbCount] == _ownRank)
			continue; // don't exchange data with the own process
		int count = 0;
		for (int i = 0; i < numMolsToRecv[neighbCount]; i++) {
			newMol = particlesRecvBufs[neighbCount][i];
			count++;
			// change coordinates (especially needed if particle was moved across boundary)
			if (newMol.r[0] < lowLimit[0])
				newMol.r[0] += domain->getGlobalLength(0);
			else if (newMol.r[0] >= highLimit[0])
				newMol.r[0] -= domain->getGlobalLength(0);

			if (newMol.r[1] < lowLimit[1])
				newMol.r[1] += domain->getGlobalLength(1);
			else if (newMol.r[1] >= highLimit[1])
				newMol.r[1] -= domain->getGlobalLength(1);

			if (newMol.r[2] < lowLimit[2])
				newMol.r[2] += domain->getGlobalLength(2);
			else if (newMol.r[2] >= highLimit[2])
				newMol.r[2] -= domain->getGlobalLength(2);
			
			Component *component = _simulation.getEnsemble()->component(newMol.cid);
			Molecule m1 = Molecule(newMol.id, component, newMol.r[0], newMol.r[1], newMol.r[2], newMol.v[0], newMol.v[1], newMol.v[2], newMol.q[0], newMol.q[1], newMol.q[2], newMol.q[3], newMol.D[0], newMol.D[1], newMol.D[2]);
			moleculeContainer->addParticle(m1);
		}
	}
	// create the copies of local molecules due to periodic boundaries
	// (only for procs covering the whole domain in one dimension)
	// (If there was a balance, all procs have to be checked)
	createLocalCopies(moleculeContainer, domain);

	for (int neighbCount = 0; neighbCount < (int) procsToRecvFrom.size(); neighbCount++) {
		if (procsToRecvFrom[neighbCount] == _ownRank)
			continue; // don't exchange data with the own process
		delete[] particlesRecvBufs[neighbCount];
	}
	if (_steps % _frequency == 0 || _steps <= 1) {
		moleculeContainer->update();
	}
	particlesRecvBufs.resize(0);
	_steps++;
}

bool KDDecomposition2::procOwnsPos(double x, double y, double z, Domain* domain) {
	if (x < getBoundingBoxMin(0, domain) || x >= getBoundingBoxMax(0, domain)) {
		return false;
	}
	else if (y < getBoundingBoxMin(1, domain) || y >= getBoundingBoxMax(1, domain)) {
		return false;
	}
	else if (z < getBoundingBoxMin(2, domain) || z >= getBoundingBoxMax(2, domain)) {
		return false;
	}
	else {
		return true;
	}
}


double KDDecomposition2::getBoundingBoxMin(int dimension, Domain* domain) {
	double globalLength = domain->getGlobalLength(dimension);
	double pos = (_ownArea->_lowCorner[dimension]) * _cellSize[dimension];
	if (pos < 0) {
		return 0;
	}
	else if (pos > globalLength) {
		return globalLength;
	}
	else {
		return pos;
	}
}

double KDDecomposition2::getBoundingBoxMax(int dimension, Domain* domain) {
	double globalLength = domain->getGlobalLength(dimension);
	double pos = (_ownArea->_highCorner[dimension] + 1) * _cellSize[dimension];
	if (pos < 0) {
		return 0;
	}
	else if (pos > globalLength) {
		return globalLength;
	}
	else {
		return pos;
	}
}

void KDDecomposition2::printDecomp(string filename, Domain* domain) {
	if (_ownRank == 0) {
		ofstream povcfgstrm(filename.c_str());
		povcfgstrm << "size " << domain->getGlobalLength(0) << " " << domain->getGlobalLength(1) << " " << domain->getGlobalLength(2) << endl;
		povcfgstrm << "decompData Regions" << endl;
		povcfgstrm.close();
	}

	for (int process = 0; process < _numProcs; process++) {
		if (_ownRank == process) {
			ofstream povcfgstrm(filename.c_str(), ios::app);
			povcfgstrm << getBoundingBoxMin(0,domain) << " " << getBoundingBoxMin(1,domain) << " "
			           << getBoundingBoxMin(2,domain) << " " << getBoundingBoxMax(0,domain) << " "
			           << getBoundingBoxMax(1,domain) << " " << getBoundingBoxMax(2,domain) << endl;
			povcfgstrm.close();
		}
		barrier();
	}
}


unsigned KDDecomposition2::Ndistribution(unsigned localN, float* minrnd, float* maxrnd) {
	int num_procs;
	MPI_CHECK( MPI_Comm_size(this->_collComm.getTopology(), &num_procs) );
	unsigned* moldistribution = new unsigned[num_procs];
	MPI_CHECK( MPI_Allgather(&localN, 1, MPI_UNSIGNED, moldistribution, 1, MPI_UNSIGNED, this->_collComm.getTopology()) );
	unsigned globalN = 0;
	for (int r = 0; r < this->_ownRank; r++)
		globalN += moldistribution[r];
	unsigned localNbottom = globalN;
	globalN += moldistribution[this->_ownRank];
	unsigned localNtop = globalN;
	for (int r = this->_ownRank + 1; r < num_procs; r++)
		globalN += moldistribution[r];
	delete[] moldistribution;
	*minrnd = (float) localNbottom / globalN;
	*maxrnd = (float) localNtop / globalN;
	return globalN;
}

void KDDecomposition2::assertIntIdentity(int IX) {
	if (this->_ownRank) {
		MPI_CHECK( MPI_Send(&IX, 1, MPI_INT, 0, 2 * _ownRank + 17, this->_collComm.getTopology()) );
	}
	else {
		int recv;
		int num_procs;
		MPI_CHECK( MPI_Comm_size(this->_collComm.getTopology(), &num_procs) );
		MPI_Status s;
		for (int i = 1; i < num_procs; i++) {
			MPI_CHECK( MPI_Recv(&recv, 1, MPI_INT, i, 2 * i + 17, this->_collComm.getTopology(), &s) );
			if (recv != IX) {
				global_log->error() << "IX is " << IX << " for rank 0, but " << recv << " for rank " << i << "." << endl;
				MPI_Abort(MPI_COMM_WORLD, 911);
			}
		}
		global_log->debug() << "IX = " << recv << " for all " << num_procs << " ranks." << endl;
	}
}

void KDDecomposition2::assertDisjunctivity(TMoleculeContainer* mm) {
	Molecule* m;

	if (_ownRank) {
		int num_molecules = mm->getNumberOfParticles();
		unsigned long *tids;
		tids = new unsigned long[num_molecules];

		int i = 0;
		for (m = mm->begin(); m != mm->end(); m = mm->next()) {
			tids[i] = m->id();
			i++;
		}
		MPI_CHECK( MPI_Send(tids, num_molecules, MPI_UNSIGNED_LONG, 0, 2674 + _ownRank, this->_collComm.getTopology()) );
		delete[] tids;
		global_log->info() << "Data consistency checked: for results see rank 0." << endl;
	}
	else {
		map<unsigned long, int> check;
		int num_procs;
		MPI_CHECK( MPI_Comm_size(this->_collComm.getTopology(), &num_procs) );

		for (m = mm->begin(); m != mm->end(); m = mm->next())
			check[m->id()] = 0;

		MPI_Status status;
		for (int i = 1; i < num_procs; i++) {
			int num_recv = 0;
			unsigned long *recv;
			MPI_CHECK( MPI_Probe(i, 2674 + i, this->_collComm.getTopology(), &status) );
			MPI_CHECK( MPI_Get_count(&status, MPI_UNSIGNED_LONG, &num_recv) );
			recv = new unsigned long[num_recv];

			MPI_CHECK( MPI_Recv(recv, num_recv, MPI_UNSIGNED_LONG, i, 2674 + i, this->_collComm.getTopology(), &status) );
			for (int j = 0; j < num_recv; j++) {
				if (check.find(recv[j]) != check.end()) {
					global_log->error() << "Ranks " << check[recv[j]] << " and " << i << " both propagate ID " << recv[j] << endl;
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
				else
					check[recv[j]] = i;
			}
			delete[] recv;
		}
		global_log->info() << "Data consistency checked: No duplicate IDs detected among " << check.size() << " entries." << endl;
	}
}


//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//$ private Methoden, die von exchangeMolecule benvtigt werden $
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

// TODO use vector<vector> instead of vector<list<Molecule*> >& particlesToSend ?
void KDDecomposition2::getPartsToSend(KDNode* sourceArea, KDNode* decompTree, ParticleContainer* moleculeContainer, Domain* domain, vector<int>& procIDs, vector<int>& numMolsToSend, vector<list<Molecule*> >& particlesToSend) {
	int haloCellIdxMin[3]; // Assuming a global 3D Cell index, haloCellIdxMin[3] gives the position
	                       // of the low local domain corner within this global 3D cell index
	int haloCellIdxMax[3]; // same as heloCellIdxMax, only high instead of low Corner
	for (int dim = 0; dim < 3; dim++) {
		haloCellIdxMin[dim] = sourceArea->_lowCorner[dim] - 1;
		haloCellIdxMax[dim] = sourceArea->_highCorner[dim] + 1;
	}
	vector<int> neighbHaloAreas; // The areas (unit: cells, including halo) of the neighbouring procs
	                             // For each proc, 6 int values are reserved (xlow, ylow, zlow, xhigh,...)
	getOwningProcs(haloCellIdxMin, haloCellIdxMax, decompTree, decompTree, &procIDs, &neighbHaloAreas);
	particlesToSend.resize(procIDs.size());
	numMolsToSend.resize(procIDs.size());

	double regToSendLow[3];  // Region that belongs to a neighbouring process
	double regToSendHigh[3]; // -> regToSendLow
	double shiftRegLow[3];   // As not only the overlap with a neighbours region, but also with some
	                         // of the periodic copies of that region has to be considered, those
	                         // periodic copies have to be calculated. shiftRegLow is used to store that
	double shiftRegHigh[3];  // -> shiftRegLow

	// For all neighbours, find the particles that they need from this process
	// and store pointers to those particles in particlesToSend[proc]
	for (int neighbCount = 0; neighbCount < (int) procIDs.size(); neighbCount++) {
		// No need to send particles to the own proc
		if (procIDs[neighbCount] == _ownRank)
			continue;

		regToSendLow[0] = neighbHaloAreas[6 * neighbCount + 0] * _cellSize[0];
		regToSendLow[1] = neighbHaloAreas[6 * neighbCount + 1] * _cellSize[1];
		regToSendLow[2] = neighbHaloAreas[6 * neighbCount + 2] * _cellSize[2];
		regToSendHigh[0] = (neighbHaloAreas[6 * neighbCount + 3] + 1) * _cellSize[0];
		regToSendHigh[1] = (neighbHaloAreas[6 * neighbCount + 4] + 1) * _cellSize[1];
		regToSendHigh[2] = (neighbHaloAreas[6 * neighbCount + 5] + 1) * _cellSize[2];
		// Not only the neighbouring region itself, but also the periodic copies have to be checked
		int shift[3] = { 0, 0, 0 };
		for (int dim = 0; dim < 3; dim++) {
			// If the neighbouring region overlaps the left side of the global domain,
			// a copy of the neighb. region which is shifted to the right has to be examined
			if (regToSendLow[dim] < 0.0 && regToSendHigh[dim] <= domain->getGlobalLength(dim)) {
				shift[dim] = 1;
			}
			// same as before, but with shift to the left
			else if (regToSendLow[dim] >= 0.0 && regToSendHigh[dim] > domain->getGlobalLength(dim)) {
				shift[dim] = -1;
			}
			// The other cases:
			// neither overlap left or right --> no copies necessary
			// overlap on both sides --> The neighbouring area already covers the whole domain
			//                           (in that dimension), so shifted copies could not
			//                           cover more
		}

		particlesToSend[neighbCount].clear();
		for (int iz = 0; iz <= 1; iz++) {
			if (iz == 1 && shift[2] == 0) break; // no shift in z-direction
			for (int iy = 0; iy <= 1; iy++) {
				if (iy == 1 && shift[1] == 0) break; // no shift in y-direction
				for (int ix = 0; ix <= 1; ix++) {
					if (ix == 1 && shift[0] == 0) break; // no shift in x-direction
					shiftRegLow[0] = regToSendLow[0] + ix * shift[0] * domain->getGlobalLength(0);
					shiftRegLow[1] = regToSendLow[1] + iy * shift[1] * domain->getGlobalLength(1);
					shiftRegLow[2] = regToSendLow[2] + iz * shift[2] * domain->getGlobalLength(2);
					shiftRegHigh[0] = regToSendHigh[0] + ix * shift[0] * domain->getGlobalLength(0);
					shiftRegHigh[1] = regToSendHigh[1] + iy * shift[1] * domain->getGlobalLength(1);
					shiftRegHigh[2] = regToSendHigh[2] + iz * shift[2] * domain->getGlobalLength(2);
					moleculeContainer->getRegion(shiftRegLow, shiftRegHigh, particlesToSend[neighbCount]);
				}
			}
		}
		// store number of particles to be sent to the neighbour
		numMolsToSend[neighbCount] = particlesToSend[neighbCount].size();
	}
}


void KDDecomposition2::sendReceiveParticleData(vector<int>& procsToSendTo, vector<int>& procsToRecvFrom, vector<int>& numMolsToSend, vector<int>& numMolsToRecv, /*vector<ParticleData*>& particlesSendBufs*/ std::vector<std::list<Molecule*> >& particlePtrsToSend, vector<ParticleData*>& particlesRecvBufs) {

	particlesRecvBufs.resize(procsToRecvFrom.size());
	numMolsToRecv.resize(procsToRecvFrom.size());
	vector<ParticleData*> particlesSendBufs;

	// Initialise send and recieve buffers
	particlesSendBufs.resize(procsToSendTo.size());
	for (int neighbCount = 0; neighbCount < (int) procsToSendTo.size(); neighbCount++) {
		if (procsToSendTo[neighbCount] == _ownRank)
			continue; // don't exchange data with the own process
		particlesSendBufs[neighbCount] = new ParticleData[numMolsToSend[neighbCount]];
	}

	// Fill send buffer with particle data
	for (int neighbCount = 0; neighbCount < (int) procsToSendTo.size(); neighbCount++) {
		if (procsToSendTo[neighbCount] == _ownRank)
			continue; // don't exchange data with the own process
		list<Molecule*>::iterator particleIter;
		int partCount = 0;

		for (particleIter = particlePtrsToSend[neighbCount].begin(); particleIter != particlePtrsToSend[neighbCount].end(); particleIter++) {
			ParticleData::MoleculeToParticleData(particlesSendBufs[neighbCount][partCount], **particleIter);
			partCount++;
		}
	}

	vector<MPI_Request> request(procsToRecvFrom.size());
	vector<MPI_Request> sendRequests(procsToSendTo.size());

	// send all particles
	for (int neighbCount = 0; neighbCount < (int) procsToSendTo.size(); neighbCount++) {
		if (procsToSendTo[neighbCount] == _ownRank)
			continue; // don't exchange data with the own process
		MPI_CHECK( MPI_Isend(particlesSendBufs[neighbCount], numMolsToSend[neighbCount], _mpi_Particle_data, procsToSendTo[neighbCount], 0, MPI_COMM_WORLD, &sendRequests[neighbCount]) );
//		cout << "[" << _ownRank << "] send " << numMolsToSend[neighbCount] << " particles to rank " << procsToSendTo[neighbCount] << endl;
	}

	// get number of particles to receive
	for (int neighbCount = 0; neighbCount < (int) procsToRecvFrom.size(); neighbCount++) {
		MPI_Status status;
		if (procsToRecvFrom[neighbCount] == _ownRank)
			continue; // don't exchange data with the own process
		MPI_CHECK( MPI_Probe(procsToRecvFrom[neighbCount], 0, MPI_COMM_WORLD, &status) );
		MPI_CHECK( MPI_Get_count(&status, _mpi_Particle_data, &(numMolsToRecv[neighbCount])) );
//		cout << "[" << _ownRank << "] rcv " << numMolsToRecv[neighbCount] << " particles from rank " << procsToRecvFrom[neighbCount] << endl;
	}

	// create receive buffers
	for (int neighbCount = 0; neighbCount < (int) procsToRecvFrom.size(); neighbCount++) {
		if (procsToRecvFrom[neighbCount] == _ownRank)
			continue; // don't exchange data with the own process
		particlesRecvBufs[neighbCount] = new ParticleData[numMolsToRecv[neighbCount]];
	}

	// initiate recv calls
	for (int neighbCount = 0; neighbCount < (int) procsToRecvFrom.size(); neighbCount++) {
		if (procsToRecvFrom[neighbCount] == _ownRank)
			continue; // don't exchange data with the own process
		MPI_CHECK( MPI_Irecv(particlesRecvBufs[neighbCount], numMolsToRecv[neighbCount], _mpi_Particle_data, procsToRecvFrom[neighbCount], 0, MPI_COMM_WORLD, &request[neighbCount]) );
	}

	// wait for the completion of all send calls
	for (int neighbCount = 0; neighbCount < (int) procsToSendTo.size(); neighbCount++) {
		MPI_Status status;
		if (procsToSendTo[neighbCount] == _ownRank)
			continue; // don't exchange data with the own process
		MPI_CHECK( MPI_Wait(&sendRequests[neighbCount], &status) );
	}

	// wait for the completion of all recv calls
	for (int neighbCount = 0; neighbCount < (int) procsToRecvFrom.size(); neighbCount++) {
		MPI_Status status;
		if (procsToRecvFrom[neighbCount] == _ownRank)
			continue; // don't exchange data with the own process
		MPI_CHECK( MPI_Wait(&request[neighbCount], &status) );
	}

	// free memory of send buffers
	for (int neighbCount = 0; neighbCount < (int) procsToSendTo.size(); neighbCount++) {
		if (procsToSendTo[neighbCount] == _ownRank)
			continue; // don't exchange data with the own process
		delete[] particlesSendBufs[neighbCount];
	}
}


void KDDecomposition2::adjustOuterParticles(KDNode*& newOwnArea, ParticleContainer* moleculeContainer, Domain* domain) {

	Molecule* molPtr;

	double haloBBoxMin[3];
	double haloBBoxMax[3];
	for (int dim = 0; dim < 3; dim++) {
		haloBBoxMin[dim] = (newOwnArea->_lowCorner[dim] - 1) * _cellSize[dim];
		haloBBoxMax[dim] = (newOwnArea->_highCorner[dim] + 2) * _cellSize[dim];
	}

	molPtr = moleculeContainer->begin();
	while (molPtr != moleculeContainer->end()) {
		const double& x = molPtr->r(0);
		const double& y = molPtr->r(1);
		const double& z = molPtr->r(2);

		if (x >= haloBBoxMax[0] && x >= domain->getGlobalLength(0)) molPtr->setr(0, x - domain->getGlobalLength(0));
		else if (x < haloBBoxMin[0] && x < 0.0)                     molPtr->setr(0, x + domain->getGlobalLength(0));
		if (y >= haloBBoxMax[1] && y >= domain->getGlobalLength(1)) molPtr->setr(1, y - domain->getGlobalLength(1));
		else if (y < haloBBoxMin[1] && y < 0.0)                     molPtr->setr(1, y + domain->getGlobalLength(1));
		if (z >= haloBBoxMax[2] && z >= domain->getGlobalLength(2)) molPtr->setr(2, z - domain->getGlobalLength(2));
		else if (z < haloBBoxMin[2] && z < 0.0)                     molPtr->setr(2, z + domain->getGlobalLength(2));

		molPtr = moleculeContainer->next();
	}
}

void KDDecomposition2::createLocalCopies(ParticleContainer* moleculeContainer, Domain* domain) {
	Molecule* molPtr;
	// molecules that have to be copied, get a new position
	double newPosition[3];

	for (unsigned short d = 0; d < 3; ++d) {
		// only if the local area covers the whole domain in one dimension, copies have to
		// be created. Otherwise, the copies are on other procs and are send there
		// by the method exchangeMolecules
		//if(not(_ownArea->_coversWholeDomain[d]) and not(balance)) continue;

		molPtr = moleculeContainer->begin();
		while (molPtr != moleculeContainer->end()) {
			const double& rd = molPtr->r(d);
			int copy = 0; // -1: copy to left, 1: copy to right, 0: don't copy
			if (rd < moleculeContainer->get_halo_L(d))
				copy = 1;
			if (rd >= domain->getGlobalLength(d) - moleculeContainer->get_halo_L(d))
				copy = -1;
			if (copy != 0) {
				// determine the position for the copy of the molecule
				for (unsigned short d2 = 0; d2 < 3; d2++) {
					// when moving parallel to the coordinate d2 to another process, the
					// local coordinates in d2 change
					if (d2 == d)
						newPosition[d2] = rd + copy * domain->getGlobalLength(d2);
					else
						newPosition[d2] = molPtr->r(d2);
				}
				Component* component = _simulation.getEnsemble()->component(molPtr->componentid());
				Molecule m1 = Molecule(molPtr->id(),component,
				                       newPosition[0], newPosition[1], newPosition[2],
				                       molPtr->v(0),molPtr->v(1),molPtr->v(2),
				                       molPtr->q().qw(),molPtr->q().qx(),molPtr->q().qy(),molPtr->q().qz(),
				                       molPtr->D(0),molPtr->D(1),molPtr->D(2));
				moleculeContainer->addParticle(m1);
			}
			molPtr = moleculeContainer->next();
		}
	}
}



void KDDecomposition2::completeTreeInfo(KDNode*& root, KDNode*& ownArea) {

	int numElementsToRecv = root->_numProcs * 2 - 1;
	vector<KDNode*> ptrToAllNodes;
	vector<int> child1;
	vector<int> child2;
	ptrToAllNodes.resize(numElementsToRecv);
	child1.resize(numElementsToRecv);
	child2.resize(numElementsToRecv);

	KDNode* oldNode = root;
	// each process walks down the tree until an owned element is reached
	while (oldNode->_owningProc != _ownRank) {
		if (_ownRank < oldNode->_owningProc + oldNode->_child1->_numProcs) {
			oldNode = oldNode->_child1;
		}
		else {
			oldNode = oldNode->_child2;
		}
	}

	int nextSendingProcess = 0;
	for (int nodeID = 0; nodeID < numElementsToRecv; nodeID++) {

		KDNode::MPIKDNode mpiKDNode;
		if (oldNode->_nodeID == nodeID) {
			mpiKDNode = oldNode->getMPIKDNode();

			if (oldNode->_numProcs > 1) {
				oldNode = oldNode->_child1;
			}
		}

		MPI_CHECK( MPI_Bcast(&mpiKDNode, 1, KDNode::MPIKDNode::Datatype, nextSendingProcess, MPI_COMM_WORLD) );
		bool coversAll[3];

		coversAll[0] = mpiKDNode.getCoversWholeDomain(0);
		coversAll[1] = mpiKDNode.getCoversWholeDomain(1);
		coversAll[2] = mpiKDNode.getCoversWholeDomain(2);

		ptrToAllNodes[nodeID] = new KDNode(mpiKDNode.getNumProcs(), mpiKDNode.getLowCorner(),
				mpiKDNode.getHighCorner(), mpiKDNode.getNodeID(), mpiKDNode.getOwningProc(), coversAll, mpiKDNode.getLevel());
		ptrToAllNodes[nodeID]->_load = mpiKDNode.getLoad();
		ptrToAllNodes[nodeID]->_optimalLoadPerProcess = mpiKDNode.getOptimalLoadPerProcess();
		ptrToAllNodes[nodeID]->_expectedDeviation = mpiKDNode.getExpectedDeviation();
		ptrToAllNodes[nodeID]->_deviation = mpiKDNode.getDeviation();
		child1[nodeID] = mpiKDNode.getFirstChildID();
		child2[nodeID] = mpiKDNode.getSecondChildID();
		nextSendingProcess = mpiKDNode.getNextSendingProcess();
	}

	// connect all Nodes of the tree
	for (int nodeID = 0; nodeID < numElementsToRecv; nodeID++) {
		// if node is not leaf node, connect to the two children
		if (child1[nodeID] >= 0) {
			ptrToAllNodes[nodeID]->_child1 = ptrToAllNodes[child1[nodeID]];
			ptrToAllNodes[nodeID]->_child2 = ptrToAllNodes[child2[nodeID]];
		}
	}

	root = ptrToAllNodes[0];
	ownArea = ptrToAllNodes[ownArea->_nodeID];
}

void KDDecomposition2::printDecompTrees(KDNode* root) {
// use std::cout as I want to have all nodes at all processes printed
	for (int process = 0; process < _numProcs; process++) {
		if (_ownRank == process) {
			std::cout << "DecompTree at process " << process << endl;
			_decompTree->printTree();
		}
		barrier();
	}
}

#ifdef DEBUG_DECOMP
void printChildrenInfo(std::ofstream& filestream, KDNode* node, double minDev) {
	for (int i = 0; i < node->_level; i++) { filestream << "   ";}
	filestream << " * " << "load=" << node->_load << " optLoad=" << node->_optimalLoadPerProcess << " expDev=" << node->_expectedDeviation << " minDev=" << minDev << endl;
	for (int i = 0; i < node->_level; i++) { filestream << "   ";}
	filestream << "   [" << node->_child1->_lowCorner[0] << "," << node->_child1->_lowCorner[1] << "," << node->_child1->_lowCorner[2] << "]"
			<< "[" << node->_child1->_highCorner[0] << "," << node->_child1->_highCorner[1] << "," << node->_child1->_highCorner[2] << "]"
			<< " load=" << node->_child1->_load << " #procs=" << node->_child1->_numProcs << " avgLoad=" << node->_child1->calculateAvgLoadPerProc() << endl;
	for (int i = 0; i < node->_level; i++) { filestream << "   ";}
	filestream << "   [" << node->_child2->_lowCorner[0] << "," << node->_child2->_lowCorner[1] << "," << node->_child2->_lowCorner[2] << "]"
				<< "[" << node->_child2->_highCorner[0] << "," << node->_child2->_highCorner[1] << "," << node->_child2->_highCorner[2] << "]"
				<< " load=" << node->_child2->_load << " #procs=" << node->_child2->_numProcs << " avgLoad=" << node->_child2->calculateAvgLoadPerProc() << endl;
}
#endif

bool KDDecomposition2::decompose(KDNode* fatherNode, KDNode*& ownArea, MPI_Comm commGroup) {
	return decompose(fatherNode, ownArea, commGroup, FLT_MAX);
}

bool KDDecomposition2::decompose(KDNode* fatherNode, KDNode*& ownArea, MPI_Comm commGroup, const double globalMinimalDeviation) {
	bool domainTooSmall = false;
	// recursion termination criterion
	if (fatherNode->_numProcs == 1) {
		// own area must belong to this process!
		assert(fatherNode->_owningProc == _ownRank);
		ownArea = fatherNode;
		fatherNode->calculateDeviation();
		return domainTooSmall;
	}

	std::list<KDNode*> subdivisions;
	domainTooSmall = calculateAllSubdivisions(fatherNode, subdivisions, commGroup);
	assert(subdivisions.size() > 0);

	KDNode* bestSubdivision = NULL;
	double minimalDeviation = globalMinimalDeviation;
	list<KDNode*>::iterator iter = subdivisions.begin();
	int iterations = 0;
	int log2 = 0;
	while ((_numProcs >> log2) > 1) {
		log2++;
	}
	// if we are near the root of the tree, we just take the first best subdivision
	int maxIterations = 1;
	if (fatherNode->_level > (log2 - _fullSearchThreshold)) {
		maxIterations = INT_MAX;
	}

#ifdef DEBUG_DECOMP
	std::stringstream fname;
	fname << "Div_proc_" << _ownRank << "_step_" << _steps << ".txt";
	std::ofstream filestream(fname.str().c_str(), ios::app);
	filestream.precision(8);
	for (int i = 0; i < fatherNode->_level; i++) { filestream << "   ";}
	filestream << "Division at rank=" << _ownRank << " for [" << fatherNode->_lowCorner[0]
               << ","<< fatherNode->_lowCorner[1] << "," << fatherNode->_lowCorner[2] << "] [" << fatherNode->_highCorner[0]
               << ","<< fatherNode->_highCorner[1] << "," << fatherNode->_highCorner[2] << "] " <<
               "level=" << fatherNode->_level << " #divisions=" << subdivisions.size() << endl;
#endif

	while (iter !=  subdivisions.end() && (iterations < maxIterations) && (*iter)->_expectedDeviation < minimalDeviation) {
		iterations++;
#ifdef DEBUG_DECOMP
		printChildrenInfo(filestream, *iter, minimalDeviation);
#endif
		vector<int> origRanks;
		int newNumProcs;
		if (_ownRank < (*iter)->_child2->_owningProc) {
			origRanks.resize((*iter)->_child1->_numProcs);
			for (int i = 0; i < (*iter)->_child1->_numProcs; i++) {
				origRanks[i] = i;
			}
			newNumProcs = (*iter)->_child1->_numProcs;
		}
		else {
			origRanks.resize((*iter)->_child2->_numProcs);
			for (int i = 0; i < (*iter)->_child2->_numProcs; i++) {
				origRanks[i] = i + (*iter)->_child1->_numProcs;
			}
			newNumProcs = (*iter)->_child2->_numProcs;
		}

		MPI_Comm newComm;
		MPI_Group origGroup, newGroup;

		MPI_CHECK( MPI_Comm_group(commGroup, &origGroup) );
		MPI_CHECK( MPI_Group_incl(origGroup, newNumProcs, &origRanks[0], &newGroup) );
		MPI_CHECK( MPI_Comm_create(commGroup, newGroup, &newComm) );

		KDNode* newOwnArea = NULL;
		double deviationChildren[] = {0.0, 0.0};

		if (_ownRank < (*iter)->_child2->_owningProc) {
			// do not use the function call directly in the logical expression, as it may
			// not be executed due to conditional / short-circuit evaluation!
			bool subdomainTooSmall = decompose((*iter)->_child1, newOwnArea, newComm, minimalDeviation);
			deviationChildren[0] = (*iter)->_child1->_deviation;
			domainTooSmall = (domainTooSmall || subdomainTooSmall);
		} else {
			assert(_ownRank >= (*iter)->_child2->_owningProc);
			bool subdomainTooSmall = decompose((*iter)->_child2, newOwnArea, newComm, minimalDeviation);
			deviationChildren[1] = (*iter)->_child2->_deviation;
			domainTooSmall = (domainTooSmall || subdomainTooSmall);
		}

		MPI_CHECK( MPI_Group_free(&newGroup));
		MPI_CHECK( MPI_Comm_free(&newComm) );
		MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, deviationChildren, 2, MPI_DOUBLE, MPI_SUM, commGroup));
		(*iter)->_child1->_deviation = deviationChildren[0];
		(*iter)->_child2->_deviation = deviationChildren[1];
		(*iter)->calculateDeviation();

#ifdef DEBUG_DECOMP
		for (int i = 0; i < fatherNode->_level; i++) { filestream << "   ";}
		filestream << "   deviation=" << (*iter)->_deviation << " (ch1:" << deviationChildren[0] << "ch2:" << deviationChildren[1] << endl;
#endif
		if ((*iter)->_deviation < minimalDeviation) {
			delete bestSubdivision;
			bestSubdivision = *iter;
			minimalDeviation = (*iter)->_deviation;
			ownArea = newOwnArea;
		} else {
			delete *iter;
		}
		iter++;
	}

	while (iter !=  subdivisions.end()) {
		delete *iter;
		iter++;
	}

	// reassign children and delete cloned node, if a better solution
	// was found in this subtree.
	if (bestSubdivision == NULL) {
		fatherNode->_deviation = FLT_MAX;
	} else {
		*fatherNode = *bestSubdivision;
		bestSubdivision->_child1 = NULL;
		bestSubdivision->_child2 = NULL;
		delete bestSubdivision;
	}
#ifdef DEBUG_DECOMP
	filestream.close();
#endif
	return domainTooSmall;
}


bool KDDecomposition2::calculateAllSubdivisions(KDNode* node, std::list<KDNode*>& subdivededNodes, MPI_Comm commGroup) {
	bool domainTooSmall = false;
	vector<vector<double> > costsLeft(3);
	vector<vector<double> > costsRight(3);
	calculateCostsPar(node, costsLeft, costsRight, commGroup);

	for (int dim = 0; dim < 3; dim++) {

		// if a node has some more processors, we probably don't have to find the
		// "best" partitioning, but it's sufficient to divide the node in the middle
		// and shift the processor count according to the load imbalance.

		// loop only from 1 to max-1 (instead 0 to max) to avoid 1-cell-regions
		int startIndex = 1;
		int maxEndIndex = node->_highCorner[dim] - node->_lowCorner[dim] - 1;
		int endIndex = maxEndIndex;
		if (node->_numProcs > _fullSearchThreshold) {
			startIndex = max(startIndex, (node->_highCorner[dim] - node->_lowCorner[dim] - 1) / 2);
			endIndex = min(endIndex, startIndex + 1);
		}

		int i = startIndex;
		while ( (i < endIndex) || (subdivededNodes.size() == 0 && i < maxEndIndex) ) {
		//for (int i = startIndex; i < endIndex; i++) {
			double optCostPerProc = (costsLeft[dim][i] + costsRight[dim][i]) / ((double) node->_numProcs);
			int optNumProcsLeft = min(round(costsLeft[dim][i] / optCostPerProc), (double) (node->_numProcs - 1));
			int numProcsLeft = (int) max(1, optNumProcsLeft);

			KDNode* clone = new KDNode(*node);
			if (clone->_level == 0) {
				clone->_optimalLoadPerProcess = optCostPerProc;
			}

			clone->split(dim, node->_lowCorner[dim] + i, numProcsLeft);
			if ( (unsigned int) (clone->_child1->_numProcs + clone->_child2->_numProcs) >
			        (clone->_child1->getNumMaxProcs() + clone->_child2->getNumMaxProcs())) {
				domainTooSmall = true;
				delete clone;
				// domain is not resolvable at all -> continue
				// I think, this case should not happen, but Martin Buchholz coded it in his code...
				// OK, it happens...!
				i++;
				continue;
			}

			while ( (! clone->_child1->isResolvable()) && clone->_child2->isResolvable()) {
				// shift procs to child 2, adapt owner of child 2
				clone->_child1->_numProcs--;
				clone->_child2->_numProcs++;
				clone->_child2->_owningProc--;
				clone->_child2->_nodeID = clone->_nodeID + 2 * clone->_child1->_numProcs;
				domainTooSmall = true;
			}

			while ( clone->_child1->isResolvable() && (! clone->_child2->isResolvable())) {
				// shift procs to child 1, , adapt owner of child 2
				clone->_child1->_numProcs++;
				clone->_child2->_numProcs--;
				clone->_child2->_owningProc++;
				clone->_child2->_nodeID = clone->_nodeID + 2 * clone->_child1->_numProcs;
				domainTooSmall = true;
			}

			// In this place, MBu had some processor shifting in his algorithm.
			// I believe it to be unneccessary, as the ratio of left and right processors
			// is chosen according to the load ratio (I use round instead of floor).
			if ((clone->_child1->_numProcs <= 0 || clone->_child1->_numProcs >= node->_numProcs) ||
					(clone->_child2->_numProcs <= 0 || clone->_child2->_numProcs >= node->_numProcs) ){
				global_log->error() << "ERROR in calculateAllSubdivisions(), part of the domain was not assigned to a proc" << endl;
				global_simulation->exit(1);
			}
			assert( clone->_child1->isResolvable() && clone->_child2->isResolvable() );

			clone->_child1->_load = costsLeft[dim][i];
			clone->_child2->_load = costsRight[dim][i];
			clone->_load = costsLeft[dim][i] + costsRight[dim][i];
			clone->calculateExpectedDeviation();

			// sort node according to expected deviation
			list<KDNode*>::iterator iter = subdivededNodes.begin();
			while (iter != subdivededNodes.end() && ((*iter)->_expectedDeviation < clone->_expectedDeviation)) {
				iter++;
			}
			subdivededNodes.insert(iter, clone);
			i++;
		}
	}
	return domainTooSmall;
}

/**
 * Calculate the division cost for all possible divisions of the node area.
 *
 * Is done in the following way:
 * - get the number of particles per cell globally (i.e. for all cell in the domain)
 * - get the number of particle pairs per cell globally
 * - then calculate the costs for all possible subdivisions.
 */
void KDDecomposition2::calculateCostsPar(KDNode* area, vector<vector<double> >& costsLeft, vector<vector<double> >& costsRight, MPI_Comm commGroup) {

	vector<vector<double> > cellCosts;
	cellCosts.resize(3);

	int newRank;
	MPI_Group group;

	MPI_CHECK( MPI_Comm_group(commGroup, &group) );
	MPI_CHECK( MPI_Group_rank(group, &newRank) );

	int dimStartIndex[3];
	int dimStopIndex[3];
	dimStartIndex[0] = 0;
	dimStopIndex[0] = (area->_highCorner[0] - area->_lowCorner[0]);
	for (int dim = 1; dim < 3; dim++) {
		dimStartIndex[dim] = dimStopIndex[dim - 1] + 1;
		dimStopIndex[dim] = dimStopIndex[dim - 1] + (area->_highCorner[dim] - area->_lowCorner[dim] + 1);
	}
	int sumNumLayers = dimStopIndex[2] + 1;
	int startIndex = sumNumLayers * newRank / area->_numProcs;
	int stopIndex = sumNumLayers * (newRank + 1) / area->_numProcs - 1;
	for (int dim = 0; dim < 3; dim++) {

		cellCosts[dim].resize(area->_highCorner[dim] - area->_lowCorner[dim] + 1, 0.0);
		costsLeft[dim].resize(area->_highCorner[dim] - area->_lowCorner[dim] + 1, 0.0);
		costsRight[dim].resize(area->_highCorner[dim] - area->_lowCorner[dim] + 1, 0.0);

		// if this process doesn't has to do anything in this dimension, continue
		if (startIndex > dimStopIndex[dim] || stopIndex < dimStartIndex[dim]) continue;

		int loopstart = max(0, startIndex - dimStartIndex[dim]);
		int loopend = min(area->_highCorner[dim] - area->_lowCorner[dim], (area->_highCorner[dim] - area->_lowCorner[dim]) + stopIndex - dimStopIndex[dim]);

		bool sendCostValue = false;
		bool recvCostValue = false;
		if (loopstart > 0)
			recvCostValue = true;
		if (loopend < area->_highCorner[dim] - area->_lowCorner[dim])
			sendCostValue = true;

		for (int i_dim = loopstart; i_dim <= loopend; i_dim++) {

			if (i_dim == 0)
				cellCosts[dim][i_dim] = 0;
			else
				cellCosts[dim][i_dim] = cellCosts[dim][i_dim - 1];

			int dim1, dim2;
			if (dim == 0) {
				dim1 = 1;
				dim2 = 2;
			}
			else if (dim == 1) {
				dim1 = 0;
				dim2 = 2;
			}
			else {
				dim1 = 0;
				dim2 = 1;
			}
			for (int i_dim1 = 0; i_dim1 <= area->_highCorner[dim1] - area->_lowCorner[dim1]; i_dim1++) {
				for (int i_dim2 = 0; i_dim2 <= area->_highCorner[dim2] - area->_lowCorner[dim2]; i_dim2++) {
					int numParts = (int) _numParticlesPerCell[getGlobalIndex(dim, dim1, dim2, i_dim, i_dim1, i_dim2, area)];

					// #######################
					// ## Cell Costs        ##
					// #######################
					cellCosts[dim][i_dim] += (double) (numParts * numParts);

					// all Neighbours
					for (int neighInd_divDim = i_dim - 1; neighInd_divDim <= i_dim + 1; neighInd_divDim++) {
						int nI_dim = neighInd_divDim;
						if (nI_dim + area->_lowCorner[dim] < 0)
							nI_dim = _globalCellsPerDim[dim] - 1;
						if (nI_dim + area->_lowCorner[dim] >= _globalCellsPerDim[dim])
							nI_dim = 0;

						for (int neighInd_dim1 = i_dim1 - 1; neighInd_dim1 <= i_dim1 + 1; neighInd_dim1++) {
							// adjust index in case of periodic boundary
							int nI_dim1 = neighInd_dim1;
							if (nI_dim1 + area->_lowCorner[dim1] < 0)
								nI_dim1 = _globalCellsPerDim[dim1] - 1;
							if (nI_dim1 + area->_lowCorner[dim1] >= _globalCellsPerDim[dim1])
								nI_dim1 = 0;

							for (int neighInd_dim2 = i_dim2 - 1; neighInd_dim2 <= i_dim2 + 1; neighInd_dim2++) {
								// adjust index in case of periodic boundary
								int nI_dim2 = neighInd_dim2;
								if (nI_dim2 + area->_lowCorner[dim2] < 0)
									nI_dim2 = _globalCellsPerDim[dim2] - 1;
								if (nI_dim2 + area->_lowCorner[dim2] >= _globalCellsPerDim[dim2])
									nI_dim2 = 0;
								// count only forward neighbours
								int numPartsNeigh = (int) _numParticlesPerCell[getGlobalIndex(dim, dim1, dim2, nI_dim, nI_dim1, nI_dim2, area)];
								cellCosts[dim][i_dim] += 0.5 * (double) (numParts * numPartsNeigh);
							}
						}
					}
				}
			}
		}

		// exchange intermediate calc costs
		MPI_Status recvStat;
		double tempRecvCosts, tempSendCosts;
		tempSendCosts = 0;
		if (recvCostValue) {
			MPI_CHECK( MPI_Recv(&tempRecvCosts, 1, MPI_DOUBLE, _ownRank - 1, 123, MPI_COMM_WORLD, &recvStat) );
			if (sendCostValue) {
				tempSendCosts = tempRecvCosts;
			}
		}
		if (sendCostValue) {
			tempSendCosts += cellCosts[dim][loopend];
			MPI_CHECK( MPI_Send(&tempSendCosts, 1, MPI_DOUBLE, _ownRank + 1, 123, MPI_COMM_WORLD) );
		}
		if (recvCostValue) {
			for (int i_dim = loopstart; i_dim <= loopend; i_dim++) {
				cellCosts[dim][i_dim] += tempRecvCosts;
			}
		}
	}
	for (int dim = 0; dim < 3; dim++) {
		vector<double> cellCostsSum;
		cellCostsSum.resize(area->_highCorner[dim] - area->_lowCorner[dim] + 1, 0.0);

		int size2 = cellCostsSum.size();
		MPI_CHECK( MPI_Allreduce(&cellCosts[dim][0], &cellCostsSum[0], size2, MPI_DOUBLE, MPI_SUM, commGroup) );

		for (int i_dim = 0; i_dim <= area->_highCorner[dim] - area->_lowCorner[dim]; i_dim++) {
			costsLeft[dim][i_dim] = cellCostsSum[i_dim];
			costsRight[dim][i_dim] = cellCostsSum[area->_highCorner[dim] - area->_lowCorner[dim]] - cellCostsSum[i_dim];
		}
	}
}


unsigned int KDDecomposition2::getGlobalIndex(int divDim, int dim1, int dim2, int index_divDim, int index_dim1, int index_dim2, KDNode* localArea) {
	int xIndex, yIndex, zIndex;
	if (divDim == 0) {
		xIndex = index_divDim + localArea->_lowCorner[divDim];
		yIndex = index_dim1 + localArea->_lowCorner[dim1];
		zIndex = index_dim2 + localArea->_lowCorner[dim2];
	}
	else if (divDim == 1) {
		xIndex = index_dim1 + localArea->_lowCorner[dim1];
		yIndex = index_divDim + localArea->_lowCorner[divDim];
		zIndex = index_dim2 + localArea->_lowCorner[dim2];
	}
	else {
		xIndex = index_dim1 + localArea->_lowCorner[dim1];
		yIndex = index_dim2 + localArea->_lowCorner[dim2];
		zIndex = index_divDim + localArea->_lowCorner[divDim];
	}
	return _globalCellsPerDim[0] * (zIndex * _globalCellsPerDim[1] + yIndex) + xIndex;
}


//##########################################################################
//##########################################################################
//###                                                                    ###
//###                  evtl unnoetige private Methoden                   ###
//###                                                                    ###
//##########################################################################
//##########################################################################

int KDDecomposition2::mod(int number, int modulo) {
	int result = number % modulo;
	if (result < 0)
		result += modulo;
	return result;
}

// TODO: this method could or should be moved to KDNode.
void KDDecomposition2::getOwningProcs(int low[KDDIM], int high[KDDIM], KDNode* decompTree, KDNode* testNode, vector<int>* procIDs, vector<int>* neighbHaloAreas) {
	// For areas overlapping the domain given by decompTree, the overlapping part is
	// mapped to the corresponding area on the other side of the domain (periodic boundary)
	// The boolean variable overlap stores for each coordinate direction whether the area overlaps.
	bool overlap[KDDIM];
	bool coversWholeDomain[KDDIM];
	for (int dim = 0; dim < KDDIM; dim++) {
		// if the domain is e.g. 10 cells long, and one proc has 0-8, halo will go from -1 to 9,
		// which equals (periodic boundary) 9 to 9, which would be the same as a proc which has onle
		// this one cell. Therefore, the flag coversWholeRegion is introduced
		// @todo Do this nicer
		if (high[dim] - low[dim] >= decompTree->_highCorner[dim] - decompTree->_lowCorner[dim])
			coversWholeDomain[dim] = true;
		else
			coversWholeDomain[dim] = false;

		low[dim] = mod(low[dim], (decompTree->_highCorner[dim] - decompTree->_lowCorner[dim] + 1));
		high[dim] = mod(high[dim], (decompTree->_highCorner[dim] - decompTree->_lowCorner[dim] + 1));

		if (low[dim] > high[dim])
			overlap[dim] = true;
		else
			overlap[dim] = false;

		if (coversWholeDomain[dim])
			low[dim] = high[dim] + 1;
	}

	// First find out whether the area intersects the area given by testNode
	bool areasIntersect = true;
	for (int dim = 0; dim < KDDIM; dim++) {
		if (not (coversWholeDomain[dim])) { // TODO: check wheter this can be removed || _ownArea->_coversWholeDomain[dim])){
			if (overlap[dim]) {
				if (high[dim] < testNode->_lowCorner[dim] && low[dim] > testNode->_highCorner[dim])
					areasIntersect = false;
			}
			else {
				if (high[dim] < testNode->_lowCorner[dim] || low[dim] > testNode->_highCorner[dim])
					areasIntersect = false;
			}
		}
	}
	if (areasIntersect) {
		// while testNode is still a inner node (more than one proc), recursively call
		// this method for the children
		if (testNode->_numProcs > 1) {
			getOwningProcs(low, high, decompTree, testNode->_child1, procIDs, neighbHaloAreas);
			getOwningProcs(low, high, decompTree, testNode->_child2, procIDs, neighbHaloAreas);
		}
		else { // leaf node found, add it (and it's area)
			procIDs->push_back(testNode->_owningProc);
			neighbHaloAreas->push_back(testNode->_lowCorner[0] - 1);
			neighbHaloAreas->push_back(testNode->_lowCorner[1] - 1);
			neighbHaloAreas->push_back(testNode->_lowCorner[2] - 1);
			neighbHaloAreas->push_back(testNode->_highCorner[0] + 1);
			neighbHaloAreas->push_back(testNode->_highCorner[1] + 1);
			neighbHaloAreas->push_back(testNode->_highCorner[2] + 1);
		}
	}
}

void KDDecomposition2::getNumParticles(ParticleContainer* moleculeContainer) {
	for (int i = 0; i < _globalNumCells; i++)
		_numParticlesPerCell[i] = 0;

	int count = 0;
	double bBMin[3]; // haloBoundingBoxMin
    /* TODO: We do not use values form bBMax anywhere ... */
	//double bBMax[3]; // haloBoundingBoxMax
	for (int dim = 0; dim < 3; dim++) {
		bBMin[dim] = moleculeContainer->getBoundingBoxMin(dim);// - moleculeContainer->get_halo_L(dim);
		//bBMax[dim] = moleculeContainer->getBoundingBoxMax(dim);// + moleculeContainer->get_halo_L(dim);
	}
	Molecule* molPtr = moleculeContainer->begin();
	while (molPtr != moleculeContainer->end()) {
		int cellIndex[3]; // 3D Cell index (local)
		int globalCellIdx[3]; // 3D Cell index (global)

		for (int dim = 0; dim < 3; dim++) {
			cellIndex[dim] = (int) floor((molPtr->r(dim) - bBMin[dim]) / _cellSize[dim]);
			globalCellIdx[dim] = _ownArea->_lowCorner[dim] + cellIndex[dim];
			if (globalCellIdx[dim] < 0)
				globalCellIdx[dim] += _globalCellsPerDim[dim];
			if (globalCellIdx[dim] >= _globalCellsPerDim[dim])
				globalCellIdx[dim] -= _globalCellsPerDim[dim];
		}

		_numParticlesPerCell[_globalCellsPerDim[0] * (globalCellIdx[2] * _globalCellsPerDim[1] + globalCellIdx[1]) + globalCellIdx[0]]++;
		molPtr = moleculeContainer->next();
		count++;
	}
	MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, _numParticlesPerCell, _globalNumCells, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD) );

}

