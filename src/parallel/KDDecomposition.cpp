#include "KDDecomposition.h" //$ Pfad aendern (parallel/...)

#include "../molecules/Molecule.h"
#include "../Domain.h"
#include "../particleContainer/ParticleContainer.h"
#include "../particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "../particleContainer/adapter/ParticlePairs2LoadCalcAdapter.h"
#include "ParticleData.h"
#include "KDNode.h"

#include <float.h>
#include <sstream>
#include <fstream>
#include <climits>

using namespace std;

KDDecomposition::KDDecomposition(double cutoffRadius, Domain* domain, double alpha, double beta)
		: _steps(0), _frequency(10), _alpha(alpha), _beta(beta) {

	MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &_ownRank) );
	MPI_CHECK( MPI_Comm_size(MPI_COMM_WORLD, &_numProcs) );

	int lowCorner[KDDIM];
	int highCorner[KDDIM];
	bool coversWholeDomain[KDDIM];

	for (int dim = 0; dim < KDDIM; dim++) {
		_globalCellsPerDim[dim] = (int) floor(domain->getGlobalLength(dim) / cutoffRadius);
		lowCorner[dim] = 0;
		highCorner[dim] = _globalCellsPerDim[dim] - 1;
		_cellSize[dim] = domain->getGlobalLength(dim) / ((double) _globalCellsPerDim[dim]);
		coversWholeDomain[dim] = true;
	}
	_globalNumCells = (_globalCellsPerDim[0]) * (_globalCellsPerDim[1]) * (_globalCellsPerDim[2]);
	_numParticlesPerCell = new unsigned char[_globalNumCells];
	for (int i = 0; i < _globalNumCells; i++)
		_numParticlesPerCell[i] = 0;

	// create initial decomposition
	// ensure that enough cells for the number of procs are avaialble
	int maxProcs = 1;
	for (int dim = 0; dim < 3; dim++) {
		maxProcs *= (highCorner[dim] - lowCorner[dim] + 1) / 2;
	}
	if (_ownRank == 0 && maxProcs < _numProcs) {
		cerr << "KDDecompsition not possible. Each Proc needs at least 8 cells." << endl;
		cerr << "the number of Cells is only sufficient for " << maxProcs << " Procs!" << endl;
		exit(1);
	}
	_decompTree = new KDNode(_numProcs, lowCorner, highCorner, 0, 0, coversWholeDomain);
	if (_numProcs > 1)
		recInitialDecomp(_decompTree, _ownArea);
	else
		_ownArea = _decompTree;
}

KDDecomposition::~KDDecomposition() {
	delete[] _numParticlesPerCell;
	MPI_Finalize();
}

void KDDecomposition::exchangeMolecules(ParticleContainer* moleculeContainer, const vector<Component>& components, Domain* domain) {
	balanceAndExchange(false, moleculeContainer, components, domain);
}

void KDDecomposition::balanceAndExchange(bool balance, ParticleContainer* moleculeContainer, const vector<Component>& components, Domain* domain) {
	_moleculeContainer = moleculeContainer;
	KDNode* newDecompTree = NULL;
	KDNode* newOwnArea = NULL;

	if (_steps % _frequency == 0 || _steps <= 1) {
		getNumParticles(moleculeContainer);
		newDecompTree = new KDNode(_numProcs, _decompTree->_lowCorner, _decompTree->_highCorner, 0, 0, _decompTree->_coversWholeDomain);
		ParticlePairs2LoadCalcAdapter* loadHandler;
		loadHandler = new ParticlePairs2LoadCalcAdapter(_globalCellsPerDim, _ownArea->_lowCorner, _cellSize, _moleculeContainer);
		ParticlePairsHandler* tempHandler = _moleculeContainer->getPairHandler();
		_moleculeContainer->setPairHandler(loadHandler);
		_moleculeContainer->traversePairs();
		_moleculeContainer->setPairHandler(tempHandler);
		_globalLoadPerCell = loadHandler->getLoad();
		if (recDecompPar(newDecompTree, newOwnArea, MPI_COMM_WORLD)) {
			cerr << "Domain too small to achieve a perfect load balancing" << endl;
		}
		completeTreeInfo(newDecompTree, newOwnArea);
		delete loadHandler;
	}

	vector<int> procsToSendTo; // all processes to which this process has to send data
	vector<int> procsToRecvFrom; // all processes from which this process has to recv data
	vector<list<Molecule*> > particlePtrsToSend; // pointer to particles to be send
	vector<ParticleData*> particlesSendBufs; // buffer used by my send call
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
	numMolsToRecv.resize(procsToRecvFrom.size());
	exchangeNumToSend(procsToSendTo, procsToRecvFrom, numMolsToSend, numMolsToRecv);

	// Initialise send and recieve buffers
	particlesSendBufs.resize(procsToSendTo.size());
	particlesRecvBufs.resize(procsToRecvFrom.size());
	for (int neighbCount = 0; neighbCount < (int) procsToSendTo.size(); neighbCount++) {
		if (procsToSendTo[neighbCount] == _ownRank)
			continue; // don't exchange data with the own process
		// TODO: replace particlePtrsToSend[neighbCount].size() by numMolsToSend[neighbCount]
		//particlesSendBufs[neighbCount] = new ParticleData[particlePtrsToSend[neighbCount].size()];
		particlesSendBufs[neighbCount] = new ParticleData[numMolsToSend[neighbCount]];
	}
	for (int neighbCount = 0; neighbCount < (int) procsToRecvFrom.size(); neighbCount++) {
		if (procsToRecvFrom[neighbCount] == _ownRank)
			continue; // don't exchange data with the own process
		particlesRecvBufs[neighbCount] = new ParticleData[numMolsToRecv[neighbCount]];
	}
	// Fill send buffer with particle data
	for (int neighbCount = 0; neighbCount < (int) procsToSendTo.size(); neighbCount++) {
		if (procsToSendTo[neighbCount] == _ownRank)
			continue; // don't exchange data with the own process
		list<Molecule*>::iterator particleIter;
		int partCount = 0;

		for (particleIter = particlePtrsToSend[neighbCount].begin(); particleIter != particlePtrsToSend[neighbCount].end(); particleIter++) {
			ParticleData::MoleculeToParticleData(particlesSendBufs[neighbCount][partCount], **particleIter);
			//  particlesSendBufs[neighbCount][partCount].id = (*particleIter)->id();
			//particlesSendBufs[neighbCount][partCount].cid = (*particleIter)->componentid();
			//particlesSendBufs[neighbCount][partCount].rx = (*particleIter)->r(0);
			//particlesSendBufs[neighbCount][partCount].ry = (*particleIter)->r(1);
			//particlesSendBufs[neighbCount][partCount].rz = (*particleIter)->r(2);
			//particlesSendBufs[neighbCount][partCount].vx = (*particleIter)->v(0);
			//particlesSendBufs[neighbCount][partCount].vy = (*particleIter)->v(1);
			//particlesSendBufs[neighbCount][partCount].vz = (*particleIter)->v(2);
			//particlesSendBufs[neighbCount][partCount].qw = (*particleIter)->q().qw();
			//particlesSendBufs[neighbCount][partCount].qx = (*particleIter)->q().qx();
			//particlesSendBufs[neighbCount][partCount].qy = (*particleIter)->q().qy();
			//particlesSendBufs[neighbCount][partCount].qz = (*particleIter)->q().qz();
			//particlesSendBufs[neighbCount][partCount].Dx = (*particleIter)->D(0);
			//particlesSendBufs[neighbCount][partCount].Dy = (*particleIter)->D(1);
			//particlesSendBufs[neighbCount][partCount].Dz = (*particleIter)->D(2);
			partCount++;

		}
	}
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

	// Transfer data
	transferMolData(procsToSendTo, procsToRecvFrom, numMolsToSend, numMolsToRecv, particlesSendBufs, particlesRecvBufs);
	double lowLimit[3];
	double highLimit[3];
	for (int dim = 0; dim < 3; dim++) {
		lowLimit[dim] = moleculeContainer->getBoundingBoxMin(dim) - moleculeContainer->get_halo_L(dim);
		highLimit[dim] = moleculeContainer->getBoundingBoxMax(dim) + moleculeContainer->get_halo_L(dim);
	}
	// store recieved molecules in the molecule container
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

			Molecule m1 = Molecule(newMol.id, newMol.cid, newMol.r[0], newMol.r[1], newMol.r[2], newMol.v[0], newMol.v[1], newMol.v[2], newMol.q[0], newMol.q[1], newMol.q[2], newMol.q[3], newMol.D[0], newMol.D[1], newMol.D[2], &components);
			//if(newMol.rx < lowLimit[0]) newMol.rx += domain->getGlobalLength(0);
			//else if(newMol.rx >= highLimit[0]) newMol.rx -= domain->getGlobalLength(0);
			//if(newMol.ry < lowLimit[1]) newMol.ry += domain->getGlobalLength(1);
			//else if(newMol.ry >= highLimit[1]) newMol.ry -= domain->getGlobalLength(1);
			//if(newMol.rz < lowLimit[2]) newMol.rz += domain->getGlobalLength(2);
			//else if(newMol.rz >= highLimit[2]) newMol.rz -= domain->getGlobalLength(2);
			//Molecule m1 = Molecule(newMol.id, newMol.cid, newMol.rx, newMol.ry, newMol.rz, newMol.vx, newMol.vy, newMol.vz, newMol.qw, newMol.qx, newMol.qy, newMol.qz, newMol.Dx, newMol.Dy, newMol.Dz, &components);

			moleculeContainer->addParticle(m1);
		}
	}
	// create the copies of local molecules due to periodic boundaries
	// (only for procs covering the whole domain in one dimension)
	// (If there was a balance, all procs have to be checked)
	createLocalCopies(moleculeContainer, domain, components);

	// free memory of send and recv buffers
	for (int neighbCount = 0; neighbCount < (int) procsToSendTo.size(); neighbCount++) {
		if (procsToSendTo[neighbCount] == _ownRank)
			continue; // don't exchange data with the own process
		delete[] particlesSendBufs[neighbCount];
	}
	for (int neighbCount = 0; neighbCount < (int) procsToRecvFrom.size(); neighbCount++) {
		if (procsToRecvFrom[neighbCount] == _ownRank)
			continue; // don't exchange data with the own process
		delete[] particlesRecvBufs[neighbCount];
	}
	if (_steps % _frequency == 0 || _steps <= 1) {
		moleculeContainer->update();
	}
	particlesSendBufs.resize(0);
	particlesRecvBufs.resize(0);
	_steps++;
}

bool KDDecomposition::procOwnsPos(double x, double y, double z, Domain* domain) {
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

double KDDecomposition::guaranteedDistance(double x, double y, double z, Domain* domain) {
	double xdist = 0;
	double ydist = 0;
	double zdist = 0;
	if (x < getBoundingBoxMin(0, domain)) {
		xdist = getBoundingBoxMin(0, domain) - x;
	}
	else if (x >= getBoundingBoxMax(0, domain)) {
		xdist = x - getBoundingBoxMax(0, domain);
	}

	if (y < getBoundingBoxMin(1, domain)) {
		ydist = getBoundingBoxMin(1, domain) - y;
	}
	else if (y >= getBoundingBoxMax(1, domain)) {
		ydist = y - getBoundingBoxMax(1, domain);
	}

	if (z < getBoundingBoxMin(2, domain)) {
		zdist = getBoundingBoxMin(2, domain) - z;
	}
	else if (z >= getBoundingBoxMax(2, domain)) {
		zdist = z - getBoundingBoxMax(2, domain);
	}
	return sqrt(xdist * xdist + ydist * ydist + zdist * zdist);
}

unsigned long KDDecomposition::countMolecules(ParticleContainer* moleculeContainer, vector<unsigned long> &compCount) {
	vector<unsigned long> localCompCount(compCount.size(), 0);
	Molecule* tempMolecule;
	int temp = 0;
	for (tempMolecule = moleculeContainer->begin(); tempMolecule != moleculeContainer->end(); tempMolecule = moleculeContainer->next()) {
		localCompCount[tempMolecule->componentid()] += 1;
		temp++;
	}
	int numMolecules = 0;
	for (int i = 0; i < (int) localCompCount.size(); i++) {
		MPI_CHECK( MPI_Allreduce(&localCompCount[i], &compCount[i], 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD) );
		numMolecules += compCount[i];
	}
	return numMolecules;
}

double KDDecomposition::getBoundingBoxMin(int dimension, Domain* domain) {
	double globalLength = domain->getGlobalLength(dimension);
	double pos = (_ownArea->_lowCorner[dimension]) * _cellSize[dimension];
	if (pos < 0)
		return pos + globalLength;
	else if (pos > globalLength)
		return pos - globalLength;
	else
		return pos;
}

double KDDecomposition::getBoundingBoxMax(int dimension, Domain* domain) {
	double globalLength = domain->getGlobalLength(dimension);
	double pos = (_ownArea->_highCorner[dimension] + 1) * _cellSize[dimension];
	if (pos < 0)
		return pos + globalLength;
	else if (pos > globalLength)
		return pos - globalLength;
	else
		return pos;
}

void KDDecomposition::printDecomp(string filename, Domain* domain) {
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

void KDDecomposition::writeMoleculesToFile(string filename, ParticleContainer* moleculeContainer) {

	for (int process = 0; process < _numProcs; process++) {
		if (_ownRank == process) {
			ofstream checkpointfilestream(filename.c_str(), ios::app);
			checkpointfilestream.precision(20);
			Molecule* tempMolecule;
			for (tempMolecule = moleculeContainer->begin(); tempMolecule != moleculeContainer->end(); tempMolecule = moleculeContainer->next()) {
				tempMolecule->write(checkpointfilestream);
			}
			checkpointfilestream.close();
		}
		barrier();
	}
}

unsigned KDDecomposition::Ndistribution(unsigned localN, float* minrnd, float* maxrnd) {
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

void KDDecomposition::assertIntIdentity(int IX) {
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
				cout << "SEVERE ERROR: IX is " << IX << " for rank 0, but " << recv << " for rank " << i << ".\n";
				MPI_Abort(MPI_COMM_WORLD, 911);
			}
		}
		cout << "IX = " << recv << " for all " << num_procs << " ranks.\n";
	}
}

void KDDecomposition::assertDisjunctivity(TMoleculeContainer* mm) {
	Molecule* m;
	if (this->_ownRank) {
		unsigned long int tid;
		for (m = mm->begin(); m != mm->end(); m = mm->next()) {
			tid = m->id();
			MPI_CHECK( MPI_Send(&tid, 1, MPI_UNSIGNED_LONG, 0, 2674 + _ownRank, this->_collComm.getTopology()) );
		}
		// use ULONG_MAX to tell the receiving process that there are no more molecules
		tid = ULONG_MAX;
		MPI_CHECK( MPI_Send(&tid, 1, MPI_UNSIGNED_LONG, 0, 2674 + _ownRank, this->_collComm.getTopology()) );
	}
	else {
		unsigned long int recv;
		map<int, int> check;
		for (m = mm->begin(); m != mm->end(); m = mm->next()) {
			check[m->id()] = 0;
		}

		int num_procs;
		MPI_CHECK( MPI_Comm_size(this->_collComm.getTopology(), &num_procs) );
		MPI_Status s;
		for (int i = 1; i < num_procs; i++) {
			bool receiveMoreMolecules = true;
			while (receiveMoreMolecules) {
				MPI_CHECK( MPI_Recv(&recv, 1, MPI_UNSIGNED_LONG, i, 2674 + i, this->_collComm.getTopology(), &s) );
				if (recv == ULONG_MAX)
					receiveMoreMolecules = false;
				else {
					if (check.find(recv) != check.end()) {
						cout << "\nSEVERE ERROR. Ranks " << check[recv] << " and "
						     << i << " both propagate ID " << recv << ". Aborting.\n";
						cout << "recv=" << recv << endl;
						//MPI_Abort(MPI_COMM_WORLD, 2674);
					}
					else
						check[recv] = i;
				}
			}
		}
		cout << "\nData consistency checked. No duplicate IDs detected among " << check.size()
		     << " entries.\n";
	}
}


//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//$ private Methoden, die von exchangeMolecule benvtigt werden $
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

void KDDecomposition::getPartsToSend(KDNode* sourceArea, KDNode* decompTree, ParticleContainer* moleculeContainer, Domain* domain, vector<int>& procIDs, vector<int>& numMolsToSend, vector<list<Molecule*> >& particlesToSend) {
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

void KDDecomposition::exchangeNumToSend(vector<int>& procsToSendTo, vector<int>& procsToRecvFrom, vector<int>& numMolsToSend, vector<int>& numMolsToRecv) {
	numMolsToRecv.resize(procsToRecvFrom.size());
	MPI_Request request[procsToRecvFrom.size()];
	MPI_Status status;
	// initiate recv calls
	for (int neighbCount = 0; neighbCount < (int) procsToRecvFrom.size(); neighbCount++) {
		if (procsToRecvFrom[neighbCount] == _ownRank) continue; // don't exchange data with the own process
		MPI_CHECK( MPI_Irecv(&numMolsToRecv[neighbCount], 1, MPI_INT, procsToRecvFrom[neighbCount], 0, MPI_COMM_WORLD, &request[neighbCount]) );
	}

	// send all numbers
	for (int neighbCount = 0; neighbCount < (int) procsToSendTo.size(); neighbCount++) {
		if (procsToSendTo[neighbCount] == _ownRank) continue; // don't exchange data with the own process
		MPI_CHECK( MPI_Send(&numMolsToSend[neighbCount], 1, MPI_INT, procsToSendTo[neighbCount], 0, MPI_COMM_WORLD) );
	}
	// wait for the completion of all recv calls
	for (int neighbCount = 0; neighbCount < (int) procsToRecvFrom.size(); neighbCount++) {
		if (procsToRecvFrom[neighbCount] == _ownRank) continue; // don't exchange data with the own process
		MPI_CHECK( MPI_Wait(&request[neighbCount], &status) );
	}
}

void KDDecomposition::transferMolData(vector<int>& procsToSendTo, vector<int>& procsToRecvFrom, vector<int>& numMolsToSend, vector<int>& numMolsToRecv, vector<ParticleData*>& particlesSendBufs, vector<ParticleData*>& particlesRecvBufs) {

	// create a MPI Datatype which can store that molecule-data that has to be sent
	MPI_Datatype sendPartType;
	ParticleData::setMPIType(sendPartType);

	MPI_Request request[procsToRecvFrom.size()];
	MPI_Status status;

	// initiate recv calls
	for (int neighbCount = 0; neighbCount < (int) procsToRecvFrom.size(); neighbCount++) {
		if (procsToRecvFrom[neighbCount] == _ownRank) continue; // don't exchange data with the own process
		MPI_CHECK( MPI_Irecv(particlesRecvBufs[neighbCount], numMolsToRecv[neighbCount], sendPartType, procsToRecvFrom[neighbCount], 0, MPI_COMM_WORLD, &request[neighbCount]) );
	}
	// send all numbers
	for (int neighbCount = 0; neighbCount < (int) procsToSendTo.size(); neighbCount++) {
		if (procsToSendTo[neighbCount] == _ownRank) continue; // don't exchange data with the own process
		MPI_CHECK( MPI_Send(particlesSendBufs[neighbCount], numMolsToSend[neighbCount], sendPartType, procsToSendTo[neighbCount], 0, MPI_COMM_WORLD) );
	}
	// wait for the completion of all recv calls
	for (int neighbCount = 0; neighbCount < (int) procsToRecvFrom.size(); neighbCount++) {
		if (procsToRecvFrom[neighbCount] == _ownRank) continue; // don't exchange data with the own process
		MPI_CHECK( MPI_Wait(&request[neighbCount], &status) );
	}
}

void KDDecomposition::adjustOuterParticles(KDNode*& newOwnArea, ParticleContainer* moleculeContainer, Domain* domain) {

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

void KDDecomposition::createLocalCopies(ParticleContainer* moleculeContainer, Domain* domain, const vector<Component>& components) {
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
				Molecule m1 = Molecule(molPtr->id(),molPtr->componentid(),
				                       newPosition[0], newPosition[1], newPosition[2],
				                       molPtr->v(0),molPtr->v(1),molPtr->v(2),
				                       molPtr->q().qw(),molPtr->q().qx(),molPtr->q().qy(),molPtr->q().qz(),
				                       molPtr->D(0),molPtr->D(1),molPtr->D(2), &components);
				moleculeContainer->addParticle(m1);
			}
			molPtr = moleculeContainer->next();
		}
	}
}

void KDDecomposition::recInitialDecomp(KDNode* fatherNode, KDNode*& ownArea) {
	bool coversAll[KDDIM];
	int cellsPerDim[KDDIM];
	for (int dim = 0; dim < KDDIM; dim++) {
		coversAll[dim] = fatherNode->_coversWholeDomain[dim];
		cellsPerDim[dim] = fatherNode->_highCorner[dim] - fatherNode->_lowCorner[dim] + 1;
	}
	int divDir = 0;
	int maxCells = cellsPerDim[0];
	if (cellsPerDim[1] > maxCells) {
		divDir = 1;
		maxCells = cellsPerDim[1];
	}
	if (cellsPerDim[2] > maxCells) {
		divDir = 2;
		maxCells = cellsPerDim[2];
	}

	coversAll[divDir] = false;

	int low1[KDDIM];
	int low2[KDDIM];
	int high1[KDDIM];
	int high2[KDDIM];
	int numProcs1;
	int numProcs2;
	int id1;
	int id2;
	int owner1;
	int owner2;

	for (int dim = 0; dim < KDDIM; dim++) {
		low1[dim] = fatherNode->_lowCorner[dim];
		low2[dim] = fatherNode->_lowCorner[dim];
		high1[dim] = fatherNode->_highCorner[dim];
		high2[dim] = fatherNode->_highCorner[dim];
	}

	low1[divDir] = fatherNode->_lowCorner[divDir];
	high1[divDir] = (fatherNode->_highCorner[divDir] + fatherNode->_lowCorner[divDir]) / 2;
	low2[divDir] = high1[divDir] + 1;
	high2[divDir] = fatherNode->_highCorner[divDir];

	numProcs1 = fatherNode->_numProcs / 2;
	numProcs2 = fatherNode->_numProcs - numProcs1;
	id1 = fatherNode->_nodeID + 1;
	id2 = fatherNode->_nodeID + 2 * numProcs1;
	owner1 = fatherNode->_owningProc;
	owner2 = owner1 + numProcs1;
	fatherNode->_child1 = new KDNode(numProcs1, low1, high1, id1, owner1, coversAll);
	fatherNode->_child2 = new KDNode(numProcs2, low2, high2, id2, owner2, coversAll);

	if (numProcs1 > 1)
		recInitialDecomp(fatherNode->_child1, ownArea);
	else if (owner1 == _ownRank)
		ownArea = fatherNode->_child1;

	if (numProcs2 > 1)
		recInitialDecomp(fatherNode->_child2, ownArea);
	else if (owner2 == _ownRank)
		ownArea = fatherNode->_child2;
}

bool KDDecomposition::recDecompPar(KDNode* fatherNode, KDNode*& ownArea, MPI_Comm commGroup) {
	bool domainTooSmall = false;
	// This happens if the whole simulation is run with a single process
	if (fatherNode->_numProcs == 1) {
		if (fatherNode->_owningProc == _ownRank) {
			ownArea = fatherNode;
			return domainTooSmall;
		}
		else {
			cerr << "ERROR in recDecompPar: called with a leaf node" << endl;
			exit(1);
		}
	}
	bool coversAll[KDDIM];
    
    /* TODO: We do not use this values anywhere ... */
	int cellsPerDim[KDDIM];
	for (int dim = 0; dim < KDDIM; dim++) {
		coversAll[dim] = fatherNode->_coversWholeDomain[dim];
		cellsPerDim[dim] = fatherNode->_highCorner[dim] - fatherNode->_lowCorner[dim] + 1;
	}

	int divDir = 0;
	int divIdx = 0;
	double maxProcCost = DBL_MAX;
	int numProcsLeft = 0;
	int numProcsRight = 0;

	vector<vector<double> > costsLeft(3);
	vector<vector<double> > costsRight(3);
	calculateCostsPar(fatherNode, costsLeft, costsRight, true, commGroup);

	for (int dim = 0; dim < 3; dim++) {

		// loop only from 1 to max-1 (instead 0 to max) to avoid 1-cell-regions
		for (int i = 1; i < fatherNode->_highCorner[dim] - fatherNode->_lowCorner[dim] - 1; i++) {
			//#################################
			// New Version: Costs for "future" boundary cells
			// Num Cells in separation layer
			unsigned long numCellsInLayer = 1;
			for (int temp = 0; temp < 3; temp++) {
				if (temp != dim) {
					numCellsInLayer *= fatherNode->_highCorner[temp] - fatherNode->_lowCorner[temp] + 1;
				}
			}
			unsigned long numCellsLeft = numCellsInLayer * (i + 1);
			unsigned long numCellsRight = numCellsInLayer * (fatherNode->_highCorner[dim] - fatherNode->_lowCorner[dim] - i);

			double optCostPerProc = (costsLeft[dim][i] + costsRight[dim][i]) / ((double) fatherNode->_numProcs);

			// costs for initial guessed process distribution
			int numProcsLeftTest = (int) max(1.0, min(floor(costsLeft[dim][i] / optCostPerProc), (double) (fatherNode->_numProcs - 1)));
			int numProcsRightTest = fatherNode->_numProcs - numProcsLeftTest;
			// Ensure that each sub-area is large enough to be distributed to the chosen number of processes
			bool ok = false;
			bool procShiftAllowed = true;
			int maxProcsLeft = (i + 1) / 2;
			int maxProcsRight = (fatherNode->_highCorner[dim] - fatherNode->_lowCorner[dim] - i) / 2;
			for (int innerDim = 0; innerDim < 3; innerDim++) {
				if (innerDim != dim) {
					maxProcsLeft *= (fatherNode->_highCorner[innerDim] - fatherNode->_lowCorner[innerDim] + 1) / 2;
					maxProcsRight *= (fatherNode->_highCorner[innerDim] - fatherNode->_lowCorner[innerDim] + 1) / 2;
				}
			}
			// if not resolvable, continue with next loop
			if ((numProcsLeftTest + numProcsRightTest) > (maxProcsLeft + maxProcsRight)) {
				domainTooSmall = true;
				continue;
			}
			if (numProcsLeftTest <= maxProcsLeft && numProcsRightTest <= maxProcsRight)
				ok = true;
			if (not (ok)) {
				domainTooSmall = true;
				procShiftAllowed = false;
				int shiftDir;
				if (maxProcsLeft < numProcsLeftTest)
					shiftDir = 1; // move procs to right
				else
					shiftDir = -1; // move procs to left
				while (not (ok)) {
					numProcsLeftTest -= shiftDir;
					numProcsRightTest += shiftDir;
					if (maxProcsLeft >= numProcsLeftTest && maxProcsRight >= numProcsRightTest)
						ok = true;
				}
			}
			unsigned long numBoundCellsLeft = numCellsLeft - (unsigned long) (((double) numProcsLeftTest) * pow(pow(numCellsLeft / numProcsLeftTest, 1. / 3.) - 2, 3));
			unsigned long numBoundCellsRight = numCellsRight - (unsigned long) (((double) numProcsRightTest) * pow(pow(numCellsRight / numProcsRightTest, 1. / 3.) - 2, 3));
			double maxProcCostOld = max(costsLeft[dim][i] / (double) numProcsLeftTest, costsRight[dim][i] / (double) numProcsRightTest);
			// Find out in which direction process distribution has to be shifted
			int procShift = 1;
			if (costsLeft[dim][i] / numProcsLeftTest > costsRight[dim][i] / numProcsRightTest) {
				procShift = -1;
			}
			double maxProcCostNew = maxProcCostOld;
			// if shifted distribution better the old one, continue shifting
			if (procShiftAllowed) {
				while (maxProcCostNew <= maxProcCostOld) {
					maxProcCostOld = maxProcCostNew;
					numProcsLeftTest = numProcsLeftTest + procShift;
					numProcsRightTest = numProcsRightTest - procShift;
					if (numProcsLeftTest > maxProcsLeft || numProcsRightTest > maxProcsRight) {
						break;
					}
					if (numProcsLeftTest >= fatherNode->_numProcs || numProcsLeftTest <= 0 || numProcsRightTest >= fatherNode->_numProcs || numProcsRightTest <= 0) {
						break;
					}
					else {
						numBoundCellsLeft = numCellsLeft - (unsigned long) (((double) numProcsLeftTest) * pow(pow(numCellsLeft / numProcsLeftTest, 1. / 3.) - 2, 3));
						numBoundCellsRight = numCellsRight - (unsigned long) (((double) numProcsRightTest) * pow(pow(numCellsRight / numProcsRightTest, 1. / 3.) - 2, 3));
						maxProcCostNew = max(costsLeft[dim][i] / (double) numProcsLeftTest, costsRight[dim][i] / (double) numProcsRightTest);
					}
				}
				// shift back one element, as old dist was better
				numProcsLeftTest = numProcsLeftTest - procShift;
				numProcsRightTest = numProcsRightTest + procShift;
			}

			if (maxProcCostOld < maxProcCost) {
				divDir = dim;
				divIdx = i + fatherNode->_lowCorner[dim];
				maxProcCost = maxProcCostOld;
				numProcsLeft = numProcsLeftTest;
				numProcsRight = fatherNode->_numProcs - numProcsLeftTest;
			}
			// End new version

			if (numProcsLeft <= 0 || numProcsLeft >= fatherNode->_numProcs) {
				if (_ownRank == 0) {
					cerr << "ERROR in recDecompPar, part of the domain was not assigned to a proc" << endl;
				}
				exit(1);
			}
		}
	}

	coversAll[divDir] = false;

	int low1[KDDIM];
	int low2[KDDIM];
	int high1[KDDIM];
	int high2[KDDIM];
	int id1;
	int id2;
	int owner1;
	int owner2;

	for (int dim = 0; dim < KDDIM; dim++) {
		low1[dim] = fatherNode->_lowCorner[dim];
		low2[dim] = fatherNode->_lowCorner[dim];
		high1[dim] = fatherNode->_highCorner[dim];
		high2[dim] = fatherNode->_highCorner[dim];
	}

	low1[divDir] = fatherNode->_lowCorner[divDir];
	high1[divDir] = divIdx;
	low2[divDir] = high1[divDir] + 1;
	high2[divDir] = fatherNode->_highCorner[divDir];

	id1 = fatherNode->_nodeID + 1;
	id2 = fatherNode->_nodeID + 2 * numProcsLeft;
	owner1 = fatherNode->_owningProc;
	owner2 = owner1 + numProcsLeft;

	fatherNode->_child1 = new KDNode(numProcsLeft, low1, high1, id1, owner1, coversAll);
	fatherNode->_child2 = new KDNode(numProcsRight, low2, high2, id2, owner2, coversAll);

	vector<int> origRanks;
	int newNumProcs;
	if (_ownRank < owner2) {
		origRanks.resize(numProcsLeft);
		for (int i = 0; i < numProcsLeft; i++) {
			//origRanks[i] = i + owner1;
			origRanks[i] = i;
		}
		newNumProcs = numProcsLeft;
	}
	else {
		origRanks.resize(numProcsRight);
		for (int i = 0; i < numProcsRight; i++) {
			origRanks[i] = i + numProcsLeft;
		}
		newNumProcs = numProcsRight;
	}

	MPI_Comm newComm;
	MPI_Group origGroup, newGroup;

	MPI_CHECK( MPI_Comm_group(commGroup, &origGroup) );
	MPI_CHECK( MPI_Group_incl(origGroup, newNumProcs, &origRanks[0], &newGroup) );
	MPI_CHECK( MPI_Comm_create(commGroup, newGroup, &newComm) );

	if (numProcsLeft > 1 && _ownRank < owner2) {
		domainTooSmall = (domainTooSmall || recDecompPar(fatherNode->_child1, ownArea, newComm));
	}
	else if (owner1 == _ownRank)
		ownArea = fatherNode->_child1;

	if (numProcsRight > 1 && _ownRank >= owner2) {
		domainTooSmall = (domainTooSmall || recDecompPar(fatherNode->_child2, ownArea, newComm));
	}
	else if (owner2 == _ownRank)
		ownArea = fatherNode->_child2;

	return domainTooSmall;

}

void KDDecomposition::completeTreeInfo(KDNode*& root, KDNode*& ownArea) {

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
		//pack:
		vector<int> data;
		data.resize(13, 0);
		if (oldNode->_nodeID == nodeID) {

			data[0] = oldNode->_numProcs;
			for (int i = 0; i < 3; i++) {
				data[1 + i] = oldNode->_lowCorner[i];
				data[4 + i] = oldNode->_highCorner[i];
			}
			if (oldNode->_coversWholeDomain[0]) data[7] += 4;
			if (oldNode->_coversWholeDomain[1]) data[7] += 2;
			if (oldNode->_coversWholeDomain[2]) data[7] += 1;
			data[8] = oldNode->_nodeID;
			data[9] = oldNode->_owningProc;
			// determine owner of next data set and IDs of children (-1 of no children)
			if (oldNode->_numProcs > 1) {
				data[10] = oldNode->_child1->_nodeID;
				data[11] = oldNode->_child2->_nodeID;
				data[12] = _ownRank;
				// forward pointer to children
				oldNode = oldNode->_child1;
			}
			else {
				data[10] = -1;
				data[11] = -1;
				data[12] = _ownRank + 1;
			}

		}
		MPI_CHECK( MPI_Bcast(&data[0], 13, MPI_INT, nextSendingProcess, MPI_COMM_WORLD) );
		bool coversAll[3];
		for (int i = 0; i < 3; i++)
			coversAll[i] = false;
		if (data[7] >= 4) {
			coversAll[0] = true;
			data[7] -= 4;
		}
		if (data[7] >= 2) {
			coversAll[0] = true;
			data[7] -= 2;
		}
		if (data[7] >= 1) {
			coversAll[0] = true;
		}
		ptrToAllNodes[nodeID] = new KDNode(data[0], &data[1], &data[4], data[8], data[9], coversAll);
		child1[nodeID] = data[10];
		child2[nodeID] = data[11];
		nextSendingProcess = data[12];
	}

	// connect all Nodes of the tree
	for (int nodeID = 0; nodeID < numElementsToRecv; nodeID++) {
		// if node is not leaf node, connect to the two children
		if (child1[nodeID] >= 0) {
			ptrToAllNodes[nodeID]->_child1 = ptrToAllNodes[child1[nodeID]];
			ptrToAllNodes[nodeID]->_child2 = ptrToAllNodes[child2[nodeID]];
		}
	}

	//@TODO: the old root has to be deleted
	//KDNode* oldRoot = root;
	root = ptrToAllNodes[0];
	ownArea = ptrToAllNodes[ownArea->_nodeID];

	//@TODO: the old root has to be deleted
	//delete oldRoot;

}

void KDDecomposition::printDecompTree(KDNode* root, string prefix) {
	if (root->_numProcs == 1) {
		cout << prefix << "LEAF: " << root->_nodeID << ", Owner: " << root->_owningProc << ", Corners: (" << root->_lowCorner[0] << ", " << root->_lowCorner[1] << ", " << root->_lowCorner[2] << ") / (" << root->_highCorner[0] << ", " << root->_highCorner[1] << ", " << root->_highCorner[2] << ")" << endl;
	}
	else {
		cout << prefix << "INNER: " << root->_nodeID << ", Owner: " << root->_owningProc << "(" << root->_numProcs << " procs)" << ", Corners: (" << root->_lowCorner[0] << ", " << root->_lowCorner[1] << ", " << root->_lowCorner[2] << ") / (" << root->_highCorner[0] << ", " << root->_highCorner[1] << ", " << root->_highCorner[2] << ")" << endl;
		stringstream childprefix;
		childprefix << prefix << "  ";
		printDecompTree(root->_child1, childprefix.str());
		printDecompTree(root->_child2, childprefix.str());
	}
}

void KDDecomposition::calculateCostsPar(KDNode* area, vector<vector<double> >& costsLeft, vector<vector<double> >& costsRight, bool calcDivisionCosts, MPI_Comm commGroup) {

	double areaCosts = 0.0;

	vector<vector<double> > cellCosts;
	vector<vector<double> > sepCostLeft;
	vector<vector<double> > sepCostRight;
	cellCosts.resize(3);
	sepCostLeft.resize(3);
	sepCostRight.resize(3);

	int newRank;
	MPI_Group group;

	MPI_CHECK( MPI_Comm_group(commGroup, &group) );
	MPI_CHECK( MPI_Group_rank(group, &newRank) );

	vector<vector<double> > costsLeftTemp;
	vector<vector<double> > costsRightTemp;
	costsLeftTemp.resize(3);
	costsRightTemp.resize(3);
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
		// If only total costs shall be calculated, and not the different possible divisions,
		// the loop runs only once
		if (not calcDivisionCosts && dim != 0) break;

		cellCosts[dim].resize(area->_highCorner[dim] - area->_lowCorner[dim] + 1, 0.0);
		sepCostLeft[dim].resize(area->_highCorner[dim] - area->_lowCorner[dim] + 1, 0.0);
		sepCostRight[dim].resize(area->_highCorner[dim] - area->_lowCorner[dim] + 1, 0.0);

		costsLeftTemp[dim].resize(area->_highCorner[dim] - area->_lowCorner[dim] + 1, 0.0);
		costsRightTemp[dim].resize(area->_highCorner[dim] - area->_lowCorner[dim] + 1, 0.0);
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

			sepCostLeft[dim][i_dim] = 0;
			sepCostRight[dim][i_dim] = 0;
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
					float numForcePairs = _globalLoadPerCell[getGlobalIndex(dim, dim1, dim2, i_dim, i_dim1, i_dim2, area)];

					// #######################
					// ## Cell Costs        ##
					// #######################
					cellCosts[dim][i_dim] += _alpha * (double) (numParts * numParts);
					cellCosts[dim][i_dim] += (1.0 - _alpha) * (double) numForcePairs;

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
								if (getGlobalIndex(dim, dim1, dim2, nI_dim, nI_dim1, nI_dim2, area) > getGlobalIndex(dim, dim1, dim2, i_dim, i_dim1, i_dim2, area)) {
									int numPartsNeigh = (int) _numParticlesPerCell[getGlobalIndex(dim, dim1, dim2, nI_dim, nI_dim1, nI_dim2, area)];
									cellCosts[dim][i_dim] += _alpha * (double) (numParts * numPartsNeigh);
								}
							}
						}
					}

					// #######################
					// ## Separation Costs  ##
					// #######################

					// Neighbours on the other (right) side of the dividing plane (consider periodic boundary)
					int neighInd_divDim = i_dim + 1;
					if (neighInd_divDim + area->_lowCorner[dim] == _globalCellsPerDim[dim])
						neighInd_divDim = 0;

					// loop over all nine neighbours right of that plane
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
							int numPartsNeigh = (int) _numParticlesPerCell[getGlobalIndex(dim, dim1, dim2, neighInd_divDim, nI_dim1, nI_dim2, area)];
							if (getGlobalIndex(dim, dim1, dim2, neighInd_divDim, nI_dim1, nI_dim2, area) > getGlobalIndex(dim, dim1, dim2, i_dim, i_dim1, i_dim2, area)) {
								sepCostRight[dim][i_dim] += _alpha * (double) (numParts * numPartsNeigh);
							}
							else {
								sepCostLeft[dim][i_dim] += _alpha * (double) (numParts * numPartsNeigh);
							}

							// beta times direct neighbour
							if (neighInd_dim1 == i_dim1 && neighInd_dim2 == i_dim2) {
								sepCostLeft[dim][i_dim] += _beta * (double) (numParts + numPartsNeigh);
								sepCostRight[dim][i_dim] += _beta * (double) (numParts + numPartsNeigh);
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
		vector<double> sepCostSumLeft;
		vector<double> sepCostSumRight;
		cellCostsSum.resize(area->_highCorner[dim] - area->_lowCorner[dim] + 1, 0.0);
		sepCostSumLeft.resize(area->_highCorner[dim] - area->_lowCorner[dim] + 1, 0.0);
		sepCostSumRight.resize(area->_highCorner[dim] - area->_lowCorner[dim] + 1, 0.0);

		int size2 = cellCostsSum.size();
		MPI_CHECK( MPI_Allreduce(&cellCosts[dim][0], &cellCostsSum[0], size2, MPI_DOUBLE, MPI_SUM, commGroup) );
		MPI_CHECK( MPI_Allreduce(&sepCostLeft[dim][0], &sepCostSumLeft[0], size2, MPI_DOUBLE, MPI_SUM, commGroup) );
		MPI_CHECK( MPI_Allreduce(&sepCostRight[dim][0], &sepCostSumRight[0], size2, MPI_DOUBLE, MPI_SUM, commGroup) );

		for (int i_dim = 0; i_dim <= area->_highCorner[dim] - area->_lowCorner[dim]; i_dim++) {
			areaCosts += cellCosts[dim][i_dim];
			// Num Cells in separation layer
			unsigned long sepCells = 1;
			for (int i = 0; i < 3; i++) {
				if (i != dim) {
					sepCells *= area->_highCorner[i] - area->_lowCorner[i] + 1;
				}
			}
			// guessed costs for future boundary cells
			costsLeft[dim][i_dim] = cellCostsSum[i_dim] + sepCostSumLeft[i_dim];
			costsRight[dim][i_dim] = cellCostsSum[area->_highCorner[dim] - area->_lowCorner[dim]] - cellCostsSum[i_dim] + sepCostSumRight[i_dim];
		}
	}

}

unsigned int KDDecomposition::getGlobalIndex(int divDim, int dim1, int dim2, int index_divDim, int index_dim1, int index_dim2, KDNode* localArea) {
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

int KDDecomposition::mod(int number, int modulo) {
	int result = number % modulo;
	if (result < 0)
		result += modulo;
	return result;
}

//##########################################################################
//##########################################################################
//###                                                                    ###
//###                  evtl unnoetige private Methoden                   ###
//###                                                                    ###
//##########################################################################
//##########################################################################

void KDDecomposition::getOwningProcs(int low[KDDIM], int high[KDDIM], KDNode* decompTree, KDNode* testNode, vector<int>* procIDs, vector<int>* neighbHaloAreas) {
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

void KDDecomposition::getNumParticles(ParticleContainer* moleculeContainer) {
	for (int i = 0; i < _globalNumCells; i++)
		_numParticlesPerCell[i] = 0;

	int count = 0;
	double bBMin[3]; // haloBoundingBoxMin
    /* TODO: We do not use values form bBMax anywhere ... */
	double bBMax[3]; // haloBoundingBoxMax
	for (int dim = 0; dim < 3; dim++) {
		bBMin[dim] = moleculeContainer->getBoundingBoxMin(dim);// - moleculeContainer->get_halo_L(dim);
		bBMax[dim] = moleculeContainer->getBoundingBoxMax(dim);// + moleculeContainer->get_halo_L(dim);
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
	// TODO
	// memory problem:
	// Some MPI-Implementations demand that in a Allreduce command, send and recv buffers
	// must be different. The consequence is, that a temporary array has to be created for
	// the values to be recieved, which that has to be copied back to the original array
	unsigned char* numParticlesPerCellTemp = new unsigned char[_globalNumCells];
	MPI_CHECK( MPI_Allreduce(_numParticlesPerCell, numParticlesPerCellTemp, _globalNumCells, MPI_UNSIGNED_CHAR, MPI_SUM, MPI_COMM_WORLD) );
	delete[] _numParticlesPerCell;
	_numParticlesPerCell = numParticlesPerCellTemp;

}

