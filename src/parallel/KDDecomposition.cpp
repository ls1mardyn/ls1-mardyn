#include "parallel/KDDecomposition.h" //$ Pfad dndern (parallel/...)

#include <sstream>
#include "molecules/Molecule.h"
#include "Domain.h"
#include "datastructures/ParticleContainer.h"
#include "parallel/ParticleData.h"
#include "parallel/KDNode.h"
#include <float.h>

KDDecomposition::KDDecomposition(int *argc, char ***argv, double cutoffRadius, Domain* domain){
  int procnamelen;

  MPI_Comm_rank(MPI_COMM_WORLD, &_ownRank);
  MPI_Comm_size(MPI_COMM_WORLD, &_numProcs);
  MPI_Get_processor_name(_processorName, &procnamelen);

  int lowCorner[KDDIM];
  int highCorner[KDDIM];
  bool coversWholeDomain[KDDIM];

  for(int dim=0; dim<KDDIM; dim++){
    _globalCellsPerDim[dim] = (int) floor(domain->getGlobalLength(dim)/cutoffRadius);
    lowCorner[dim] = 0;
    highCorner[dim] = _globalCellsPerDim[dim] - 1;
    _cellSize[dim] = domain->getGlobalLength(dim)/((double)_globalCellsPerDim[dim]);
    coversWholeDomain[dim] = true;
  }
  _globalNumCells = (_globalCellsPerDim[0])*(_globalCellsPerDim[1])*(_globalCellsPerDim[2]);
  _numParticlesPerCell = new int[_globalNumCells];
  for(int i=0; i<_globalNumCells; i++) _numParticlesPerCell[i] = 0;

  // create initial decomposition
  _decompTree = new KDNode(_numProcs, lowCorner, highCorner, 0, 0, coversWholeDomain);
  if(_numProcs>1) recInitialDecomp(_decompTree, _ownArea);
  else _ownArea = _decompTree;

  _decompValid = true;

}

KDDecomposition::~KDDecomposition(){
  delete []_numParticlesPerCell;
  MPI_Finalize();
}

void KDDecomposition::exchangeMolecules(ParticleContainer* moleculeContainer, const vector<Component>& components, Domain* domain, double rc){
  balanceAndExchange(false, moleculeContainer, components, domain, rc);
}

void KDDecomposition::balanceAndExchange(bool balance, ParticleContainer* moleculeContainer, const vector<Component>& components, Domain* domain, double rc){
  _decompValid=true;
  KDNode* newDecompTree = NULL;
  KDNode* newOwnArea = NULL;
  if(balance){
    getNumParticles(moleculeContainer);
    newDecompTree = new KDNode(_numProcs, _decompTree->_lowCorner, _decompTree->_highCorner, 0, 0, _decompTree->_coversWholeDomain);
    _decompValid = false;
    recDecomp(newDecompTree, newOwnArea);
  }
  vector<int> procsToSendTo; // all processes to which this process has to send data
  vector<int> procsToRecvFrom; // all processes from which this process has to recv data
  vector<list<Molecule*> > particlePtrsToSend; // pointer to particles to be send
  //vector<struct ParticleData*> particlesSendBufs; // buffer used by my send call
  //vector<struct ParticleData*> particlesRecvBufs; // buffer used by my recv call
  vector<ParticleData*> particlesSendBufs; // buffer used by my send call
  vector<ParticleData*> particlesRecvBufs; // buffer used by my recv call
  vector<int> numMolsToSend; // number of particles to be send to other procs
  vector<int> numMolsToRecv; // number of particles to be recieved from other procs
  // collect particles to be send and find out number of particles to be recieved
  if(balance) {
    int haloCellIdxMin[3]; // Assuming a global 3D Cell index, haloCellIdxMin[3] gives the position
    // of the low local domain corner within this global 3D cell index
    int haloCellIdxMax[3]; // same as heloCellIdxMax, only high instead of low Corner
    for(int dim=0; dim<3; dim++){
      haloCellIdxMin[dim] = newOwnArea->_lowCorner[dim]-1;
      haloCellIdxMax[dim] = newOwnArea->_highCorner[dim]+1;
    }
    vector<int> neighbHaloAreas; // The areas (unit: cells, including halo) of the neighbouring procs
    // For each proc, 6 int values are reserved (xlow, ylow, zlow, xhigh,...)
    // These values are not used in this context.
    getOwningProcs(haloCellIdxMin, haloCellIdxMax, _decompTree, _decompTree, &procsToRecvFrom, &neighbHaloAreas);
    getPartsToSend(_ownArea, newDecompTree, moleculeContainer, domain, procsToSendTo, numMolsToSend, particlePtrsToSend);
  }
  else{
    getPartsToSend(_ownArea, _decompTree, moleculeContainer, domain, procsToSendTo, numMolsToSend, particlePtrsToSend);
    procsToRecvFrom = procsToSendTo;
  }
  numMolsToRecv.resize(procsToRecvFrom.size());
  exchangeNumToSend(procsToSendTo, procsToRecvFrom, numMolsToSend, numMolsToRecv);
  // Initialise send and recieve buffers
  particlesSendBufs.resize(procsToSendTo.size());
  particlesRecvBufs.resize(procsToRecvFrom.size());
  for(int neighbCount=0; neighbCount<(int)procsToSendTo.size();neighbCount++) {
    if(procsToSendTo[neighbCount]==_ownRank) continue; // don't exchange data with the own process
    // TODO: replace particlePtrsToSend[neighbCount].size() by numMolsToSend[neighbCount]
    particlesSendBufs[neighbCount] = new ParticleData[particlePtrsToSend[neighbCount].size()];
  }
  for(int neighbCount=0; neighbCount<(int)procsToRecvFrom.size();neighbCount++) {
    if(procsToRecvFrom[neighbCount]==_ownRank) continue; // don't exchange data with the own process
    particlesRecvBufs[neighbCount] = new ParticleData[numMolsToRecv[neighbCount]];
  }

  // Fill send buffer with particle data
  for(int neighbCount = 0; neighbCount < (int) procsToSendTo.size(); neighbCount++ ) {
    if(procsToSendTo[neighbCount]==_ownRank) continue; // don't exchange data with the own process
    list<Molecule*>::iterator particleIter;
    int partCount = 0;
    for( particleIter =  particlePtrsToSend[neighbCount].begin(); particleIter != particlePtrsToSend[neighbCount].end(); particleIter++ ){
      ParticleData::MoleculeToParticleData( particlesSendBufs[neighbCount][partCount], **particleIter );
      partCount++;
    }
  }
  if(balance){
    // find out new bounding boxes (of newOwnArea)
    double bBoxMin[3];
    double bBoxMax[3];
    for(int dim=0; dim<3; dim++){
      bBoxMin[dim] = (newOwnArea->_lowCorner[dim])*_cellSize[dim];
      bBoxMax[dim] = (newOwnArea->_highCorner[dim]+1)*_cellSize[dim];
    }
    // shift the region of the moleculeContainer and delete all particles which are no
    // longer in the new region
    // TODO Rebuild problem with particles leaving the domain
    moleculeContainer->rebuild(bBoxMin, bBoxMax);

    _ownArea = newOwnArea;
    _decompTree = newDecompTree;
    _decompValid = true;
  }
  // Transfer data

  transferMolData(procsToSendTo, procsToRecvFrom, numMolsToSend, numMolsToRecv, particlesSendBufs, particlesRecvBufs);

  double lowLimit[3];
  double highLimit[3];
  for(int dim=0; dim<3; dim++){
    lowLimit[dim] = moleculeContainer->getBoundingBoxMin(dim) - moleculeContainer->get_halo_L(dim);
    highLimit[dim] = moleculeContainer->getBoundingBoxMax(dim) + moleculeContainer->get_halo_L(dim);
  }
  // store recieved molecules in the molecule container
  ParticleData newMol;
  for(int neighbCount=0;neighbCount<(int)procsToRecvFrom.size();neighbCount++) {
    if(procsToRecvFrom[neighbCount]==_ownRank) continue; // don't exchange data with the own process
    for(int i=0; i<numMolsToRecv[neighbCount]; i++){
    	newMol = particlesRecvBufs[neighbCount][i];
      // change coordinates (especially needed if particle was moved across boundary)
      if(newMol.r[0] < lowLimit[0]) newMol.r[0] += domain->getGlobalLength(0);
      else if(newMol.r[0] >= highLimit[0]) newMol.r[0] -= domain->getGlobalLength(0);
      if(newMol.r[1] < lowLimit[1]) newMol.r[1] += domain->getGlobalLength(1);
      else if(newMol.r[1] >= highLimit[1]) newMol.r[1] -= domain->getGlobalLength(1);
      if(newMol.r[2] < lowLimit[2]) newMol.r[2] += domain->getGlobalLength(2);
      else if(newMol.r[2] >= highLimit[2]) newMol.r[2] -= domain->getGlobalLength(2);

      Molecule m1 = Molecule( newMol.id, newMol.cid, newMol.r[0], newMol.r[1], newMol.r[2], newMol.v[0], newMol.v[1], newMol.v[2],
                              newMol.q[0], newMol.q[1], newMol.q[2], newMol.q[3], newMol.D[0], newMol.D[1], newMol.D[2], &components );
      moleculeContainer->addParticle(m1);

    }
  }
  // create the copies of local molecules due to periodic boundaries
  // (only for procs covering the whole domain in one dimension)
  createLocalCopies(moleculeContainer, domain, components);

  // free memory of send and recv buffers
  for(int neighbCount=0; neighbCount<(int)procsToSendTo.size();neighbCount++) {
    if(procsToSendTo[neighbCount]==_ownRank) continue; // don't exchange data with the own process
    delete [] particlesSendBufs[neighbCount];
  }
  for(int neighbCount=0; neighbCount<(int)procsToRecvFrom.size();neighbCount++) {
    if(procsToRecvFrom[neighbCount]==_ownRank) continue; // don't exchange data with the own process
    delete [] particlesRecvBufs[neighbCount];
  }
  if(balance){
    moleculeContainer->update();
  }
  particlesSendBufs.resize(0);
  particlesRecvBufs.resize(0);
}

bool KDDecomposition::procOwnsPos(double x, double y, double z, Domain* domain){
  if(x < getBoundingBoxMin(0, domain) || x >= getBoundingBoxMax(0, domain)){
    return false;
  }
  else if(y < getBoundingBoxMin(1, domain) || y >= getBoundingBoxMax(1, domain)){
    return false;
  }
  else if(z < getBoundingBoxMin(2, domain) || z >= getBoundingBoxMax(2, domain)){
    return false;
  }
  else {
    return true;
  }
}

double KDDecomposition::guaranteedDistance(double x, double y, double z, Domain* domain){
  double xdist = 0;
  double ydist = 0;
  double zdist = 0;
  if (x < getBoundingBoxMin(0, domain)){
    xdist = getBoundingBoxMin(0, domain) - x;
  }
  else if(x >= getBoundingBoxMax(0, domain)){
    xdist = x - getBoundingBoxMin(0, domain);
  }
  if (y < getBoundingBoxMin(1, domain)){
    ydist = getBoundingBoxMin(1, domain) - y;
  }
  else if(y >= getBoundingBoxMax(1, domain)){
    ydist = y - getBoundingBoxMin(1, domain);
  }
  if (z < getBoundingBoxMin(2, domain)){
    zdist = getBoundingBoxMin(2, domain) - z;
  }
  else if(z >= getBoundingBoxMax(2, domain)){
    zdist = z - getBoundingBoxMin(2, domain);
  }
  return sqrt(xdist*xdist + ydist*ydist + zdist*zdist);
}

int KDDecomposition::countMolecules(ParticleContainer* moleculeContainer, vector<int> &compCount){
  vector<int> localCompCount;
  localCompCount.resize(compCount.size());
  for(int i=0; i<(int)localCompCount.size(); i++){
    localCompCount[i] = 0;
  }
  Molecule* tempMolecule;
  for(tempMolecule = moleculeContainer->begin(); tempMolecule != moleculeContainer->end(); tempMolecule = moleculeContainer->next()){
    localCompCount[tempMolecule->componentid()] += 1;
  }
  int numMolecules = 0;
  for(int i=0; i<(int)localCompCount.size(); i++){
    MPI_Allreduce(&localCompCount[i], &compCount[i], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    numMolecules += compCount[i];
  }
  return numMolecules;
}

double KDDecomposition::getBoundingBoxMin(int dimension, Domain* domain){
  if(_decompValid==false){
    cerr << "KDDecomposition::getBoundingBoxMin failed. Decomposition is invalid" << endl;
    exit(1);
  }
  double globalLength = domain->getGlobalLength(dimension);
  double pos = (_ownArea->_lowCorner[dimension])*_cellSize[dimension];
  if(pos < 0) return pos + globalLength;
  else if(pos>globalLength) return pos - globalLength;
  else return pos;
}

double KDDecomposition::getBoundingBoxMax(int dimension, Domain* domain){
  if(_decompValid==false){
    cerr << "KDDecomposition::getBoundingBoxMax failed. Decomposition is invalid" << endl;
    exit(1);
  }
  double globalLength = domain->getGlobalLength(dimension);
  double pos = (_ownArea->_highCorner[dimension] + 1)*_cellSize[dimension];
  if(pos < 0) return pos + globalLength;
  else if(pos>globalLength) return pos - globalLength;
  else return pos;
}

void KDDecomposition::printDecomp(string filename, Domain* domain){

  if(_ownRank==0) {
    ofstream povcfgstrm(filename.c_str());
    povcfgstrm << "size " << domain->getGlobalLength(0) << " " << domain->getGlobalLength(1) << " " << domain->getGlobalLength(2) << endl;
    povcfgstrm << "decompData Regions" << endl;
    povcfgstrm.close();
  }

  for(int process = 0; process < _numProcs; process++){
    if(_ownRank==process){
      ofstream povcfgstrm(filename.c_str(), ios::app);
      povcfgstrm << getBoundingBoxMin(0,domain) << " " << getBoundingBoxMin(1,domain) << " "
                 << getBoundingBoxMin(2,domain) << " " << getBoundingBoxMax(0,domain) << " "
                 << getBoundingBoxMax(1,domain) << " " << getBoundingBoxMax(2,domain) << endl;
      povcfgstrm.close();
    }
    barrier();
  }
}

void KDDecomposition::writeMoleculesToFile(string filename, ParticleContainer* moleculeContainer){

  for(int process = 0; process < _numProcs; process++){
    if(_ownRank==process){
      ofstream checkpointfilestream(filename.c_str(), ios::app);
      checkpointfilestream.precision(20);
      Molecule* tempMolecule;
      for(tempMolecule = moleculeContainer->begin(); tempMolecule != moleculeContainer->end(); tempMolecule = moleculeContainer->next()){
        tempMolecule->write(checkpointfilestream);
      }
      checkpointfilestream.close();
    }
    barrier();
  }
}


//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//$ private Methoden, die von exchangeMolecule benvtigt werden $
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

void KDDecomposition::getPartsToSend(KDNode* sourceArea, KDNode* decompTree, ParticleContainer* moleculeContainer, Domain* domain, vector<int>& procIDs, vector<int>& numMolsToSend, vector<list<Molecule*> >& particlesToSend){
  int haloCellIdxMin[3]; // Assuming a global 3D Cell index, haloCellIdxMin[3] gives the position
                         // of the low local domain corner within this global 3D cell index
  int haloCellIdxMax[3]; // same as heloCellIdxMax, only high instead of low Corner
  for(int dim=0; dim<3; dim++){
    haloCellIdxMin[dim] = sourceArea->_lowCorner[dim]-1;
    haloCellIdxMax[dim] = sourceArea->_highCorner[dim]+1;
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
  for(int neighbCount=0;neighbCount<(int)procIDs.size();neighbCount++) {
    // No need to send particles to the own proc
    if(procIDs[neighbCount]==_ownRank) continue;

    regToSendLow[0] = neighbHaloAreas[6*neighbCount+0]*_cellSize[0];
    regToSendLow[1] = neighbHaloAreas[6*neighbCount+1]*_cellSize[1];
    regToSendLow[2] = neighbHaloAreas[6*neighbCount+2]*_cellSize[2];
    regToSendHigh[0] = (neighbHaloAreas[6*neighbCount+3]+1)*_cellSize[0];
    regToSendHigh[1] = (neighbHaloAreas[6*neighbCount+4]+1)*_cellSize[1];
    regToSendHigh[2] = (neighbHaloAreas[6*neighbCount+5]+1)*_cellSize[2];

    // Not only the neighbouring region itself, but also the periodic copies have to be checked
    int shift[3] = {0,0,0};
    for(int dim=0; dim<3; dim++){
      // If the neighbouring region overlaps the left side of the global domain,
      // a copy of the neighb. region which is shifted to the right has to be examined
      if(regToSendLow[dim] < 0.0 && regToSendHigh[dim] <= domain->getGlobalLength(dim)){
        shift[dim] = 1;
      }
      // same as before, but with shift to the left
      else if(regToSendLow[dim] >= 0.0 && regToSendHigh[dim] > domain->getGlobalLength(dim)){
        shift[dim] = -1;
      }
      // The other cases:
      // neither overlap left or right --> no copies necessary
      // overlap on both sides --> The neighbouring area already covers the whole domain
      //                           (in that dimension), so shifted copies could not
      //                           cover more
    }

    particlesToSend[neighbCount].clear();
    for (int iz = 0; iz<=1; iz++){
      if(iz==1 && shift[2]==0) break; // no shift in z-direction
      for (int iy = 0; iy<=1; iy++){
        if(iy==1 && shift[1]==0) break; // no shift in y-direction
        for (int ix = 0; ix<=1; ix++){
          if(ix==1 && shift[0]==0) break; // no shift in x-direction
          shiftRegLow[0] = regToSendLow[0] + ix*shift[0]*domain->getGlobalLength(0);
          shiftRegLow[1] = regToSendLow[1] + iy*shift[1]*domain->getGlobalLength(1);
          shiftRegLow[2] = regToSendLow[2] + iz*shift[2]*domain->getGlobalLength(2);
          shiftRegHigh[0] = regToSendHigh[0] + ix*shift[0]*domain->getGlobalLength(0);
          shiftRegHigh[1] = regToSendHigh[1] + iy*shift[1]*domain->getGlobalLength(1);
          shiftRegHigh[2] = regToSendHigh[2] + iz*shift[2]*domain->getGlobalLength(2);
          moleculeContainer->getRegion(shiftRegLow, shiftRegHigh, particlesToSend[neighbCount]);
        }
      }
    }
    // store number of particles to be sent to the neighbour
    numMolsToSend[neighbCount] = particlesToSend[neighbCount].size();
  }
}

void KDDecomposition::exchangeNumToSend(vector<int>& procsToSendTo, vector<int>& procsToRecvFrom, vector<int>& numMolsToSend, vector<int>& numMolsToRecv){

  numMolsToRecv.resize(procsToRecvFrom.size());
  MPI_Request request[procsToRecvFrom.size()];
  MPI_Status status;
  // initiate recv calls
  for(int neighbCount=0;neighbCount<(int)procsToRecvFrom.size();neighbCount++) {
    if(procsToRecvFrom[neighbCount]==_ownRank) continue; // don't exchange data with the own process
    MPI_Irecv(&numMolsToRecv[neighbCount], 1, MPI_INT, procsToRecvFrom[neighbCount], 0, MPI_COMM_WORLD, &request[neighbCount]);
  }

  // send all numbers
  for(int neighbCount=0;neighbCount<(int)procsToSendTo.size();neighbCount++) {
    if(procsToSendTo[neighbCount]==_ownRank) continue; // don't exchange data with the own process
    MPI_Send(&numMolsToSend[neighbCount], 1, MPI_INT, procsToSendTo[neighbCount], 0, MPI_COMM_WORLD);
  }
  // wait for the completion of all recv calls
  for(int neighbCount=0;neighbCount<(int)procsToRecvFrom.size();neighbCount++) {
    if(procsToRecvFrom[neighbCount]==_ownRank) continue; // don't exchange data with the own process
    MPI_Wait(&request[neighbCount], &status);
  }
}


void KDDecomposition::transferMolData(vector<int>& procsToSendTo, vector<int>& procsToRecvFrom, vector<int>& numMolsToSend, vector<int>& numMolsToRecv, vector<ParticleData*>& particlesSendBufs, vector<ParticleData*>& particlesRecvBufs){

  // create a MPI Datatype which can store that molecule-data that has to be sent
  MPI_Datatype sendPartType;
  ParticleData::setMPIType(sendPartType);

  MPI_Request request[procsToRecvFrom.size()];
  MPI_Status status;

  // initiate recv calls
  for(int neighbCount=0;neighbCount<(int)procsToRecvFrom.size();neighbCount++) {
    if(procsToRecvFrom[neighbCount]==_ownRank) continue; // don't exchange data with the own process
    MPI_Irecv(particlesRecvBufs[neighbCount], numMolsToRecv[neighbCount], sendPartType, procsToRecvFrom[neighbCount], 0, MPI_COMM_WORLD, &request[neighbCount]);
  }

  // send all numbers
  for(int neighbCount=0;neighbCount<(int)procsToSendTo.size();neighbCount++) {
    if(procsToSendTo[neighbCount]==_ownRank) continue; // don't exchange data with the own process
    MPI_Send(particlesSendBufs[neighbCount], numMolsToSend[neighbCount], sendPartType, procsToSendTo[neighbCount], 0, MPI_COMM_WORLD);
  }

  // wait for the completion of all recv calls
  for(int neighbCount=0;neighbCount<(int)procsToRecvFrom.size();neighbCount++) {
    if(procsToRecvFrom[neighbCount]==_ownRank) continue; // don't exchange data with the own process
    MPI_Wait(&request[neighbCount], &status);
  }
}


void KDDecomposition::createLocalCopies(ParticleContainer* moleculeContainer, Domain* domain, const vector<Component>& components){
  if(_decompValid==false){
    cerr << "KDDecomposition::createLocalCopies failed. Decomposition is invalid" << endl;
    exit(1);
  }

  Molecule* molPtr;
  // molecules that have to be copied, get a new position
  double newPosition[3];

  for(unsigned short d=0;d<3;++d) {
    // only if the local area covers the whole domain in one dimension, copies have to
    // be created. Otherwise, the copies are on other procs and are send there
    // by the method exchangeMolecules
    if(not(_ownArea->_coversWholeDomain[d])) continue;

    molPtr = moleculeContainer->begin();
    while(molPtr!=moleculeContainer->end()){

      const double& rd=molPtr->r(d);
      int copy = 0; // -1: copy to left, 1: copy to right, 0: don't copy
      if(rd < moleculeContainer->get_halo_L(d)) copy = 1;
      if(rd >= domain->getGlobalLength(d)-moleculeContainer->get_halo_L(d)) copy = -1;
      if(copy != 0) {
        // determine the position for the copy of the molecule
        for(unsigned short d2=0; d2<3; d2++){
          // when moving parallel to the coordinate d2 to another process, the
          // local coordinates in d2 change
          if(d2==d) newPosition[d2] = rd+copy*domain->getGlobalLength(d2);
          else newPosition[d2] = molPtr->r(d2);
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



void KDDecomposition::recInitialDecomp(KDNode* fatherNode, KDNode*& ownArea){
  bool coversAll[KDDIM];
  int cellsPerDim[KDDIM];
  for(int dim=0; dim<KDDIM; dim++){
    coversAll[dim] = fatherNode->_coversWholeDomain[dim];
    cellsPerDim[dim] = fatherNode->_highCorner[dim]-fatherNode->_lowCorner[dim]+1;
  }
  int divDir = 0;
  int maxCells = cellsPerDim[0];
  if(cellsPerDim[1] > maxCells) {
    divDir = 1;
    maxCells = cellsPerDim[1];
  }
  if(cellsPerDim[2] > maxCells) {
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

  for(int dim=0; dim<KDDIM; dim++){
    low1[dim] = fatherNode->_lowCorner[dim];
    low2[dim] = fatherNode->_lowCorner[dim];
    high1[dim] = fatherNode->_highCorner[dim];
    high2[dim] = fatherNode->_highCorner[dim];
  }

  low1[divDir] = fatherNode->_lowCorner[divDir];
  high1[divDir] = (fatherNode->_highCorner[divDir]+fatherNode->_lowCorner[divDir])/2;
  low2[divDir] = high1[divDir] + 1;
  high2[divDir] = fatherNode->_highCorner[divDir];

  numProcs1 = fatherNode->_numProcs/2;
  numProcs2 = fatherNode->_numProcs - numProcs1;
  id1 = fatherNode->_nodeID + 1;
  id2 = fatherNode->_nodeID + 2*numProcs1;
  owner1 = fatherNode->_owningProc;
  owner2 = owner1+numProcs1;
  fatherNode->_child1 = new KDNode(numProcs1, low1, high1, id1, owner1, coversAll);
  fatherNode->_child2 = new KDNode(numProcs2, low2, high2, id2, owner2, coversAll);

  if(numProcs1 > 1) recInitialDecomp(fatherNode->_child1, ownArea);
  else if (owner1==_ownRank) ownArea = fatherNode->_child1;
  if(numProcs2 > 1) recInitialDecomp(fatherNode->_child2, ownArea);
  else if (owner2==_ownRank) ownArea = fatherNode->_child2;
}

void KDDecomposition::recDecomp(KDNode* fatherNode, KDNode*& ownArea){

  if(fatherNode->_numProcs == 1) {
    if(fatherNode->_owningProc == _ownRank) ownArea = fatherNode;
    else cerr << "ERROR in recDecomp: called with a leaf node" << endl;
    return;
  }
  bool coversAll[KDDIM];
  int cellsPerDim[KDDIM];
  for(int dim=0; dim<KDDIM; dim++){
    coversAll[dim] = fatherNode->_coversWholeDomain[dim];
    cellsPerDim[dim] = fatherNode->_highCorner[dim]-fatherNode->_lowCorner[dim]+1;
  }

  // calculate particles per layer
  vector<int> particlesPerLayer[KDDIM];
  for(int dim=0; dim<KDDIM; dim++){
    particlesPerLayer[dim].resize(fatherNode->_highCorner[dim]-fatherNode->_lowCorner[dim]+1);
    for(int i=0; i<=fatherNode->_highCorner[dim]-fatherNode->_lowCorner[dim]; i++){
      particlesPerLayer[dim][i] = 0;
    }
  }
  for(int iz=0; iz<=fatherNode->_highCorner[2]-fatherNode->_lowCorner[2]; iz++) {
    for(int iy=0; iy<=fatherNode->_highCorner[1]-fatherNode->_lowCorner[1]; iy++) {
      for(int ix=0; ix<=fatherNode->_highCorner[0]-fatherNode->_lowCorner[0]; ix++) {
        int partsInCell = _numParticlesPerCell[_globalCellsPerDim[0]*(iz*_globalCellsPerDim[1] + iy) + ix];
        particlesPerLayer[0][ix] += partsInCell;
        particlesPerLayer[1][iy] += partsInCell;
        particlesPerLayer[2][iz] += partsInCell;
      }
    }
  }


  vector<double> partCost;
  double partFactor = 1.0;
  vector<double> sepCost;
  double sepFactor = 1.0;


  int divDir = 0;
  int divIdx = 0;
  double maxProcCost = DBL_MAX;
  int numProcsLeft = 0;
  int numProcsRight = 0;


  for(int dim=0; dim<3; dim++){
    partCost.resize(fatherNode->_highCorner[dim]-fatherNode->_lowCorner[dim]+1);
    sepCost.resize(fatherNode->_highCorner[dim]-fatherNode->_lowCorner[dim]+1);
    for(int i=0; i<=fatherNode->_highCorner[dim]-fatherNode->_lowCorner[dim]; i++){
      if(i==0) {
        partCost[i] = (double) particlesPerLayer[dim][i];
        sepCost[i] = sepFactor* (double) (particlesPerLayer[dim][i] + particlesPerLayer[dim][i+1]);
      }
      else if(i==fatherNode->_highCorner[dim]) {
        partCost[i] = partFactor*(partCost[i-1] + (double) particlesPerLayer[dim][i]);
        sepCost[i] = 0.0;
      }
      else{
        partCost[i] = partFactor*(partCost[i-1] + (double) particlesPerLayer[dim][i]);
        sepCost[i] = sepFactor* (double) (particlesPerLayer[dim][i] + particlesPerLayer[dim][i+1]);
      }


    }
    for(int i=0; i<fatherNode->_highCorner[dim]-fatherNode->_lowCorner[dim]; i++){
      double partCostLeft = partCost[i];
      double partCostRight = partCost[fatherNode->_highCorner[dim]-fatherNode->_lowCorner[dim]] - partCostLeft;

      double costLeft = partCostLeft + sepCost[i];
      double costRight = partCostRight + sepCost[i];
      double costTotal = costLeft + costRight;
      double optCostPerProc = costTotal/((double) fatherNode->_numProcs);
      int numProcsLeft1 = (int) floor(costLeft/optCostPerProc);
      int numProcsLeft2 = (int) ceil(costLeft/optCostPerProc);
      int numProcsRight1 = fatherNode->_numProcs - numProcsLeft1;
      int numProcsRight2 = fatherNode->_numProcs - numProcsLeft2;
      double maxProcCost1 = max(costLeft/(double)numProcsLeft1, costRight/(double)numProcsRight1);
      double maxProcCost2 = max(costLeft/(double)numProcsLeft2, costRight/(double)numProcsRight2);

      if(maxProcCost1 < maxProcCost){
        divDir = dim;
        divIdx = i+fatherNode->_lowCorner[dim];
        maxProcCost = maxProcCost1;
        numProcsLeft = numProcsLeft1;
        numProcsRight = fatherNode->_numProcs - numProcsLeft;
      }
      if(maxProcCost2 < maxProcCost){
        divDir = dim;
        divIdx = i+fatherNode->_lowCorner[dim];
        maxProcCost = maxProcCost2;
        numProcsLeft = numProcsLeft2;
        numProcsRight = fatherNode->_numProcs - numProcsLeft;
      }

      if(numProcsLeft <= 0 || numProcsLeft >= fatherNode->_numProcs){
        cerr << "ERROR in recDecomp, part of the domain was not assigned to a proc" << endl;
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

  for(int dim=0; dim<KDDIM; dim++){
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
  id2 = fatherNode->_nodeID + 2*numProcsLeft;
  owner1 = fatherNode->_owningProc;
  owner2 = owner1+numProcsLeft;
  fatherNode->_child1 = new KDNode(numProcsLeft, low1, high1, id1, owner1, coversAll);
  fatherNode->_child2 = new KDNode(numProcsRight, low2, high2, id2, owner2, coversAll);


  if(numProcsLeft > 1) recDecomp(fatherNode->_child1, ownArea);
  else if (owner1==_ownRank) ownArea = fatherNode->_child1;
  if(numProcsRight > 1) recDecomp(fatherNode->_child2, ownArea);
  else if (owner2==_ownRank) ownArea = fatherNode->_child2;

}

void KDDecomposition::printDecompTree(KDNode* root, string prefix){
  if(root->_numProcs == 1) {
    cout << prefix << "LEAF: " << root->_nodeID << ", Owner: " << root->_owningProc << ", Corners: ("
         << root->_lowCorner[0] << ", " << root->_lowCorner[1] << ", "<< root->_lowCorner[2] << ") / ("
         << root->_highCorner[0] << ", " << root->_highCorner[1] << ", "<< root->_highCorner[2] << ")" << endl;
  }
  else {
    cout << prefix << "INNER: " << root->_nodeID << ", Owner: " << root->_owningProc << "(" << root->_numProcs << " procs)" << ", Corners: ("
         << root->_lowCorner[0] << ", " << root->_lowCorner[1] << ", "<< root->_lowCorner[2] << ") / ("
         << root->_highCorner[0] << ", " << root->_highCorner[1] << ", "<< root->_highCorner[2] << ")" << endl;
    stringstream childprefix;
    childprefix << prefix << "  ";
    printDecompTree(root->_child1, childprefix.str());
    printDecompTree(root->_child2, childprefix.str());
  }
}



int KDDecomposition::mod(int number, int modulo){
  int result = number % modulo;
  if(result < 0) result+=modulo;
  return result;
}

//##########################################################################
//##########################################################################
//###                                                                    ###
//###                  evtl unnvtige private Methoden                    ###
//###                                                                    ###
//##########################################################################
//##########################################################################

unsigned long KDDecomposition::predictBoundariesPerProc(int numcells, int numprocs){
  double cellsPerProc = numcells/numprocs;
  int cellsPerDim = (int)ceil(pow(cellsPerProc,1.0/(double) KDDIM));
  return 4*cellsPerDim;
}




void KDDecomposition::getOwningProcs(int low[KDDIM], int high[KDDIM], KDNode* decompTree, KDNode* testNode, vector<int>* procIDs, vector<int>* neighbHaloAreas){
  // For areas overlapping the domain given by decompTree, the overlapping part is
  // mapped to the corresponding area on the other side of the domain (periodic boundary)
  // The boolean variable overlap stores for each coordinate direction whether the area overlaps.
  bool overlap[KDDIM];
  for(int dim=0; dim<KDDIM; dim++){
    low[dim] = mod(low[dim],(decompTree->_highCorner[dim]-decompTree->_lowCorner[dim]+1));
    high[dim] = mod(high[dim],(decompTree->_highCorner[dim]-decompTree->_lowCorner[dim]+1));
    if (low[dim] > high[dim]) overlap[dim] = true;
    else overlap[dim] = false;
  }
  // First find out wheter the area intersects the area given by testNode
  bool areasIntersect = true;
  for(int dim=0; dim<KDDIM; dim++){
    if(overlap[dim]){
      if(high[dim]<testNode->_lowCorner[dim] && low[dim] > testNode->_highCorner[dim]) areasIntersect = false;
    }
    else {
      if(high[dim]<testNode->_lowCorner[dim] || low[dim] > testNode->_highCorner[dim]) areasIntersect = false;
    }
  }
  if (areasIntersect){
    // while testNode is still a inner node (more than one proc), recursively call
    // this method for the children
    if (testNode->_numProcs > 1){
      getOwningProcs(low, high, decompTree, testNode->_child1, procIDs, neighbHaloAreas);
      getOwningProcs(low, high, decompTree, testNode->_child2, procIDs, neighbHaloAreas);
    }
    else { // leaf node found, add it (and it's area)
      procIDs->push_back(testNode->_owningProc);
      neighbHaloAreas->push_back(testNode->_lowCorner[0]-1);
      neighbHaloAreas->push_back(testNode->_lowCorner[1]-1);
      neighbHaloAreas->push_back(testNode->_lowCorner[2]-1);
      neighbHaloAreas->push_back(testNode->_highCorner[0]+1);
      neighbHaloAreas->push_back(testNode->_highCorner[1]+1);
      neighbHaloAreas->push_back(testNode->_highCorner[2]+1);
    }
  }
}

void KDDecomposition::getNumParticles(ParticleContainer* moleculeContainer){

  if(_decompValid==false){
    cerr << "KDDecomposition::getNumParticles failed. Decomposition is invalid" << endl;
    exit(1);
  }
  int count = 0;
  double bBMin[3]; // haloBoundingBoxMin
  double bBMax[3]; // haloBoundingBoxMax
  for(int dim=0; dim<3; dim++){
    bBMin[dim] = moleculeContainer->getBoundingBoxMin(dim);// - moleculeContainer->get_halo_L(dim);
    bBMax[dim] = moleculeContainer->getBoundingBoxMax(dim);// + moleculeContainer->get_halo_L(dim);
  }
  Molecule* molPtr = moleculeContainer->begin();
  while(molPtr!=moleculeContainer->end()){
    int cellIndex[3]; // 3D Cell index (local)
    int globalCellIdx[3]; // 3D Cell index (global)

    for(int dim=0; dim<3; dim++){
      cellIndex[dim] = (int) floor((molPtr->r(dim)-bBMin[dim])/_cellSize[dim]);
      globalCellIdx[dim] = _ownArea->_lowCorner[dim] + cellIndex[dim];
      if(globalCellIdx[dim] < 0) globalCellIdx[dim] += _globalCellsPerDim[dim];
      if(globalCellIdx[dim] >=_globalCellsPerDim[dim]) globalCellIdx[dim] -= _globalCellsPerDim[dim];
    }


    _numParticlesPerCell[_globalCellsPerDim[0]*(globalCellIdx[2]*_globalCellsPerDim[1] + globalCellIdx[1]) + globalCellIdx[0]]++;
    molPtr = moleculeContainer->next();
    count++;
  }
  // TODO
  // memory problem:
  // Some MPI-Implementations demand that in a Allreduce command, send and recv buffers
  // must be different. The consequence is, that a temporary array has to be created for
  // the values to be recieved, which that has to be copied back to the original array
  int* numParticlesPerCellTemp = new int[_globalNumCells];
  MPI_Allreduce(_numParticlesPerCell, numParticlesPerCellTemp, _globalNumCells, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  delete [] _numParticlesPerCell;
  _numParticlesPerCell = numParticlesPerCellTemp;

}

unsigned KDDecomposition::Ndistribution(unsigned localN, float* minrnd, float* maxrnd)
{
   int num_procs;
   MPI_Comm_size(this->_collComm.getTopology(), &num_procs);
   unsigned* moldistribution = new unsigned[num_procs];
   MPI_Allgather(&localN, 1, MPI_UNSIGNED, moldistribution, 1, MPI_UNSIGNED, this->_collComm.getTopology());
   unsigned globalN = 0;
   for(int r=0; r < this->_ownRank; r++) globalN += moldistribution[r];
   unsigned localNbottom = globalN;
   globalN += moldistribution[this->_ownRank];
   unsigned localNtop = globalN;
   for(int r = this->_ownRank + 1; r < num_procs; r++) globalN += moldistribution[r];
   delete []moldistribution;
   *minrnd = (float)localNbottom / globalN;
   *maxrnd = (float)localNtop / globalN;
   return globalN;
}

void KDDecomposition::assertIntIdentity(int IX)
{
   if(this->_ownRank)
   {
      MPI_Send(&IX, 1, MPI_INT, 0, 2*_ownRank+17, this->_collComm.getTopology());
   }
   else
   {
      int recv;
      int num_procs;
      MPI_Comm_size(this->_collComm.getTopology(), &num_procs);
      MPI_Status s;
      for(int i=1; i<num_procs; i++)
      {
         MPI_Recv(&recv, 1, MPI_INT, i, 2*i+17, this->_collComm.getTopology(), &s);
         if(recv != IX)
         {
            cout << "SEVERE ERROR: IX is " << IX << " for rank 0, but " << recv << " for rank " << i << ".\n";
            MPI_Finalize();
            exit(911);
         }
      }
      cout << "IX = " << recv << " for all " << num_procs << " ranks.\n";
   }
}

void KDDecomposition::assertDisjunctivity(TMoleculeContainer* mm)
{
   Molecule* m;
   if(this->_ownRank)
   {
      int tid;
      for(m = mm->begin(); m != mm->end(); m = mm->next())
      {
         tid = m->id();
         MPI_Send(&tid, 1, MPI_INT, 0, 2674+_ownRank, this->_collComm.getTopology());
      }
      tid = -1;
      MPI_Send(&tid, 1, MPI_INT, 0, 2674+_ownRank, this->_collComm.getTopology());
   }
   else
   {
      int recv;
      map<int, int> check;
      for(m = mm->begin(); m != mm->end(); m = mm->next())
      {
         check[m->id()] = 0;
      }

      int num_procs;
      MPI_Comm_size(this->_collComm.getTopology(), &num_procs);
      MPI_Status s;
      for(int i=1; i < num_procs; i++)
      {
         bool cc = true;
         while(cc)
         {
            MPI_Recv(&recv, 1, MPI_INT, i, 2674+i, this->_collComm.getTopology(), &s);
            if(recv == -1) cc = false;
            else
            {
               if(check.find(recv) != check.end())
               {
                  cout << "\nSEVERE ERROR. Ranks " << check[recv] << " and "
                       << i << " both propagate ID " << recv << ". Aborting.\n";
                  MPI_Finalize();
                  exit(2674);
               }
               else check[recv] = i;
            }
         }
      }
      cout << "\nData consistency checked. No duplicate IDs detected among " << check.size()
           << " entries.\n";
   }
}


//void KDDecomposition::balance(ParticleContainer* moleculeContainer, const vector<Component>& components, Domain* domain){

//  balanceAndExchange(false, moleculeContainer, components, domain);

//moleculeContainer->rebuild(bBoxMin, bBoxMax);



// Number of force calculations for each cell
//vector<vector<vector<int> > > numForceCalcs;
// Number of particles for each cell (including halo)
//vector<vector<vector<int> > > numParticles;
// If e.g. the domain is divided by a layer perpendicular
// to the x-axis (dim 0), with 32 layers of cells on the "left" side
// and size[0]-32 layers of cells on the "right" side,
// separatedPairs[0][31] (!!! check correct index!!!) is the number
// of pairs whose particles are on two different sides
//vector<int> separatedPairs[KDDIM];
// Using the same separator as for separatedPairs, compTime[0][32]
// describes the time needed for the calculation of the force
// pairs (without parallelisation overhead, only pure calculation)
// compTime[*][size[*]-1] is therefore the calculation time for the
// whole domain
//vector<double> compTime[KDDIM];

//getNumForceCalcs(moleculeContainer);




// set costs for separatedPairs and compTime to zero
//for(int dim=0; dim<KDDIM; dim++){
//  for(int i=0; i<globalSize[dim]; i++) {
//    separatedPairs[dim][i] = 0;
//compTime[dim][i] = 0.0;
//  }
//}
// costs for separatedPairs + added comp costs
//for(int iz=1; iz<=oldSize[2]; iz++) {
//  for(int iy=1; iy<=oldSize[1]; iy++) {
//    for(int ix=1; ix<=oldSize[0]; ix++) {
//      separatedPairs[0][ix-1+_ownArea->_lowCorner[0]] += numParticles[ix][iy][iz]*numParticles[ix+1][iy][iz];
//      separatedPairs[1][iy-1+_ownArea->_lowCorner[1]] += numParticles[ix][iy][iz]*numParticles[ix][iy+1][iz];
//      separatedPairs[2][iz-1+_ownArea->_lowCorner[2]] += numParticles[ix][iy][iz]*numParticles[ix][iy][iz+1];
//      //compTime[0][ix-1+_ownArea->_lowCorner[0]] += _timePerForceCalc*(double)numForceCalcs[ix-1][iy-1][iz-1];
//      //compTime[1][iy-1+_ownArea->_lowCorner[1]] += _timePerForceCalc*(double)numForceCalcs[ix-1][iy-1][iz-1];
//      //compTime[2][iz-1+_ownArea->_lowCorner[2]] += _timePerForceCalc*(double)numForceCalcs[ix-1][iy-1][iz-1];
//    }
//  }
//}

// kosten aufaddieren -> compTime[i] enthdlt
// die Kosten aller Zellen links von numForceCalcs[i][*]
//for(int dim=0; dim<KDDIM; dim++){
//  for(int i=1; i<globalSize[dim]; i++) {
//compTime[dim][i] += compTime[dim][i-1];
//  }
//}

// Zeit pro Randzelle, setzt sich zusammen aus:
// - durchschnittliche Paare, die den Rand |berschreiten
// - * (Rechenzeit pro Paar (timePerForceCalc)
// -    + Kommunikationszeit pro Paar)
//_timePerBoundary = 5.0 * (_timePerForceCalc + 0.5);

// alle verschiedenen Aufteilungen ausprobieren:
// x-Richtung
//for(int dim=0; dim<1; dim++){
//  for(int i=1; i<_numProcs; i++){
//    int cellsPerLayer = 1;
//    for(int d2=0; d2<KDDIM; d2++){
//if (d2!=dim) cellsPerLayer *= size[dim]; //???
//    }
//       double totalload = compTime[dim][size[dim]-1];
//       int leftmaxindex = 0;
//       // Rechenzeit einer Seite ergibt sich wie folgt:
//       // - Rechenzeit der Seite (|ber xCompCosts) durch #procs teilen
//       // - Volumen durch #procs teilen
//       // - Anzahl Randseiten pro proc abschaetzen
//       // - Randseiten multiplizieren mit Zeit pro Randzeit
//       // - SeparatedPairs durch #procs teilen
//       // - Rechenzeit+Randzeit+xSeparatedPairs ergibt Gesamtzeit
//       double calctimeleft = (compTime[dim][leftmaxindex]/i +
//         _timePerBoundary*predictBoundariesPerProc(cellsPerLayer*(leftmaxindex+1),i) +
//            _timePerForceCalc*separatedPairs[dim][leftmaxindex+1]/i);
//       double calctimeright= ((totalload - compTime[dim][leftmaxindex])/(_numProcs-i) +
//            _timePerBoundary*predictBoundariesPerProc(cellsPerLayer*(size[dim]-(leftmaxindex+1)),_numProcs - i) + _timePerForceCalc*separatedPairs[dim][leftmaxindex+1]/(_numProcs-i));
//       double currenttime;// = max(calctimeleft, calctimeright);
//       double besttime = max(calctimeleft, calctimeright);
//       int bestindex = 0;


//       while((leftmaxindex < size[dim])){
//   calctimeleft = (compTime[dim][leftmaxindex]/i +
//       _timePerBoundary*predictBoundariesPerProc(cellsPerLayer*(leftmaxindex+1),i) +
//       _timePerForceCalc*separatedPairs[dim][leftmaxindex+1]/i);
//   calctimeright= ((totalload - compTime[dim][leftmaxindex])/(_numProcs-i) +
//       _timePerBoundary*predictBoundariesPerProc(cellsPerLayer*(size[dim]-(leftmaxindex+1)),_numProcs - i) + _timePerForceCalc*separatedPairs[dim][leftmaxindex+1]/(_numProcs-i));


//   currenttime = max(calctimeleft, calctimeright);
//   if(currenttime < besttime) {
//     bestindex = leftmaxindex;
//     besttime = currenttime;
//   }
//   //if(calctimeright > totaltime)  {cout << "BLA";}// loopcount = 111;
//   leftmaxindex++;

//       }
//       // !!! tut nicht immer (underflow)
//       if(_ownRank==0) cout << "i / Index / totaltime: "<< i << " / "<< bestindex << " / " << besttime << endl;
//  }
//}

//  barrier();
//   if(_ownRank==0) cout << "total costs x: " << compTime[0][size[0]-1] << endl;

//}






// void KDDecomposition::setCellOwnership(KDNode* localArea, KDNode* decompTree){
//   int size[KDDIM];
//   for(int dim=0; dim<KDDIM; dim++){
//     size[dim] = localArea->_highCorner[dim] - localArea->_lowCorner[dim] + 1;
//   }
//   int cellPos[KDDIM];
//   // untere Zeile (incl. Ecke)
//   cellPos[1] = localArea->_lowCorner[1]-1;
//   for(int ix = 0; ix < size[0]+2; ix++){
//     cellPos[0] = localArea->_lowCorner[0]+ix-1;
//     _localCells[0*(size[0]+2)+ix]._owningProc = getOwningProcPerBound(cellPos, decompTree );
//   }
//   // obere Zeile (incl. Ecke)
//   cellPos[1] = localArea->_highCorner[1]+1;
//   for(int ix = 0; ix < size[0]+2; ix++){
//     cellPos[0] = localArea->_lowCorner[0]+ix-1;
//     _localCells[(size[1]+1)*(size[0]+2)+ix]._owningProc = getOwningProcPerBound(cellPos, decompTree );
//   }
//   // linke Zeile (ohne Ecke)
//   cellPos[0] = localArea->_lowCorner[0]-1;
//   for(int iy = 1; iy < size[1]+1; iy++){
//     cellPos[1] = localArea->_lowCorner[1]+iy-1;
//     _localCells[iy*(size[0]+2)+0]._owningProc = getOwningProcPerBound(cellPos, decompTree );
//   }
//    // rechte Zeile (ohne Ecke)
//    cellPos[0] = localArea->_highCorner[0]+1;
//    for(int iy = 1; iy < size[1]+1; iy++){
//      cellPos[1] = localArea->_lowCorner[1]+iy-1;
//      _localCells[iy*(size[0]+2)+(size[0]+1)]._owningProc = getOwningProcPerBound(cellPos, decompTree );
//    }
// }

// int KDDecomposition::getOwningProc(int cellPos[KDDIM], KDNode* decompTree){
//   // Testen, ob die Zelle im Gebiet liegt
//   if(decompTree->cellBelongsToRegion(cellPos)){
//     // Falls ja, und es ein Blatt ist, dann den Besitzern zur|ckgeben
//     if (decompTree->_numProcs == 1) return decompTree->_owningProc;
//     // ansonsten bei den beiden Kindern weitersuchen
//     else {
//       int id1 = getOwningProc(cellPos, decompTree->_child1);
//       int id2 = getOwningProc(cellPos, decompTree->_child2);
//       if (id1 < 0 && id2 >= 0) return id2;
//       else if (id1 >= 0 && id2 < 0) return id1;
//       else {
//   cout << "ERROR, CELL IS EITHER OWNED BY NO, OR BY MANY PROCESSES" << endl;
//   exit(27);
//       }
//     }
//   }
//   // Wenn die Zelle nicht im Gebiet liegt, -1 zur|ckgeben
//   else return -1;
// }









// void KDDecomposition::readInData(){
//   ifstream datFS;
//   datFS.open("4balls.dat");
//   int globalSize[KDDIM];
//   int localSize[KDDIM];
//   int dummy;
//   datFS >> globalSize[0] >> globalSize[1];
//   if(_ownArea == NULL){
//     cout << "ERROR: Proc " << _ownRank << " has no own area" << endl;
//   }
//   for(int dim=0; dim<KDDIM; dim++){
//     localSize[dim] = _ownArea->_highCorner[dim]-_ownArea->_lowCorner[dim] + 1;
//   }
//   if(_ownRank==0) cout << "localSize: " << localSize[0] << " / " << localSize[1] << endl;
//   _numParticlesPerCell = new int[(localSize[0]+2)*(localSize[1]+2)];
//   _localCells.resize((localSize[0]+2)*(localSize[1]+2));

//   if(_ownRank==0){
//     cout << "LOW CORNER: "<< _ownArea->_lowCorner[0] << " / " << _ownArea->_lowCorner[1] << endl;
//     cout << "HIGH CORNER: "<< _ownArea->_highCorner[0] << " / " << _ownArea->_highCorner[1] << endl;
//   }

//   // read in _numParticlesPerCell
//   for(int iy=globalSize[1]-1; iy>=0; iy--){
//     for(int ix=0; ix<globalSize[0]; ix++){
//       int cellPos[2];
//       // evtl. Verschiebung des Referenzrahmens muss kompensiert werden
//       // Angenommen man legt ein Zellgitter |ber das "reale" Gebiet, so dass
//       // die unterste Zelle genau im untersten Eck des Gebiets liegt. Nun
//       // deckt das durch _decompTree abgedeckte Gebiet zwar genau das gleiche
//       // Gebiet ab, aber die unterste Zelle ist gegen|ber dem "realen" Zellgitter
//       // verschoben um _initOffset (mvglich durch periodischen Rand). Um nun zu pr|fen, ob eine
//       // Zelle in dem Gebiet _ownArea liegt, muss also zundchst die Koordinate
//       // angepasst werden. Um die Position bez|glich des neuen Gitters zu erhalten,
//       // muss von der urspr|nglichen Koordinate _initOffset abgezogen werden.
//       cellPos[0] = ix - _initOffset[0];
//       cellPos[1] = iy - _initOffset[1];
//       // Das Gebiet des lokalen Prozesses ist nur ein Teilgebiet des neuen Referenzrahmens.
//       // Analog zu oben muss f|r das Berechnen der "lokalen" Position das Gebiet erneut
//       // verschoben werden
//       int lix = cellPos[0] - _ownArea->_lowCorner[0]+1;
//       int liy = cellPos[1] - _ownArea->_lowCorner[1]+1;

//       if(_ownArea->cellBelongsToRegion(cellPos)){
//   datFS >> _numParticlesPerCell[liy*(localSize[0]+2)+lix];
//   _localCells[liy*(localSize[0]+2)+lix]._numPart = _numParticlesPerCell[liy*(localSize[0]+2)+lix];
//   if(_ownArea->regionCellBelongsToBoundary(cellPos)){
//     _localCells[liy*(localSize[0]+2)+lix]._status = 1;
//     _localCells[liy*(localSize[0]+2)+lix]._owningProc = _ownRank;
//   }
//   else{
//     _localCells[liy*(localSize[0]+2)+lix]._status = 2;
//     _localCells[liy*(localSize[0]+2)+lix]._owningProc = _ownRank;
//   }
//       }
//       else{
//   datFS >> dummy;
//       }
//     }
//   }
//   datFS.close();

// }


// void KDDecomposition::printArea(KDNode* area){
//   int areaSize[KDDIM];
//   for(int dim=0; dim<KDDIM; dim++){
//     areaSize[dim] = area->_highCorner[dim] - area->_lowCorner[dim] + 1;
//   }
//   for(int iy = 0; iy < areaSize[1]+2; iy++){
//     for(int ix = 0; ix < areaSize[0]+2; ix++){
//       //cout << _localCells[iy*(areaSize[0]+2)+ix]._status << " ";
//       cout << _localCells[iy*(areaSize[0]+2)+ix]._owningProc << " ";
//     }
//     cout << endl;
//   }
// }


// int KDDecomposition::getOwningProcPerBound(int cellPos[KDDIM], KDNode* decompTree){
//   for(int dim=0; dim<KDDIM; dim++){
//     if (cellPos[dim] < 0) cellPos[dim]+=(decompTree->_highCorner[dim]+1);
//     if (cellPos[dim] > decompTree->_highCorner[dim]) cellPos[dim]-=(decompTree->_highCorner[dim]+1);
//   }
//   return getOwningProc(cellPos, decompTree);
// }


