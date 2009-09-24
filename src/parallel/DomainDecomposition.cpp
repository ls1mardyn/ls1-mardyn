#include "parallel/DomainDecomposition.h"
#include "molecules/Molecule.h"
#include "datastructures/ParticleContainer.h"
#include "Domain.h"
#include "parallel/ParticleData.h"
#include "Logger.h"
using Log::global_log;

using namespace std;


DomainDecomposition::DomainDecomposition(int *argc, char ***argv){

  int period[DIM];       // 1(true) when using periodic boundary conditions in the corresponding dimension
  int reorder;           // 1(true) if the ranking may be reordered by MPI_Cart_create
  int i;                 // loop counter
  int num_procs;         // Number of processes
  int procnamelen;

  // All boundary conditions are periodic
  for(i=0; i<DIM; i++){
    period[i] = 1;
  }

  MPI_Comm world = MPI_COMM_WORLD;

  // Allow reordering of process ranks
  reorder = 1;

  //MPI_Init(argc,argv);
  MPI_Get_processor_name(_processorName, &procnamelen);
  MPI_Comm_size(world,&num_procs); // determine total number of procs
  setGridSize(num_procs); // determine number of procs in each dimension
  MPI_Cart_create(world,DIM,_gridSize,period,reorder,&_commTopology); // create torus topology

  MPI_Comm_rank(_commTopology, &_ownRank);
  MPI_Cart_coords(_commTopology, _ownRank, DIM, _coords);
  _neighbours[XLOWER]  = getRank(_coords[0]-1,_coords[1],_coords[2]);
  _neighbours[XHIGHER] = getRank(_coords[0]+1,_coords[1],_coords[2]);
  _neighbours[YLOWER]  = getRank(_coords[0],_coords[1]-1,_coords[2]);
  _neighbours[YHIGHER] = getRank(_coords[0],_coords[1]+1,_coords[2]);
  _neighbours[ZLOWER]  = getRank(_coords[0],_coords[1],_coords[2]-1);
  _neighbours[ZHIGHER] = getRank(_coords[0],_coords[1],_coords[2]+1);
}

DomainDecomposition::~DomainDecomposition(){
  MPI_Finalize();
}


void DomainDecomposition::exchangeMolecules(ParticleContainer* moleculeContainer, const vector<Component>& components, Domain* domain, double rc){

  double rmin[3]; // lower corner of the process-specific domain //PARALLEL
  double rmax[3];
  double halo_L[3]; // width of the halo strip //PARALLEL
  for(int i=0; i<3; i++){
    rmin[i] =  moleculeContainer->getBoundingBoxMin(i);
    rmax[i] =  moleculeContainer->getBoundingBoxMax(i);
    halo_L[i] =  moleculeContainer->get_halo_L(i);
  }
  vector<int> compCount;
  compCount.resize(1);
  countMolecules(moleculeContainer, compCount);
  //cout << "Rank " << domain->getlocalRank() << "Num molecules: " << compCount[0] << endl;

  vector<ParticleData*> particlesSendBufs;
  vector<ParticleData*> particlesRecvBufs;
  vector<int> numPartsToSend;
  vector<int> numPartsToRecv;
  particlesSendBufs.resize(6); // left + right neighbour
  particlesRecvBufs.resize(6); // left + right neighbour
  numPartsToSend.resize(6);
  numPartsToRecv.resize(6);

  MPI_Status status;
  MPI_Status  send_statuses[6];
  MPI_Status  recv_statuses[6];
  MPI_Request send_requests[6];
  MPI_Request recv_requests[6];

  // create a MPI Datatype which can store that molecule-data that has to be sent
  MPI_Datatype sendPartType;
  ParticleData::setMPIType(sendPartType);

  int direction; // direction (0=low/1=high) of molecule movement

  for(unsigned short d=0;d<3;++d){
    // when moving a particle across a periodic boundary, the molecule position has to change
    // these offset specify for each dimension (x, y and z) and each direction ("left"/lower
    // neighbour and "right"/higher neighbour, how the paritcle coordinates have to be changed.
    // e.g. for dimension x (d=0) and a process on the left boundary of the domain, particles
    // moving to the left get the length of the whole domain added to their x-value
    double offsetLower[3];
    double offsetHigher[3];
    offsetLower[d]  = 0.0;
    offsetHigher[d] = 0.0;

    // process on the left boundary
    if(_coords[d] == 0)
      offsetLower[d] = domain->getGlobalLength(d);
    // process on the right boundary
    if(_coords[d] == _gridSize[d]-1)
      offsetHigher[d] = -domain->getGlobalLength(d);

    double regToSendLow[3];  // Region that belongs to a neighbouring process
    double regToSendHigh[3]; // -> regToSendLow
    for(direction=0; direction<2; direction++){
      // find the region that each neighbour will get
      for(int i=0; i<3; i++){
        regToSendLow[i] = rmin[i]-halo_L[i];
        regToSendHigh[i] = rmax[i]+halo_L[i];
      }
      if(direction==0) regToSendHigh[d] = rmin[d] + halo_L[d];
      else             regToSendLow[d]  = rmax[d] - halo_L[d];

      list<Molecule*> particlePtrsToSend;
      moleculeContainer->getRegion(regToSendLow, regToSendHigh, particlePtrsToSend);

      // initialize send buffer
      numPartsToSend[2*d+direction] = particlePtrsToSend.size();
      particlesSendBufs[2*d+direction] = new ParticleData[numPartsToSend[2*d+direction]];

      std::list<Molecule*>::iterator particlePtrIter;
      int partCount = 0;
      for(particlePtrIter = particlePtrsToSend.begin(); particlePtrIter!=particlePtrsToSend.end(); particlePtrIter++){
        // copy relevant data from the Molecule to ParticleData type
        ParticleData::setParticleData(particlesSendBufs[2*d+direction][partCount], **particlePtrIter);
        // add offsets for particles transfered over the periodic boundary
        if(d==0 && direction==0) particlesSendBufs[2*d+direction][partCount].rx += offsetLower[0];
        if(d==1 && direction==0) particlesSendBufs[2*d+direction][partCount].ry += offsetLower[1];
        if(d==2 && direction==0) particlesSendBufs[2*d+direction][partCount].rz += offsetLower[2];
        if(d==0 && direction==1) particlesSendBufs[2*d+direction][partCount].rx += offsetHigher[0];
        if(d==1 && direction==1) particlesSendBufs[2*d+direction][partCount].ry += offsetHigher[1];
        if(d==2 && direction==1) particlesSendBufs[2*d+direction][partCount].rz += offsetHigher[2];
        partCount++;
      }
    }

    // Communicate to lower and higher neighbour
    for(direction=0; direction <=1; direction++){
      // Send to lower, receive from upper
      // Send number of values that have to be sent
      unsigned long numsend = numPartsToSend[2*d+direction];
      unsigned long numrecv;
      MPI_Sendrecv(&numsend, 1, MPI_UNSIGNED_LONG, _neighbours[2*d+direction], 99,
                   &numrecv, 1, MPI_UNSIGNED_LONG, _neighbours[2*d+(direction+1)%2], 99, _commTopology, &status);

      // initialize receive buffer
      particlesRecvBufs[2*d+direction] = new ParticleData[numrecv];
      numPartsToRecv[2*d+direction] = numrecv;

      // Send values to lower/upper and receive values from upper/lower
      MPI_Isend(particlesSendBufs[2*d+direction], numsend, sendPartType, _neighbours[2*d+direction], 99, _commTopology, &send_requests[2*d+direction]);
      MPI_Irecv(particlesRecvBufs[2*d+direction], numrecv, sendPartType, _neighbours[2*d+(direction+1)%2], 99, _commTopology, &recv_requests[2*d+direction]);

    }
    for(direction=0; direction <=1; direction++){
      unsigned long numrecv = numPartsToRecv[2*d+direction];
      MPI_Wait(&send_requests[2*d+direction], &send_statuses[2*d+direction]);
      MPI_Wait(&recv_requests[2*d+direction], &recv_statuses[2*d+direction]);
      // insert received molecules into list of molecules
      for(unsigned i=0; i<numrecv; i++){
        ParticleData newMol = particlesRecvBufs[2*d+direction][i];
        Molecule m1 = Molecule(newMol.id, newMol.cid, newMol.rx, newMol.ry, newMol.rz, newMol.vx, newMol.vy, newMol.vz,
                               newMol.qw, newMol.qx, newMol.qy, newMol.qz, newMol.Dx, newMol.Dy, newMol.Dz, &components);
        moleculeContainer->addParticle(m1);
      }
      // free memory
      delete [] particlesRecvBufs[2*d+direction];
      delete [] particlesSendBufs[2*d+direction];
    }
  }
}


void DomainDecomposition::balanceAndExchange(bool balance, ParticleContainer* moleculeContainer, const vector<Component>& components, Domain* domain, double rc){
  exchangeMolecules(moleculeContainer, components, domain, rc);
}


bool DomainDecomposition::procOwnsPos(double x, double y, double z, Domain* domain){
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


double DomainDecomposition::guaranteedDistance(double x, double y, double z, Domain* domain){
  double xdist = 0;
  double ydist = 0;
  double zdist = 0;
  if (x < getBoundingBoxMin(0, domain)){
    xdist = getBoundingBoxMin(0, domain) - x;
  }
  else if(x >= getBoundingBoxMax(0, domain)){
    xdist = x - getBoundingBoxMax(0, domain);
  }
  if (y < getBoundingBoxMin(1, domain)){
    ydist = getBoundingBoxMin(1, domain) - y;
  }
  else if(y >= getBoundingBoxMax(1, domain)){
    ydist = y - getBoundingBoxMax(1, domain);
  }
  if (z < getBoundingBoxMin(2, domain)){
    zdist = getBoundingBoxMin(2, domain) - z;
  }
  else if(z >= getBoundingBoxMax(2, domain)){
    zdist = z - getBoundingBoxMax(2, domain);
  }
  return sqrt(xdist*xdist + ydist*ydist + zdist*zdist);
}


int DomainDecomposition::countMolecules(ParticleContainer* moleculeContainer, vector<int> &compCount){
  vector<int> localCompCount;
  localCompCount.resize(compCount.size());
  for(unsigned int i=0; i<localCompCount.size(); i++){
    localCompCount[i] = 0;
  }
  Molecule* tempMolecule;
  for(tempMolecule = moleculeContainer->begin();
      tempMolecule != moleculeContainer->end();
      tempMolecule = moleculeContainer->next()){
    localCompCount[tempMolecule->componentid()] += 1;
  }
  int numMolecules = 0;
  for(unsigned int i=0; i<localCompCount.size(); i++){
    MPI_Allreduce(&localCompCount[i], &compCount[i], 1, MPI_INT, MPI_SUM, _commTopology);
    numMolecules += compCount[i];
  }
  return numMolecules;
}


double DomainDecomposition::getBoundingBoxMin(int dimension, Domain* domain) {
  return _coords[dimension]*domain->getGlobalLength(dimension)/_gridSize[dimension];
}

double DomainDecomposition::getBoundingBoxMax(int dimension, Domain* domain) {
  return (_coords[dimension]+1)*domain->getGlobalLength(dimension)/_gridSize[dimension];
}


void DomainDecomposition::printDecomp(string filename, Domain* domain){
  int numprocs;
  MPI_Comm_size(_commTopology, &numprocs);

  if(_ownRank==0) {
    ofstream povcfgstrm(filename.c_str());
    povcfgstrm << "size " << domain->getGlobalLength(0) << " " << domain->getGlobalLength(1) << " " << domain->getGlobalLength(2) << endl;
    povcfgstrm << "cells " << _gridSize[0] << " " << _gridSize[1] << " " << _gridSize[2] << endl;
    povcfgstrm << "procs " << numprocs << endl;
    povcfgstrm << "data DomainDecomp" << endl;
    povcfgstrm.close();
  }


  for(int process = 0; process < numprocs; process++){
    if(_ownRank==process){
      ofstream povcfgstrm(filename.c_str(), ios::app);
      povcfgstrm << _coords[2]*_gridSize[0]*_gridSize[1] + _coords[1]*_gridSize[0] + _coords[0] << " " << _ownRank << endl;
      povcfgstrm.close();
    }
    barrier();
  }
}


void DomainDecomposition::writeMoleculesToFile(string filename, ParticleContainer* moleculeContainer){

  int numprocs;
  MPI_Comm_size(_commTopology, &numprocs);

  for(int process = 0; process < numprocs; process++){
    if(_ownRank==process){
      ofstream checkpointfilestream(filename.c_str(), ios::app);
      Molecule* tempMolecule;
      for(tempMolecule = moleculeContainer->begin(); tempMolecule != moleculeContainer->end(); tempMolecule = moleculeContainer->next()){
        tempMolecule->write(checkpointfilestream);
      }
      checkpointfilestream.close();
    }
    barrier();
  }
}


inline int DomainDecomposition::getRank(int x, int y, int z){
  int neigh_coords[DIM]; // Array for the coordinates
  int neigh_rank;        // Rank of the neighbour
  neigh_coords[0] = x;
  neigh_coords[1] = y;
  neigh_coords[2] = z;
  MPI_Cart_rank(_commTopology,neigh_coords,&neigh_rank);
  return(neigh_rank);
}

int DomainDecomposition::getNumProcs(){
  int numProcs;
  MPI_Comm_size(_commTopology, &numProcs);
  return numProcs;
}


const char* DomainDecomposition::getProcessorName() const {
  return _processorName;
}


double DomainDecomposition::getTime(){
  return MPI_Wtime();
}



void DomainDecomposition::setGridSize(int num_procs) {
  int remainder;      // remainder during the calculation of the prime factors
  int i;              // counter
  int num_factors;    // number of prime factors
  int *prime_factors; // array for the prime factors

  // Set the initial number of processes in each dimension to zero
  for(i=0;i<DIM;i++){
    _gridSize[i] = 1;
  }

  remainder = num_procs;

  // The maximal number of prime factors of a number n is log2(n)
  prime_factors = new int[int(log2(num_procs))];

  num_factors = 0;
  // calculate prime numbers
  for(i=2; i<=remainder;i++){
    while(remainder%i == 0){ // -> i is prime factor
      remainder = remainder/i;
      prime_factors[num_factors] = i;
      num_factors++;
    }
  }

  i = num_factors-1;
  while(i>=0){
    if (_gridSize[0] <= _gridSize[1] && _gridSize[0] <= _gridSize[2]){
      _gridSize[0] *= prime_factors[i];
    }
    else if(_gridSize[1] <= _gridSize[0] && _gridSize[1] <= _gridSize[2]){
      _gridSize[1] *= prime_factors[i];
    }
    else{
      _gridSize[2] *= prime_factors[i];
    }
    i--;
  }
}

unsigned DomainDecomposition::Ndistribution(unsigned localN, float* minrnd, float* maxrnd)
{
   int num_procs;
   MPI_Comm_size(this->_commTopology, &num_procs);
   unsigned* moldistribution = new unsigned[num_procs];
   MPI_Allgather(&localN, 1, MPI_UNSIGNED, moldistribution, 1, MPI_UNSIGNED, this->_commTopology);
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

void DomainDecomposition::assertIntIdentity(int IX)
{
   if(this->_ownRank)
   {
      MPI_Send(&IX, 1, MPI_INT, 0, 2*_ownRank+17, this->_commTopology);
   }
   else
   {
      int recv;
      int num_procs;
      MPI_Comm_size(this->_commTopology, &num_procs);
      MPI_Status s;
      for(int i=1; i<num_procs; i++)
      {
         MPI_Recv(&recv, 1, MPI_INT, i, 2*i+17, this->_commTopology, &s);
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

void DomainDecomposition::assertDisjunctivity(TMoleculeContainer* mm)
{
   Molecule* m;
   global_log->set_mpi_output_root(0);

   if(this->_ownRank)
   {
      int tid;
      for(m = mm->begin(); m != mm->end(); m = mm->next())
      {
         tid = m->id();
         MPI_Send(&tid, 1, MPI_INT, 0, 2674+_ownRank, this->_commTopology);
      }
      tid = -1;
      MPI_Send(&tid, 1, MPI_INT, 0, 2674+_ownRank, this->_commTopology);
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
      MPI_Comm_size(this->_commTopology, &num_procs);
      MPI_Status s;
      for(int i=1; i < num_procs; i++)
      {
         bool cc = true;
         while(cc)
         {
            MPI_Recv(&recv, 1, MPI_INT, i, 2674+i, this->_commTopology, &s);
            if(recv == -1) cc = false;
            else
            {
               if(check.find(recv) != check.end()) {
                  global_log->error() << "Ranks " << check[recv] << " and " << i 
                              << " both propagate ID " << recv << endl;
                  exit(1);
               }
               else check[recv] = i;
            }
         }
      }
      global_log->info() << "Data consistency checked: No duplicate IDs detected among " << check.size()
           << " entries.\n";
   }
}
