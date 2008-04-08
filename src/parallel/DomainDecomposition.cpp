#include "parallel/DomainDecomposition.h"
#include "molecules/Molecule.h"
#include "datastructures/ParticleContainer.h"
#include "Domain.h"
using namespace std;


parallel::DomainDecomposition::DomainDecomposition(int *argc, char ***argv){

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
  
  MPI_Init(argc,argv);
  MPI_Get_processor_name(processorName, &procnamelen);
  MPI_Comm_size(world,&num_procs); // determine total number of procs 
  setGridSize(num_procs); // determine number of procs in each dimension 
  MPI_Cart_create(world,DIM,gridSize,period,reorder,&comm_topology); // create torus topology 
  
  MPI_Comm_rank(comm_topology, &ownrank);
  MPI_Cart_coords(comm_topology, ownrank, DIM, coords);
  neighbours[XLOWER]  = getRank(coords[0]-1,coords[1],coords[2]); 
  neighbours[XHIGHER] = getRank(coords[0]+1,coords[1],coords[2]); 
  neighbours[YLOWER]  = getRank(coords[0],coords[1]-1,coords[2]); 
  neighbours[YHIGHER] = getRank(coords[0],coords[1]+1,coords[2]); 
  neighbours[ZLOWER]  = getRank(coords[0],coords[1],coords[2]-1); 
  neighbours[ZHIGHER] = getRank(coords[0],coords[1],coords[2]+1);
}

parallel::DomainDecomposition::~DomainDecomposition(){
  MPI_Finalize(); 
}

const char* parallel::DomainDecomposition::getProcessorName() const {
  return processorName; 
}

void parallel::DomainDecomposition::exchangeMolecules(datastructures::ParticleContainer<Molecule>* moleculeContainer, 
                                 const vector<Component>& components, Domain* domain){

  double rmin[3]; // lower corner of the process-specific domain //PARALLEL
  double rmax[3];
  double halo_L[3]; // width of the halo strip //PARALLEL
  for(int i=0; i<3; i++){
    rmin[i] =  moleculeContainer->getBoundingBoxMin(i);
    rmax[i] =  moleculeContainer->getBoundingBoxMax(i);
    halo_L[i] =  moleculeContainer->get_halo_L(i);
  }


  Molecule* moleculePtr;
  std::vector<double> mol_rvqD[2]; // 13 elements of the vektor are one molecule: 3 pos, 3 v, 4 q, 3 D 
  std::vector<unsigned long> mol_id[2]; // 2 elements of the vektor are id and cid of a molecule
  std::vector<int> mol_cid[2]; // 2 elements of the vektor are id and cid of a molecule
  unsigned long numsend; // number of values that have to be sent 
  unsigned long numrecv; // number of values that have to be received 

  double *sendbuf_rvqD; // Array with values to be sent 
  double *recvbuf_rvqD; // Array for received values 
  unsigned long *sendbuf_id;
  unsigned long *recvbuf_id;
  int *sendbuf_cid;
  int *recvbuf_cid;

  double low_limit; // particles below this limit have to be copied or moved to the lower process 
  double high_limit; // particles above(or equal) this limit have to be copied or moved to the higher process 

  MPI_Status status;
  unsigned long i; // loop counter for all received molecules
  int direction; // direction (0=low/1=high) of molecule movement 

  for(unsigned short d=0;d<3;++d)
  {
    // set limits (outside "inner" region) 
    low_limit = rmin[d]+halo_L[d];
    high_limit = rmax[d]-halo_L[d];

    // when moving a particle across a periodic boundary, the molecule position has to change
    // these offset specify for each dimension (x, y and z) and each direction ("left"/lower 
    // neighbour and "right"/higher neighbour, how the paritcle coordinates have to be changed.
    // e.g. for dimension x (d=0) and a process on the left boundary of the domain, particles
    // moving to the left get the length of the whole domain added to their x-value
    double offsetLower[3];
    double offsetHigher[3];
    offsetLower[d] = 0.0;
    offsetHigher[d] = 0.0;
    
    // process on the left boundary
    if(coords[d] == 0){
      offsetLower[d] = domain->getGlobalLength(d);
    }
    // process on the right boundary
    if(coords[d] == gridSize[d]-1){
      offsetHigher[d] = -domain->getGlobalLength(d); 
    }

    moleculePtr = moleculeContainer->begin();
    while(moleculePtr!=moleculeContainer->end()){
      const double& rd=moleculePtr->r(d);
      if(rd<low_limit){ 
        // store the position of the molecule in a buffer
        for(unsigned short d2=0; d2<3; d2++){
          // when moving parallel to the coordinate d2 to another process, the
          // local coordinates in d2 change
          // @todo replace if-statement
          if(d2==d) mol_rvqD[0].push_back(rd+offsetLower[d]);
          else mol_rvqD[0].push_back(moleculePtr->r(d2));
        }
        // store velocity of the molecule in a buffer 
        for(unsigned short d2=0; d2<3; d2++) mol_rvqD[0].push_back(moleculePtr->v(d2));

        // store q's of the molecule in a buffer 
        mol_rvqD[0].push_back(moleculePtr->q().qw());
        mol_rvqD[0].push_back(moleculePtr->q().qx());
        mol_rvqD[0].push_back(moleculePtr->q().qy());
        mol_rvqD[0].push_back(moleculePtr->q().qz());


        // store D's of the molecule in a buffer 
        for(unsigned short d2=0; d2<3; d2++) mol_rvqD[0].push_back(moleculePtr->D(d2));

        // store id and componentid 
        mol_id[0].push_back(moleculePtr->id());
        mol_cid[0].push_back(moleculePtr->componentid());
        
        moleculePtr = moleculeContainer->next();
      }
      else if(rd>=high_limit){ 
        // store the position of the molecule in a buffer
        for(unsigned short d2=0; d2<3; d2++){
          // when moving parallel to the coordinate d2 to another process, the
          // local coordinates in d2 change    
           // @todo replace if-statement                               
          if(d2==d) mol_rvqD[1].push_back(rd+offsetHigher[d]);
          else mol_rvqD[1].push_back(moleculePtr->r(d2));
        }
        // store velocity of the molecule in a buffer 
        for(unsigned short d2=0; d2<3; d2++) mol_rvqD[1].push_back(moleculePtr->v(d2));

        // store q's of the molecule in a buffer
        mol_rvqD[1].push_back(moleculePtr->q().qw());
        mol_rvqD[1].push_back(moleculePtr->q().qx());
        mol_rvqD[1].push_back(moleculePtr->q().qy());
        mol_rvqD[1].push_back(moleculePtr->q().qz());

        // store D's of the molecule in a buffer 
        for(unsigned short d2=0; d2<3; d2++) mol_rvqD[1].push_back(moleculePtr->D(d2));

        // store id and componentid 
        mol_id[1].push_back(moleculePtr->id());
        mol_cid[1].push_back(moleculePtr->componentid());

        moleculePtr = moleculeContainer->next();
      }
      else moleculePtr = moleculeContainer->next();
    }


    // Communicate to lower and higher neighbour 
    for(direction=0; direction <=1; direction++){
      // Send to lower, receive from upper 
      // Send number of values that have to be sent 
      numsend = mol_id[direction].size();
      MPI_Sendrecv(&numsend, 1, MPI_UNSIGNED_LONG, neighbours[2*d+direction], 99, 
             &numrecv, 1, MPI_UNSIGNED_LONG, neighbours[2*d+(direction+1)%2], 99, comm_topology, &status);
    
      // initialize and fill send buffers 
      sendbuf_rvqD = new double[numsend*13];
      sendbuf_id = new unsigned long[numsend];
      sendbuf_cid = new int[numsend];
      for(i=0; i<(numsend*13); i++){
        sendbuf_rvqD[i] = mol_rvqD[direction][i];
      }
      for(i=0; i<numsend; i++){
        sendbuf_id[i] = mol_id[direction][i];
        sendbuf_cid[i] = mol_cid[direction][i];
      }
      // initialize receive buffer 
      recvbuf_rvqD = new double[numrecv*13];
      recvbuf_id = new unsigned long[numrecv];
      recvbuf_cid = new int[numrecv];
      
      // Send values to lower/upper and receive values from upper/lower 
      MPI_Sendrecv(sendbuf_rvqD, numsend*13, MPI_DOUBLE, neighbours[2*d+direction], 99, 
             recvbuf_rvqD, numrecv*13, MPI_DOUBLE, neighbours[2*d+(direction+1)%2], 99, comm_topology, &status);
      
      MPI_Sendrecv(sendbuf_id, numsend, MPI_UNSIGNED_LONG, neighbours[2*d+direction], 99, 
             recvbuf_id, numrecv, MPI_UNSIGNED_LONG, neighbours[2*d+(direction+1)%2], 99, comm_topology, &status);
      
      MPI_Sendrecv(sendbuf_cid, numsend, MPI_INT, neighbours[2*d+direction], 99, 
             recvbuf_cid, numrecv, MPI_INT, neighbours[2*d+(direction+1)%2], 99, comm_topology, &status);
      
      // insert received molecules into list of molecules 
      for(i=0; i<(numrecv); i++){
        Molecule m1 = Molecule(recvbuf_id[i],recvbuf_cid[i],recvbuf_rvqD[13*i],recvbuf_rvqD[13*i+1], recvbuf_rvqD[13*i+2], 
                    recvbuf_rvqD[13*i+3],recvbuf_rvqD[13*i+4],recvbuf_rvqD[13*i+5],
                    recvbuf_rvqD[13*i+6], recvbuf_rvqD[13*i+7],recvbuf_rvqD[13*i+8],recvbuf_rvqD[13*i+9],
                    recvbuf_rvqD[13*i+10], recvbuf_rvqD[13*i+11],recvbuf_rvqD[13*i+12], &components);
        moleculeContainer->addParticle(m1);
        //(m_molecules.back()).setFM(0.,0.,0.,0.,0.,0.);
      }

      // free memory 
      delete []sendbuf_rvqD;
      delete []sendbuf_id;
      delete []sendbuf_cid;
      delete []recvbuf_rvqD;
      delete []recvbuf_id;
      delete []recvbuf_cid;
      mol_rvqD[direction].clear();
      mol_id[direction].clear();
      mol_cid[direction].clear();
    }
  }
}


void parallel::DomainDecomposition::writeMoleculesToFile(string filename, datastructures::ParticleContainer<Molecule>* moleculeContainer){
  
  int numprocs;
  MPI_Comm_size(comm_topology, &numprocs);
  
  for(int process = 0; process < numprocs; process++){
    if(ownrank==process){
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



void parallel::DomainDecomposition::reducevalues(double *Upot, double *Virial, double *summv2, double *sumIw2, unsigned long *num_mol){
  double Upot_glob;
  double Virial_glob;
  double summv2_glob;
  double sumIw2_glob;
  unsigned long num_mol_glob;
  MPI_Allreduce(num_mol, &num_mol_glob, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm_topology);
  MPI_Allreduce(Upot, &Upot_glob, 1, MPI_DOUBLE, MPI_SUM, comm_topology);
  MPI_Allreduce(Virial, &Virial_glob, 1, MPI_DOUBLE, MPI_SUM, comm_topology);
  MPI_Allreduce(summv2, &summv2_glob, 1, MPI_DOUBLE, MPI_SUM, comm_topology);
  MPI_Allreduce(sumIw2, &sumIw2_glob, 1, MPI_DOUBLE, MPI_SUM, comm_topology);
  
  *Upot = Upot_glob;
  *Virial = Virial_glob;
  *summv2 = summv2_glob;
  *sumIw2 = sumIw2_glob;
  *num_mol = num_mol_glob;
}

void parallel::DomainDecomposition::setGridSize(int num_procs) {
  int remainder;      // remainder during the calculation of the prime factors 
  int i;              // counter                                               
  int num_factors;    // number of prime factors                               
  int *prime_factors; // array for the prime factors                           
  
  // Set the initial number of processes in each dimension to zero
  for(i=0;i<DIM;i++){
    gridSize[i] = 1;
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
  
  i = 0;
  while(i<num_factors){
    // The largest floor(num_factors/3) factors are the number of procs in x 
    if(i<num_factors/DIM) gridSize[0] *= prime_factors[num_factors-1-i];
    
    // The next factors are the number of procs in y 
    else if(i<2*num_factors/3) gridSize[1] *= prime_factors[num_factors-1-i];
    
    // The last ceil(n/3) factors are the number of procs in z 
    else gridSize[2] *= prime_factors[num_factors-1-i];
    
    i++;
  }
}


inline int parallel::DomainDecomposition::getRank(int x, int y, int z){
  int neigh_coords[DIM]; // Array for the coordinates 
  int neigh_rank;        // Rank of the neighbour     
  neigh_coords[0] = x;
  neigh_coords[1] = y;
  neigh_coords[2] = z;
  MPI_Cart_rank(comm_topology,neigh_coords,&neigh_rank);
  return(neigh_rank);
}

double parallel::DomainDecomposition::getTime(){
  return MPI_Wtime();
}

