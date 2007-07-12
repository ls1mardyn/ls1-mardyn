/*****************************************************
 * This class represents the whole simulation domain *
 * which is divided into n smaller domains. n is the *
 * number of processes involved in the parallel      *
 * execution. The methods perform operations which   *
 * involve more than one domain                      *
 *****************************************************/
//! The domain is the spacial area which can be reached by the molecules
//! during the simulation. At the start of the simulation, the information
//! about position, velocity,... of all molecules is stored in the 
//! phasespace file. The task of this class is to read in all information
//! and distribute it to the processes (or let each process read it's own information)
//! In a domain decomposition, the domain is divided into several parts and each process
//! is assigned one of those subdomains. For one time-step, each process works on this
//! subdomain as if it was the whole domain. The different sub-domains are
//! overlapping a little bit, because at the boundary of a domain, neighbouring molecules
//! are needed for the force calculation. After one time step, there has to be some
//! Communication along the boundary with the processes which work on neighbouring domains.
//!
//! On each process, there is an instance of the class Phasespace which has to hold all
//! information about molecules which are within the process' domain.
//!
//! The phasespace file looks like this:
//! MDProject \n
//! phasespaceSize "x (double)" "y (double)" "z (double)" \n
//! timestepLength "dt (double)" \n
//! currentTime "t (double)" \n
//! cutoffRadius "rc (double)" \n
//! NumberOfMolecules "N (int)" \n
//! [id1] [type1] [x1] [y1] [z1] [vx1] [vy1] [vz1] \n
//! [id2] [type2] [x2] [y2] [z2] [vx2] [vy2] [vz2] \n
//! ... \n
//! [idN] [typeN] [xN] [yN] [zN] [vxN] [vyN] [vzN] \n 

#ifndef DOMAINDECOMPOSITION_H_
#define DOMAINDECOMPOSITION_H_

#define DIM 3

/* indexes of the neighbours in the neighbour-array */
#define XLOWER 0
#define XHIGHER 1
#define YLOWER 2
#define YHIGHER 3
#define ZLOWER 4
#define ZHIGHER 5

#include <iostream>

// See page 27 of MPICH2 User Guide for reference
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#include <mpi.h>
#include <math.h>

#include "parallel/DomainDecompBase.h"

#include "mpi.h"
using namespace std;

namespace parallel{
  class DomainDecomposition; 
}

class Molecule;

//! @brief implements a %domain decomposition for the parallel version
//!
//! In a domain decomposition, each process gets part of the spacial domain.
//! In this implementation, the whole domain has to be a cuboid which is
//! decomposed into several smaller cuboids.
//! The main difficulty is that at the boundary, each process needs molecules
//! from neighbouring domains to be able to calculate the forces on the own
//! molecules. Another problem is that molecules are moving across the
//! boundaries of local domains. So methods are implemented to transfer those molecules
//! @todo reference to a paper describing the %domain decomposition
class parallel::DomainDecomposition: public parallel::DomainDecompBase{
 public:
  //! The Constructor sets up the topology and determines the ranks of 
  //! the neighbours
  DomainDecomposition(int *argc, char ***argv);

  //! The Destructor finalizes MPI
  ~DomainDecomposition();

  //! Processor name
  const char* getProcessorName() const;
  
  //! molecules which aren't in the domain of their process any 
  //! more are transferred to their neighbours. Additionally, the
  //! molecules for the halo-region are transferred. To reduce the number
  //! of neighbours a single process has to communicate with, particles
  //! that i.e. have to be moved to the lower right neighbour are
  //! moved to the right neighbour first and then from the right neighbour
  //! to the lower neighbour.
  void exchangeMolecules(datastructures::ParticleContainer<Molecule>* moleculeContainer, 
                          const vector<Component>& components, Domain* domain);
  
  //! @brief append the molecule date of all processes to the file
  //!
  //! Currently, parallel IO isn't used.
  //! To ensure that not more than one process writes to the file at any time,
  //! there is a loop over all processes with a barrier in between
  void writeMoleculesToFile(string filename, datastructures::ParticleContainer<Molecule>* moleculeContainer);

  //! For several values, the sum over the whole domain has to be calculated
  //! This method takes the local sums of all processes and calculates the
  //! global sum. For some values, the global sum needs to be transferred to
  //! all processes, for others only to the root process
  void reducevalues(double *Upot, double *Virial, double *summv2, double *sumIw2, unsigned long *num_mol);

 
  //! determines and returns the rank of the process at the given coordinates
  int getRank(int x, int y, int z);

  //! returns the own rank
  int getRank(void){ return ownrank;}

  //! gives the grid size of the given dimension (number of procs in one dim)
  int getGridSize(int dimension) { return gridSize[dimension];}

  //! gives the own coordinates
  int getCoords(int dimension) { return coords[dimension];}

  //! synchronises all processes
  void barrier() {MPI_Barrier(comm_topology);}

  double getTime();

 private:
  
  //! new topology after initializing the torus
  MPI_Comm comm_topology;
  //! Number of processes in each dimension (i.e. 2 for 8 processes) 
  int gridSize[DIM];
  //! Grid coordinates of process
  int coords[DIM];
  //!  rank of process
  int ownrank;
  //! rank of neighbours starting with the lowest dimension (x) and
  //! the lower value (x-1; "left), then the higher value (x+1).
  //! Then follows the next dimension (y), so all elements are:
  //! (x-1,y,z),(x+1,y,z),(x,y-1,z),(x,y+1,z),(x,y,z-1),(x,y,z+1)
  int neighbours[2*DIM];
  //! name of the local processor
  char processorName[MPI_MAX_PROCESSOR_NAME];

  //! with the given number of processes, the dimensions of the grid are calculated
  void setGridSize(int num_procs);


};

#endif /*DOMAINDECOMPOSITION_H_*/
