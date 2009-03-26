#ifndef DOMAINDECOMPBASE_H_
#define DOMAINDECOMPBASE_H_

#include "utils/Log.h"
#include <list>
#include <vector>



class Molecule;
class Component;
class Domain;

namespace datastructures{
  template<class ParticleType>
  class ParticleContainer;
}

namespace parallel{
  class DomainDecompBase; 
}
using namespace std;

//! @brief handle boundary region and multiple processes
//!
//! This program is designed to run on a HPC (High Performance Computer).
//! But sometimes one might want to execute it on a single processor, possibly even
//! without having MPI installed on that machine. One way to allow this would be to
//! have two different versions of the program, one sequential version and one parallel version.
//! But that isn't feasible, as it is hardly possible to keep them both up to date without
//! investing a lot of additional time.
//! Befor describing how this problem is solved, you'll have to know a little bit about
//! how the parallel version works.
//!
//! At the moment, domain decomposition is used. The basic idea is, that the region which
//! is simulated (usually a cube) is divided into several smaller regions. Each process
//! works on one of those regions. But as the molecules move between regions, the processes
//! have to communicate to exchange molecules. Also for the calculation of some values
//! (e.g. macroscopic values), communication is necessary. Each time processes needs
//! to communicate, a method of this interface (that means a method of a class implementing
//! this interface) is called, which then somehow does the communication.
//!
//! Assume you have an class which implements this interface and which is capable of
//! doing all the necessary stuff for the parallelisation. Further assume that
//! you have a second class which also implements this interface which is only capable
//! of handling one process but doesn't need any MPI (with only one process, there is
//! no need for message passing between processes). So the main program (or in this
//! case the class Simulation) can decide which implementation to use. When MPI is
//! available, the parallel version is used, otherwise the sequential version
class parallel::DomainDecompBase{
 public:
  //! The Constructor sets up the topology and determines the ranks of the neighbours                                                       */
  DomainDecompBase() {};
  //DomainDecompBase(int *argc, char ***argv) {};

  //! The Destructor finalizes MPI
  virtual ~DomainDecompBase() {};

  //! returns the name of the Processor
  virtual const char* getProcessorName() const = 0;

  //! molecules which aren't in the domain of their process any
  //! more are transferred to their neighbours. Additionally, the
  //! molecules for the halo-region are transferred. To reduce the number
  //! of neighbours a single process has to communicate with, particles
  //! that i.e. have to be moved to the lower right neighbour are
  //! moved to the right neighbour first and then from the right neighbour
  //! to the lower neighbour.
  virtual void exchangeMolecules(datastructures::ParticleContainer<Molecule>* moleculeContainer, 
                                 const vector<Component>& components, Domain* domain) = 0;                                 
  
  
  // gather all molecules on the root process
//  virtual void gather_molecules(list<Molecule>& m_molecules, list<Molecule>& allMolecules, 
//      const double *proc_domain_L, const double *m_rmin, const vector<Component>& components) = 0;

  //! @brief appends molecule data to the file. The format is the same as that of the input file
  virtual void writeMoleculesToFile(string filename, datastructures::ParticleContainer<Molecule>* moleculeContainer) = 0;

  //! For several values, the sum over the whole domain has to be calculated
  //! This method takes the local sums of all processes and calculates the
  //! global sum. For some values, the global sum needs to be transferred to
  //! all processes, for others only to the root process
  virtual void reducevalues(double *Upot, double *Virial, double *summv2, double *sumIw2, unsigned long *num_mol) = 0;

 
  //! determines and returns the rank of the process at the given coordinates
  virtual int getRank(int x, int y, int z) = 0;

  //! returns the own rank
  virtual int getRank(void) = 0;

  //! gives the grid size of the given dimension (number of procs in one dim)
  virtual int getGridSize(int dimension) = 0;

  //! gives the own coordinates 
  virtual int getCoords(int dimension) = 0;

  //! synchronises all processes
  virtual void barrier() = 0;

  //! returns the time in seconds since some time in the past
  virtual double getTime() = 0;

 private:
 
  //! with the given number of processes, the dimensions of the grid are calculated
  //! grid_size: a array for the number of processes in each dimension
  virtual void setGridSize(int num_procs) = 0;

  //! Logging interface
  static utils::Log _log;

};


#endif /*DOMAINDECOMPBASE_H_*/
