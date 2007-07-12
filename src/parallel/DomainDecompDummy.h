#ifndef DOMAINDECOMPDUMMY_H_
#define DOMAINDECOMPDUMMY_H_


#include "utils/Log.h"
//#include "datastructures/MoleculeContainer.h"
#include "parallel/DomainDecompBase.h"

namespace parallel{
  class DomainDecompDummy; 
}

class Molecule;

//! @brief implement the %domain decomposition for a single processor  
//! 
//! As explained in the comment to DomainDecompBase, there is only one
//! code for the sequential and the parallel version. This class 
//! accomplished the decomposition for the sequential version (1 process)
//! 
//! While this class doesn't has to implement communication between processes,
//! there are still some other things to do, especially the handling of the
//! boundary region.  
class parallel::DomainDecompDummy: public parallel::DomainDecompBase {
 public:
  //! The constructor has nothing to do
  DomainDecompDummy();
  
  //! The destructor has nothing to do
  ~DomainDecompDummy();

  //! for the sequential version, the processor name is not that important
  //! @todo implement this
  const char* getProcessorName() const;

  //! molecules which aren't in the domain of their process any
  //! are moved to the opposite side of the domain (periodic boundary).
  //! Additionally, the molecules from the boundary region are duplicated
  //! and copied into the corresponding halo-regions.
  void exchangeMolecules(datastructures::ParticleContainer<Molecule>* moleculeContainer, 
                                 const vector<Component>& components, Domain* domain);
  
  //! @brief opends the file(append), loops over all molecules and writes each into the file                          
  void writeMoleculesToFile(string filename, datastructures::ParticleContainer<Molecule>* moleculeContainer);

  //! @brief There is nothing to be reduced
  //!
  //! In the sequential version there is only one process, so the local values are
  //! equal to the global values. So there's nothing to be done here.
  void reducevalues(double *Upot, double *Virial, double *summv2, double *sumIw2, unsigned long *num_mol);

 
  //! @brief There is only one process, so this method always returns 0
  int getRank(int x, int y, int z);

  //! @brief There is only one process, so this method always returns 0
  int getRank(void);

  //! @brief As there is only one process in each dimension, this method always returns 1
  int getGridSize(int dimension);

  //! @brief As there is only one process in each dimension, this method always returns 0
  int getCoords(int dimension);

  //! @brief one process doesn't need synchronisation, so nothing is done here
  void barrier();
  
  double getTime();

 private:
 
  /* with the given number of processes, the dimensions of the grid are calculated *
   * grid_size: a array for the number of processes in each dimension              */
  void setGridSize(int num_procs);

  //! Logging interface
  static utils::Log _log;

};


#endif /*DOMAINDECOMPDUMMY_H_*/
