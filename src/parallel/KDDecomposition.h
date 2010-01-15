#ifndef KDDECOMPOSITION_H_
#define KDDECOMPOSITION_H_

#define KDDIM 3

#include <iostream>
#include <list>
#include <mpi.h>

#include "parallel/CollectiveCommunication.h"
#include "DomainDecompBase.h"

using namespace std;

class Molecule; 
class ParticleData;
class KDNode;

//#######################################################
//### KDDecomposition class                           ###
//#######################################################

//! @brief domain decomposition based on a kind of KD-Tree
//! @author Martin Buchholz
//! 
//! There are two main goals which shall be achieved by this class.
//! First, it should be able to distribute the load evenly among the
//! processes (so create a load-balanced decomposition). The second goal
//! is to do that in a way so that the domain boundaries are in region
//! which are not dense (e.g. gas). This reduces the number of 
//! necessary double computations and the amount of communication.
//! 
//! The technique used is a modified KD-Tree. First, the region is discretised
//! into cells (size equals (or larger then) cutoffradius). This region is then
//! divided into two smaller cuboids. Those cuboids are again recursively divided
//! into smaller cuboids. Three questions arise:
//! - in which dimension (x,y,z) is the cube divided in each step ?
//! - where is the "cut" (how many cell-layers on the "left" and on the "right") ?
//! - how many processes share the "left" and the "right" part ?
//! This class tries to find an answer to this three questions which balances the
//! load and and the same time tries to minimise the communication (and double 
//! calculation) costs. This is done by "guessing" the overall costs/time for
//! each possible division and selection the best. A guess basically works as follows:
//! for one division, the computation costs for each part are calculated (using
//! e.g. number of particles, number of pairs, some heuristics, ...). For each
//! part, these costs are divide by the number of procs which share this part
//! (here again the best combination is used). Additionally, the costs for the "cut"
//! are calculated (again heuristics e.g. depending on the number of pairs which are cut)
//! With all the calculated costs, one can select that division which is cheapest.
//! This is done recursively until each process has its own region.
//! @todo CURRENTLY, THIS CLASS IS IN DEVELOPMENT
class KDDecomposition: public DomainDecompBase{
 public:
  //! @brief create an initial decomposition tree
  //!
  //! The constructor has to determine the own rank and the number of neighbours.
  //! It has to determine the number of cells and create an initial decomposition
  //! of the domain (no knowledge about particles yet), which is stored in
  //! _decompTree and _ownArea
  KDDecomposition(double cutoffRadius, Domain* domain);

  // documentation see father class (DomainDecompBase.h)
  ~KDDecomposition();

  //###############################################
  //### The following methods are those of the  ### 
  //### base class which have to be implemented ###
  //###############################################

  //! @brief exchange molecules between processes
  //!
  //! molecules which aren't in the domain of their process any
  //! more are transferred to their neighbours. Additionally, the
  //! molecules for the halo-region are transferred. 
  //! In this implementation, the methods used for load-balancing and those
  //! for just exchanging particles without rebalancing are quite similar.
  //! Therefore, this method just calls balanceAndExchange(0,...), where
  //! the "0" says that only exchanging and no balancing has to be done.
  //! @param moleculeContainer needed to get those molecules which have to be exchanged
  //! @param components when creating a new Molecule-object (from the recieved data), 
  //!                   the Molecule-constructor needs this component vector
  //! @param domain is e.g. needed to get the size of the local domain  
  void exchangeMolecules(ParticleContainer* moleculeContainer, const vector<Component>& components, Domain* domain, double rc);

  //! @brief balance the load (and optimise communication) and exchange boundary particles
  //! 
  //! The workflow is as follows (some steps only if balancing shall be done):
  //! - if balance: _numParticlesPerCell is updated to be able to calculate the load of the cells
  //! - if balance: A new decomposition tree is calculated (using recDecomp) based on the particle distribution
  //! - if balance: find out neighbours
  //! - collect pointers to particles to be send to all neighbours
  //! - communicate with each neighbour how many particles will be sent and recieved
  //! - copy data to be sent into buffers
  //! - if balance: rebuild the molecule container (new size and region), replace the old decomposition tree
  //! - transfer data and insert the recieved molecules into the moleculeContainer
  //! - for processes which span the whole domain in at least one direction, ensure periodic boundary
  //! @param balance if true, a rebalancing should be performed, otherwise only exchange
  //! @param moleculeContainer needed for calculating load and to get the particles
  //! @param components when creating a new Molecule-object (from the recieved data), 
  //!                   the Molecule-constructor needs this component vector
  //! @param domain is e.g. needed to get the size of the local domain
  void balanceAndExchange(bool balance, ParticleContainer* moleculeContainer,
                          const vector<Component>& components, Domain* domain, double rc);

  // documentation see father class (DomainDecompBase.h)
  bool procOwnsPos(double x, double y, double z, Domain* domain);
  
  // documentation see father class (DomainDecompBase.h)
  double guaranteedDistance(double x, double y, double z, Domain* domain);

  // documentation see father class (DomainDecompBase.h)
  int countMolecules(ParticleContainer* moleculeContainer, vector<int> &compCount);

  //! @todo comment and thing
  double getBoundingBoxMin(int dimension, Domain* domain);
  //! @todo comment and thing
  double getBoundingBoxMax(int dimension, Domain* domain);

  //! @brief writes information about the current decomposition into the given file
  //! 
  //! This decomposition first writes a very small header with the following lines:
  //! - "size", followed by three double values for the size of the domain
  //! - "decompData Regions": just this text to state that now the header is finished and data follows
  //! - one line per proc, giving the bounding box (minx, miny, minz, maxx, maxy, maxz)
  //! 
  //! An example file for 5 procs could look like this:
  //!
  //! |size 62.0 62.0 62.0 \n
  //!  decompData Regions \n
  //!  0.0 0.0 0.0 20.0 62.0 62.0 \n
  //!  20.0 0.0 0.0 62.0 40.0 25.0 \n
  //!  20.0 40.0 0.0 62.0 62.0 25.0 \n
  //!  20.0 0.0 25.0 62.0 30.0 62.0 \n
  //!  20.0 30.0 25.0 62.0 62.0 62.0
  //! @param filename name of the file into which the data will be written
  //! @param domain e.g. needed to get the bounding boxes
  void printDecomp(string filename, Domain* domain);

  //! @brief append the molecule date of all processes to the file
  //!
  //! Currently, parallel IO isn't used.
  //! To ensure that not more than one process writes to the file at any time,
  //! there is a loop over all processes with a barrier in between
  //! @param filename name of the file into which the data will be written
  //! @param moleculeContainer all Particles from this container will be written to the file
  void writeMoleculesToFile(string filename, ParticleContainer* moleculeContainer);

  // documentation see father class (DomainDecompBase.h)
  int getRank(void){ return _ownRank;}
  
  // documentation see father class (DomainDecompBase.h)
  int getNumProcs(){ return _numProcs;}

  // documentation see father class (DomainDecompBase.h)
  const char* getProcessorName() const { return _processorName;};

  // documentation see father class (DomainDecompBase.h)
  void barrier() {MPI_Barrier(MPI_COMM_WORLD);}

  // documentation see father class (DomainDecompBase.h)
  double getTime() { return MPI_Wtime(); };

  //##################################################################
  // The following methods with prefix "collComm" are all used 
  // in the context of collective communication. Each of the methods
  // basically has to call the corresponding method from the class
  // CollectiveCommunication (or CollectiveCommDummy in the sequential
  // case). To get information about how to use this methods, read 
  // the documentation of the class CollectiveCommunication and of the
  // father class of this class (DomainDecompBase.h)
  //##################################################################  
  void collCommInit(int numValues){ _collComm.init(MPI_COMM_WORLD, numValues); };
  void collCommFinalize(){ _collComm.finalize(); };
  void collCommAppendInt(int intValue){_collComm.appendInt(intValue);};
  void collCommAppendUnsLong(unsigned long unsLongValue){_collComm.appendUnsLong(unsLongValue);};
  void collCommAppendFloat(float floatValue){_collComm.appendFloat(floatValue);};
  void collCommAppendDouble(double doubleValue){_collComm.appendDouble(doubleValue);};
  void collCommAppendLongDouble(long double longDoubleValue){_collComm.appendLongDouble(longDoubleValue);};
  int collCommGetInt(){return _collComm.getInt(); };
  unsigned long collCommGetUnsLong(){return _collComm.getUnsLong(); };
  float collCommGetFloat(){return _collComm.getInt(); };
  double collCommGetDouble(){ return _collComm.getDouble(); };
  long double collCommGetLongDouble(){ return _collComm.getLongDouble(); };
  void collCommAllreduceSum(){ _collComm.allreduceSum(); };
  void collCommBroadcast(){ _collComm.broadcast(); };
  
  //! @brief returns total number of molecules
  unsigned Ndistribution(unsigned localN, float* minrnd, float* maxrnd);

  //! @brief checks identity of random number generators
  void assertIntIdentity(int IX);
  void assertDisjunctivity(TMoleculeContainer* mm);

 private:
  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  //$ Methoden, die von balanceAndExchange ben�tigt werden $
  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  //! @brief Each process collects the particles to be send to neighbours
  //!
  //! First, all processes (and their region, incl. halo) which need  
  //! part of the sourceArea are determined. 
  //! (Also all necessary periodic copies of the neighbours region)
  //! Then all particles (pointers to them) from those regions
  //! are stored in particlesToSend[procid]
  //! The area given by sourceArea and the area covered by the moleculeContainer
  //! should be the same!
  //! @param sourceArea area from which particles should be sent
  //! @param moleculeContainer needed to get the molecule pointers
  //! @param domain needed to get the domain size
  //! @param procIDs Here the ids of the neighbouring procs will be stored
  //! @param numMolsToSend Here the number of molecules to be sent are stored 
  //!                      The vector has to be initialised in this method
  //! @param particlesToSend Here the pointers to the particles will be stored
  void getPartsToSend(KDNode* sourceArea, KDNode* decompTree, ParticleContainer* moleculeContainer, Domain* domain, vector<int>& procIDs, vector<int>& numMolsToSend, vector<list<Molecule*> >& particlesToSend);

  //! @brief Neighbouring procs inform each other about number of particles
  //! 
  //! Before particle data can be transferred, the recieving procs have to
  //! allocate memory to store those data. Therefore they have to be informed
  //! about the number of particles which have to be recieved.
  //! @param procsToSendTo ids of all neighbours (possibly including self) which have
  //!        to recieve particles from this process
  //! @param procsToRecvFrom ids of all neighbours (possibly including self) which have
  //!        to send particles to this process  
  //! @param particlesToSend Needed to know the number of particles to send
  //! @param numMolsToRecv Array to store the number of molecules which have to be recieved.
  void exchangeNumToSend(vector<int>& procsToSendTo, vector<int>& procsToRecvFrom, vector<int>& numMolsToSend, vector<int>& numMolsToRecv);

  //! @brief transfer of the molecule data to the neighbours
  //! 
  //! After each process knows which particles have to be sent
  //! and how many particles have to be received, and all buffers
  //! have been initialised, the data transfer itself can be done
  //! using this method
  //! @param procIDs ids of all neighbours (possibly including self)
  //! @param numMolsToSend Number of molecules to be sent to each neighbour
  //! @param numMolsToRecv Number of molecules to be recieved from each neighbour
  //! @param particlesSendBufs already filled buffers of particles to be sent
  //! @param particlesRecvBufs buffers (already allocated) to store recieved particles
  void transferMolData(vector<int>& procsToSendTo, vector<int>& procsToRecvFrom, vector<int>& numMolsToSend, vector<int>& numMolsToRecv, vector<ParticleData*>& particlesSendBufs, vector<ParticleData*>& particlesRecvBufs);

  //! @brief create copies due to periodic boundary
  //!
  //! If a proc's area occupies the whole domain (per dimension), it has two periodic
  //! boundaries in that dimension. That means that is possesses at least two (up to eight)
  //! versions of the same cells/molecules. Those copies are created by this method
  //! and stored in the moleculeContainer
  //! @param moleculeContainer used to walk through the molecules, create copies and store them again
  //! @param domain needed to get the length of the domain
  //! @param components needed to create new molecules
  //! @todo make it work with overlapping decomposition trees
  //! @todo more efficiency (don't run over all molecules)
  void createLocalCopies(ParticleContainer* moleculeContainer, Domain* domain, const vector<Component>& components);

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  //$ sonstige Methoden
  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  //! @brief prints the given (sub-) tree to stdout
  //!
  //! For each node, it is printed whether it is a "LEAF" or a "INNER" node,
  //! The owner is printed as well as the region of each node
  //! The order of printing is a depth-first walk through the tree, children
  //! are always indented two spaces more than there parents
  //! @param root Node for which the subtree shall be printed
  //! @param prefix A string which is printed in front of each line
  void printDecompTree(KDNode* root, string prefix);

  //! @brief provides a "real" modulo function
  //!
  //! The C build-in modulo operation uses a symmetric modulo,
  //! This method provides a mathematical modulo (important for negative numbers)
  //! @param number Number to calculate the modulo for
  //! @param modulo obviously the modulo
  int mod(int number, int modulo);


  //! @brief create an initial decomposition of the domain
  //!
  //! At the beginning of the simulation, particles have to be read in.
  //! At this time, the distribution of particles is unknown, so an
  //! an initial domain decomposition has to be created.
  //! This method does achieve this with a very simple decomposition
  //! It has to be recursively called for all inner nodes
  //! A inner node is a node for which numProcs is larger than one.
  //! So if numProcs is larger than one, two new children are created
  //! If the children are inner nodes, recInitialDecomp is recursively
  //! called.
  //! @param fatherNode Pointer to a (already existing!) inner node
  //! @param ownArea This pointer has to be set to the leaf node 
  //!                representing the own area
  void recInitialDecomp(KDNode* fatherNode, KDNode*& ownArea);
  void recDecomp(KDNode* fatherNode, KDNode*& ownArea);
    
  unsigned long predictBoundariesPerProc(int numcells, int numprocs);
  


  //! @brief Find all processes which own a part of the given area
  //!
  //! This method is supposed to run recursively. It should usually be called
  //! with decompTree and testNode pointing both to the root of a decomposition tree.
  //! In the recursive calls, testNode will be changed to walk through the whole tree.
  //! The given area can be partly outside of the area covered by the decomposition tree
  //! (usually because of the halo region), in this case the periodic boundary conditions
  //! are applied to get the corresponging area on the other side of the region.
  //! @param low low corner of the area for which the owners shall be found
  //! @param high high corner of the area for which the owners shall be found
  //! @param decompTree root KDNode of a domain decomposition tree. The area covered by
  //!                   this KDNode must not overlap the domain
  //! @param testNode KDNode currenty tested for intersection (method is used recursively)
  //! @param procIDs here the ids of those procs which own a part will be stored
  //! @param neighbHaloAreas here the corresponging (see procIDs) areas will be stored. For
  //!                        each procID six double values are stored, first the three
  //!                        coordinates for the low corner of the area, then the high corner
  //! @todo make it work for overlapping decomposition trees
  void getOwningProcs(int low[KDDIM], int high[KDDIM], KDNode* decompTree, KDNode* testNode, vector<int>* procIDs, vector<int>* neighbHaloAreas);


  //! @brief 
  //! @todo _numParticles should perhaps not be a member variable (think about that)
  void getNumParticles(ParticleContainer* moleculeContainer);
  void newDecomp(ParticleContainer* moleculeContainer);

  //######################################
  //###    private member variables    ###
  //######################################

  //! Rank of the local process (to be set by the constructor)
  int _ownRank;
  //! total number of processes (to be set by the constructor)
  int _numProcs;
  //! name of the local processor (to be set by the constructor)
  char _processorName[MPI_MAX_PROCESSOR_NAME];



  //! Length of one cell. The length of the cells is the cutoff-radius, or a 
  //! little bit larger as there has to be a natural number of cells.
  double _cellSize[KDDIM];
  
  //! Number of cells in the global domain in each coordinate direction.
  //!  _globalCellsPerDim does not include the halo region
  int _globalCellsPerDim[KDDIM];
  //! global number of cells (without halo)
  int _globalNumCells;

  // Pointer to the root element of the initial decomposition tree
  KDNode* _decompTree;
  // each process own an area in the initial decomposition
  KDNode* _ownArea;
  //! during balancing, _decompTree and _ownArea become invalid.
  //! They might already point to the new decomposition/Area, but the
  //! moleculeContainer is not yet updated e.g.
  bool _decompValid; 



  // Number of particles for each cell (including halo?)
  int *_numParticlesPerCell;

  //! variable used for different kinds of collective operations
  CollectiveCommunication _collComm;
  
  // local linked cells list
  //vector<Cell> _localCells;
  //$ noch zu verbessern (static,?)
  //double _timePerBoundary;
  //double _timePerForceCalc;
};


#endif /*KDDECOMPOSITION_H_*/






//  void printArea(KDNode* area);


/*   // when readInData is called, each proc reads in his data from the file */
/*   // corresponding to the decomposition defined by _initDecompTree and */
/*   // _initOwnArea, which is the leaf node of _initDecompTree owned by this process; */
/*   // The Halo region stays empty and has to be filled by some other method! */
/*   void readInData(); */

//  int getOwningProcPerBound(int cellPos[KDDIM], KDNode* decompTree);

// den Prozess finden, dem die Zelle geh�rt
// die Koordinaten cellPos m�ssen bez�glich des selben Referenzrahmens
// angegeben werden, der auch von decompTree verwendet wird.
// int getOwningProc(int cellPos[KDDIM], KDNode* decompTree);


//! evtl. Verschiebung des Referenzrahmens muss kompensiert werden
//! Angenommen man legt ein Zellgitter �ber das "reale" Gebiet, so dass
//! die unterste Zelle genau im untersten Eck des Gebiets liegt. Nun
//! deckt das durch _initDecompTree abgedeckte Gebiet zwar genau das gleiche
//! Gebiet ab, aber die unterste Zelle ist gegen�ber dem "realen" Zellgitter
//! verschoben um _initOffset (m�glich durch periodischen Rand). Um nun zu pr�fen, ob eine
//! Zelle in dem Gebiet _initOwnArea liegt, muss also zun�chst die Koordinate
//! angepasst werden. Um die Position bez�glich des neuen Gitters zu erhalten,
//! muss von der urspr�nglichen Koordinate _initOffset abgezogen werden.
//int _initOffset[KDDIM];
