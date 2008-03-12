#ifndef ADAPTIVESUBCELLS_H_
#define ADAPTIVESUBCELLS_H_

#include "utils/Log.h"

namespace datastructures {
  template<class ParticleType>
  class AdaptiveSubCells; 
  
  template<class ParticleType>
  class Cell;
}

//! @brief Adaptive SubCell Data Structure
//! @author Martin Buchholz
//! 
//! To understand the AdaptiveSubCells datastructure, you should first read the documentation
//! on the LinkedCells Datastructure, as AdaptiveSubCells is an enhancement of LinkedCells
//! The difference is that not all cells have the same size but the cell size is adapted to the
//! local density. The adaption only allows to refine once, as more refinements result in to
//! much overhead. So during the initialisation, first a regular cell structure is created
//! (Basically LinkedCells with l=rc, so big coarse cells).
//! Depending on the local density, the cells are refinded once (8 subcells) or remain coarse.
template<class ParticleType>
class datastructures::AdaptiveSubCells: public datastructures::ParticleContainer<ParticleType> {
  public:
    //#########################################################################
    //############# methods common for all ParticleContainers #################
    //#########################################################################

    //! @brief initialize the Adaptive SubCell datastructure
    //! 
    //! The constructor sets the following variables:
    //! - _cutoffRadius
    //! - _haloWidthInNumCells[3]
    //! - _cellsPerDimension[3]
    //! - _cellLength[3]
    //! - _haloBoundingBoxMin and _haloBoundingBoxMax
    //!
    //! It resized the _cells vector
    //! It resized _localRho and _metaCellIndex
    //! The corner parameters for the constructor describe the bounding box
    //! of the phasespace which belongs directly to this process, so they correspond
    //! to a bounding box including inner + boundary cells but excluding halo cells. \n
    //! But the corners of this class have to include the halo cells.
    //! @param bBoxMin lower corner of the bounding box of the domain belonging to this container
    //! @param bBoxMax higher corner of the bounding box of the domain belonging to this container 
    //! @param cutoffRadius distance for which forces have to be calculated
    //! @param cellsInCutoffRadius describes the width of the coarse cells relative to the cutoffRadius. 
    //!        This value should be 1, only then a useful adaption to the particle distribution can 
    //!        take place. The actual cell size is usually slightly bigger than the cutoffRadius,
    //!        as the domain has to be divided into a natural number of cells --> round up
    //! @param partPairsHandler specified concrete action to be done for each pair
    AdaptiveSubCells(double bBoxMin[3], double bBoxMax[3], double cutoffRadius, double cellsInCutoffRadius,
                datastructures::ParticlePairsHandler<ParticleType>& partPairsHandler);
    
    //! Destructor
    ~AdaptiveSubCells();
 
    //! Pointers to the particles are put into cells depending on the spacial position
    //! of the particles. 
    //! Before the call of this method, this distribution might have become invalid.
    //! To ensure, that all Particles (pointers to them) are put into the corresponding cells,
    //! first all cells are cleared and then filled again depending on the spacial position
    //! of the molecules. After the update, exactly one pointer for each particle in this
    //! ParticleContainer is in it's corresponding cell.   
    void update();
    
    //! @brief Insert Molecules into the cells.
    //!
    //! Loop over all Molecules and call this->addMolecule(...) for each of them
    //void addMolecules(std::list<ParticleType>& molecules);
    
    //! @brief Insert a single molecule.
    //!
    //! Therefore, first the cell (the index) for the molecule has to be determined,
    //! then the molecule is inserted into that cell.
    void addParticle(ParticleType& particle);

    //! @brief calculate the forces between the molecules.
    //!
    //! Only molecules with a distance not larger than the cutoff radius are to be used. \n
    //! Only forces on the Molecules which are in the inner and boundary region have to be calculated
    //! Newton's third law should be used for faster computation:
    //! \li a loop over all inner cells calculates the forces with all forward neighbour cells
    //!     all forward cells have to be used, as none of them can be halo or outside
    //! \li a loop over the boundary cells first calculates forces with all forward cells and all
    //!     backward cells. Here it has to be checked whether the neighbour cell is halo or not.
    //!     If it is Halo, the force is calculated, if it isn't, the force is not calculated, 
    //!     because the same pair of cells has already been processed in one of the other loops. 
    void traversePairs();

    //! @return the number of particles stored in the Linked Cells
    unsigned long getNumberOfParticles();
    
    //! @brief returns a pointer to the first particle in the Linked Cells
    //!
    //! Internally, the particles are store in a std::list. To traverse this
    //! list, a iterator (_particleIter) for the list is used.
    //! This method sets this iterator to point to the begin of the list
    //! and return a pointer to the value pointed to by the iterator
    ParticleType* begin();

    //! @brief returns a pointer to the next particle in the Linked Cells
    //!
    //! The iterator _particleIter is first incremented. Then a pointer
    //! to the value pointed to by the iterator is returned. If the
    //! iterator points to the end of the list (which is one element after the last
    //! element), NULL is returned 
    ParticleType* next();

    //! @brief returns NULL
    ParticleType* end();

    //! @brief delete all Particles which are not within the bounding box
    void deleteOuterParticles();

    //! @brief gets the width of the halo region in dimension index
    //! @todo remove this method    
    double get_halo_L(int index);    

  private:
    //####################################
    //######### PRIVATE METHODS ##########
    //####################################
    
    //! @brief Initialze index vectors and subCells. 
    //! 
    //! Fill the vector with the indices of the inner and boundary subCells.
    //! Assign each subCell it's region (halo, boundary, inner).
    void initializeSubCells();

    //! @brief Calculate SubCell neighbour indices.
    //! 
    //! In contrast to the LinkedCell datastructure, the index offsets to get to
    //! neighbouring cells are different for all cells. So for each cell (coarse or fine),
    //! a vector of index offets has to be calculated. This is done by first
    //! looping over all cells and then for each cell looping over all possible neighbour
    //! cells and storing those in the vector which have a minimal distance less than
    //! the cutoff radius. As an example, let's assume we want to calculate the index
    //! offset between cell A and cell B. There are basically four different cases:
    //! - Cell A fine and Cell B fine
    //! - Cell A fine and Cell B coarse
    //! - Cell A coarse and Cell B fine
    //! - Cell A coarse and Cell B coarse
    //! In each of those cases, the distance is calculated and where required the corresponding
    //! index offset (offset in the vector _subCells) is stored in _forwardNeighbourSubOffsets
    //! of _backwardNeighbourSubOffsets;
    void calculateSubNeighbourIndices(); 

    //! This method determines for a given Molecule the corresponding coarse cell
    //! and returns the index of that cell in the cell vector. \n
    //! If the molecule is not inside the bounding box, an error is printed
    int getCellIndexOfMolecule(ParticleType* molecule);
    
    //! @brief Get the index in the subCell vector to which this Molecule belong
    //!
    //! each spacial position within the bounding box of the adaptive linked cells
    //! belongs unambiguously to one subCell. \n
    //! This method determines for a given Molecule the corresponding subCell
    //! and returns the index of that subCell in the subCell vector. \n
    //! If the molecule is not inside the bounding box, an error is printed
    int getSubCellIndexOfMolecule(ParticleType* molecule);

    //! @brief given the 3D index of a subCell, return the index in the subCell vector.
    //!
    //! A cell can be identified by a 3D index. \n
    //! This method determines for a given 3D index the corresponding subCell
    //! and returns the index of that subCell in the subCell vector. \n
    //! Attention: The method can't(!!) be used to get the offset between two cells in the cell
    //! vector when called with the 3D cell index offets (e.g. x: one cell to the left,
    //! y: two cells back, z: one cell up,...). The offsets can differ from subCell to subCell, because both the number and the position of the Cells that are containing subCells can differ. 
    int subCellIndexOf3DIndex(int xIndex, int yIndex, int zIndex);

    //! @brief calculates the metaCellIndex for each coarse Cell and store it in _metaCellIndex
    //!
    //! Each coarse cell (or "meta" cell if the cell is refined) has also an index in
    //! the vector _subCells, which contains coarse and fine cells. metaCellOffset is
    //! the index offset of two coarse of "meta" cellsin the vector _subCells
    double calculateMetaCellIndex();

    //! @brief calculates the density for each coarse Cell separately
    //! @todo: currently, not the density but the number of particles
    //!        is calculated. To change that (or even if it is not changed), 
    //!        the influence of cutoff-radius has to be considered
    double calculateLocalRho();

    //####################################
    //##### PRIVATE MEMBER VARIABLES #####
    //####################################

    //! the list contains all molecules from the phasespace
    list<ParticleType> _particles;

    //! Iterator to traverse the list of particles (_particles)
    typename std::list<ParticleType>::iterator _particleIter;

    //! Vector containing all cells (including halo)
    //! @todo This vector is probably only needed during dynamic adaption
    //!       of the datastructure to calculate the local density
    std::vector<Cell<ParticleType> > _cells;

    //! Vector containing all subCells (including halo)
    std::vector<Cell<ParticleType> > _subCells;

    //! Vector containing the indices (for the subCells vector) of all inner subCells (without boundary)
    std::vector<int> _innerSubCellIndices;  
    //! Vector containing the indices (for the subCells vector) of all boundary subCells
    std::vector<int> _boundarySubCellIndices;

    //! Neighbours that come in the total ordering after a subCell
    std::vector<vector<int> > _forwardNeighbourSubOffsets;
    //! Neighbours that come in the total ordering before a subCell
    std::vector<vector<int> > _backwardNeighbourSubOffsets;

    //! low corner of the bounding box around the linked cells (including halo)
    double _haloBoundingBoxMin[3];
    //! high corner of the bounding box around the linkeaad cells (including halo)
    double _haloBoundingBoxMax[3];  

    //! Number of coarse cells in each spacial dimension (including halo)
    int _cellsPerDimension[3];
    //! Halo width (in cells) in each dimension
    int _haloWidthInNumCells[3];
    //! width of the halo strip (in lenght-unit)
    double _haloLength[3];
    //! length of one coarse cell (for each dimension)
    double _cellLength[3];
    //! cutoff radius
    double _cutoffRadius;
    
    //! Depending on the density, a cell is refined (resulting in 8 subcells) or not.
    //! All cells (fine and coarse) are stored in one big vector (_subCells). 
    //! E.g. Element 37 is a coarse cell (named A). It's neighbouring cell (B) is refined, 
    //! so in the vector _subCells, the next eight indices (38-45) are subcells of B.
    //! But there are loops (i from 0 to #coarseCells-1) which run over all coarse cells 
    //! (which might be refined or not).
    //! Those cells (and possibly subcells) have to be accessed, therefor their index
    //! in the vector _subcells has to be known. _metaCellIndex contains for each i 
    //! the correspoding index of the vector _subCells.
    std::vector<int> _metaCellIndex;

    //! Vector containing the density value of each coarse cell
    std::vector<double> _localRho;

    //! Number of times the update method was called. This value
    //! should approx. equal the simulation step. It is used to
    //! dynamically adapt the datastructure every nth time.
    int _numberOfUpdates;
    
    //! Logging interface
    static utils::Log _log;
      
};

#include "datastructures/AdaptiveSubCells.cpph"

#endif /*ADAPTIVESUBCELLS_H_*/
