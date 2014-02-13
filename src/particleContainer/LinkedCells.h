#ifndef LINKEDCELLS_H_
#define LINKEDCELLS_H_

#include <vector>

#include "particleContainer/ParticleContainer.h"
#include "particleContainer/ParticleCell.h"

//! @brief Linked Cell Data Structure
//! @author Martin Buchholz
//!
//! Without any specialized data structure, it needs O(N*N) - where N is
//! the number of particles - time to find all neighbouring pairs of particles.
//! The linked cell data structure is a datastructure which allows to find all
//! neighbouring pairs of particles (neighbouring means particles pairs
//! which have less than a certain distance) in O(N) time.
//! The following picture shows a domain with some particles in it. The blue
//! circle shows the neighbouring area of the red particle
//! \image particles.jpg
//! The problem is that all particles have to be examined to find those within the circle
//! With the linked cell data structure, the domain is divided into cells (using a regular grid).
//! All particles are placed in those cells.
//! For a given cell, neighbouring cells can easily be calculated, so for a given particle,
//! only the particles from neighbouring cells have to be examined. The following
//! picture illustrates this
//! \image particles_with_lc.jpg
//!
//! The spacial domain covered by the linked cells is larger then
//! the bounding box of the domain. This halo region surrounding
//! the phasespace is used for (periodic) boundary conditions
//! and has to be at least as wide as the cutoff radius. \n
//! In total, there are three different cell types:
//! - halo
//! - boundary
//! - inner

class LinkedCells : public ParticleContainer {
public:
	//! @brief initialize the Linked Cell datastructure
	//!
	//! The constructor sets the following variables:
	//! - _cutoffRadius
	//! - _LJCutoffRadius
	//! - _tersoffCutoffRadius
	//! - _haloWidthInNumCells[3]
	//! - _cellsPerDimension[3]
	//! - _cellLength[3]
	//! - _haloBoundingBoxMin and _haloBoundingBoxMax
	//!
	//! It resized the cell vector and marks the cells as inner/halo \n
	//! It fills the array innerCellIndices \n
	//! It fills the array with forward and backward neighbour indices \n
	//! The corner parameters for the constructor describe the bounding box
	//! of the phasespace which belongs directly to this process, so they correspond
	//! to a bounding box including inner + boundary cells but excluding halo cells. \n
	//! But the corners of this class have to include the halo cells.
	//! @param bBoxMin lower corner of the bounding box of the domain belonging to this container
	//! @param bBoxMax higher corner of the bounding box of the domain belonging to this container
	//! @param cutoffRadius distance for which forces have to be calculated
	//! @param LJCutoffRadius distance for which lennard jones forces have to be calculated
	//! @param tersoffCutoffRadius distance for which tersoff forces have to be calculated
	//! @param cellsInCutoffRadius describes the width of cells relative to the cutoffRadius: \n
	//!        equal (or larger) to the cutoffRadius divided by the length of a cell
	//!        as for the number of cells in each dimension only natural numbers are allowed,
	//!        it can happen that it is not possible to set celllength = cutoffRadius / cellsInCutoffRadius.
	//!        In that case, the celllength is chosen to be the next larger value so that the sum of
	//!        the cell lengths in one dimension equals the length of the phasespace
	//!        Example: phasespacelength=100, cellsInCutoffRadius=2, CutoffRadius=3 \n
	//!        ==> celllength should be: cutoffRadius/cellsInCutoffRadius = 3/2 = 1.5 \n
	//!        ==> cellsPerDimension = phasespacelength/celllength = 100/1.5 = 66.67 cells \n
	//!        ==> cells have to be larger: cellsPerDimension = phasespacelength/celllength = 100/celllength = 66 cells \n
	//!        ==> celllength = 100/66 = 1.5152
	LinkedCells(
		 double bBoxMin[3], double bBoxMax[3], double cutoffRadius, double LJCutoffRadius,
		  double cellsInCutoffRadius
	);

	//! Default constructor
	LinkedCells(){}
	//! Destructor
	~LinkedCells();

	virtual void readXML(XMLfileUnits& xmlconfig);

	double getHaloWidthNumCells(){
		return _haloWidthInNumCells[0];
	}

	// documentation see father class (ParticleContainer.h)
	void rebuild(double bBoxMin[3], double bBoxMax[3]);

	//! Pointers to the particles are put into cells depending on the spacial position
	//! of the particles.
	//! Before the call of this method, this distribution might have become invalid.
	//! To ensure, that all Particles (pointers to them) are put into the corresponding cells,
	//! first all cells are cleared and then filled again depending on the spacial position
	//! of the molecules. After the update, exactly one pointer for each particle in this
	//! ParticleContainer is it's corresponding cell.
	void update();

	//! @brief Insert a single molecule.
	//!
	//! Therefore, first the cell (the index) for the molecule has to be determined,
	//! then the molecule is inserted into that cell.
	void addParticle(Molecule& particle);

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
	//! @param particlePairsHandler specified concrete action to be done for each pair
//	void traversePairs(ParticlePairsHandler* particlePairsHandler);

	void traverseCells(CellProcessor& cellProcessor);

	//! @return the number of particles stored in the Linked Cells
	unsigned long getNumberOfParticles();

	//! @brief returns a pointer to the first particle in the Linked Cells
	//!
	//! Internally, the particles are store in a std::list. To traverse this
	//! list, a iterator (_particleIter) for the list is used.
	//! This method sets this iterator to point to the begin of the list
	//! and return a pointer to the value pointed to by the iterator
	Molecule* begin();

	//! @brief returns a pointer to the next particle in the Linked Cells
	//!
	//! The iterator _particleIter is first incremented. Then a pointer
	//! to the value pointed to by the iterator is returned. If the
	//! iterator points to the end of the list (which is one element after the last
	//! element), NULL is returned
	Molecule* next();

	//! @brief returns NULL
	Molecule* end();

	//! @brief deletes the current Molecule the iterator is at and returns the iterator to the next Molecule
	Molecule* deleteCurrent();

	//! @brief delete all Particles which are not within the bounding box
	void deleteOuterParticles();

	//! @brief gets the width of the halo region in dimension index
	//! @todo remove this method, because a halo_L shouldn't be necessary for every ParticleContainer
	//!       e.g. replace it by the cutoff-radius
	double get_halo_L(int index) const;

	//! @brief appends pointers to all particles in the halo region to the list
	void getHaloParticles(std::list<Molecule*> &haloParticlePtrs);

	// documentation see father class (ParticleContainer.h)
	void getRegion(double lowCorner[3], double highCorner[3], std::list<Molecule*> &particlePtrs);

	double getCutoff() { return _cutoffRadius; }
	void setCutoff(double rc) { _cutoffRadius = rc; }

	int getCellsInCutoff() { return _cellsInCutoff; }
	void setCellsInCutoff(int n) { _cellsInCutoff = n; }

	//! @brief counts all particles inside the bounding box
	unsigned countParticles(unsigned int cid);

	/**
	 * @todo move this method to the ChemicalPotential, using a call to ParticleContainer::getRegion() !?
	 */
	//! @brief counts particles in the intersection of bounding box and control volume
	unsigned countParticles(unsigned int cid, double* cbottom, double* ctop);

	void deleteMolecule(unsigned long molid, double x, double y, double z);
	/* TODO: The particle container should not contain any physics, search a new place for this. */
	double getEnergy(ParticlePairsHandler* particlePairsHandler, Molecule* m1, CellProcessor& cellProcessor);

	int localGrandcanonicalBalance() {
		return this->_localInsertionsMinusDeletions;
	}
	int grandcanonicalBalance(DomainDecompBase* comm);
	void grandcanonicalStep(ChemicalPotential* mu, double T, Domain* domain, CellProcessor& cellProcessor);

	double* boundingBoxMax() {
		return _boundingBoxMax;
	}

	double* boundingBoxMin() {
		return _boundingBoxMin;
	}

	int* boxWidthInNumCells() {
		return _boxWidthInNumCells;
	}

	double* cellLength() {
		return _cellLength; 
	}
	
#ifdef VTK
	friend class VTKGridWriter;
#endif



	friend class StatisticsWriter;


	//! @brief Get the index in the cell vector to which this Molecule belong
	//!
	//! each spacial position within the bounding box of the linked cells
	//! belongs unambiguously to one cell. \n
	//! This method determines for a given Molecule the corresponding cell
	//! and returns the index of that cell in the cell vector. \n
	//! If the molecule is not inside the bounding box, an error is printed
	unsigned long getCellIndexOfMolecule(Molecule* molecule) const;

	ParticleCell getCell(int idx){ return _cells[idx];}

private:
	//####################################
	//######### PRIVATE METHODS ##########
	//####################################

	//! @brief Initialize index vectors and cells.
	//!
	//! Fill the vector with the indices of the inner and boundary cells.
	//! Assign each cell it's region (halo, boundary, inner).
	void initializeCells();

	//! @brief Calculate neighbour indices.
	//!
	//! This method is executed once for the molecule container and not for
	//! each cell. E.g. the index (in the cell vector) of the right neighbour of a cell
	//! always equals the index of the cell minus one. This method calculates two vectors
	//! of index offsets, one for positive offsets (forward neighbours) and one for negative
	//! offsets (backward neighbours). So given a specific cell, the neighbours can be retrieved
	//! by adding to the index of the cell the offsets in the two vectors.
	//!
	//! The method works as follows: \n
	//! The loop runs over all potential neighbour cells (bounding box which contains
	//! the cell itself, and in each dimension on the lower and on the higher side as
	//! many cells as the width of the halo strip. E.g. if the haloWidth is 2, a box
	//! of 5x5x5 cell is considered as potential neighbours
	//! for each of those cells, the minimal possible distance between that cell
	//! and the central cell is calculated (sqrt(x^2+y^2+z^2)). If that distance
	//! is larger than the cutoff radius, the cell can be neglected.
	//! The distance in one dimension is the width of a cell multiplied with the number
	//! of cells between the two cells (this is received by substracting one of the difference).
	void calculateNeighbourIndices();

	

	//! @brief given the 3D index of a cell, return the index in the cell vector.
	//!
	//! A cell can be identified by a 3D index. \n
	//! This method determines for a given 3D index the corresponding cell
	//! and returns the index of that cell in the cell vector. \n
	//! The method can also be used to get the offset between two cells in the cell
	//! vector when called with the 3D cell index offsets (e.g. x: one cell to the left,
	//! y: two cells back, z: one cell up,...)
	long int cellIndexOf3DIndex(long int xIndex, long int yIndex, long int zIndex) const;



	//####################################
	//##### PRIVATE MEMBER VARIABLES #####
	//####################################

	std::list<Molecule> _particles; //!< List containing all molecules from the phasespace

	std::list<Molecule>::iterator _particleIter; //!< Iterator to traverse the list of particles (_particles)

	std::vector<ParticleCell> _cells; //!< Vector containing all cells (including halo)

	std::vector<unsigned long> _innerCellIndices; //!< Vector containing the indices (for the cells vector) of all inner cells (without boundary)
	std::vector<unsigned long> _boundaryCellIndices; //!< Vector containing the indices (for the cells vector) of all boundary cells
	std::vector<unsigned long> _haloCellIndices; //!< Vector containing the indices (for the cells vector) of all halo cells

	std::vector<long> _forwardNeighbourOffsets; //!< Neighbours that come in the total ordering after a cell
	std::vector<long> _backwardNeighbourOffsets; //!< Neighbours that come in the total ordering before a cell
	long _maxNeighbourOffset;
	long _minNeighbourOffset;

	double _haloBoundingBoxMin[3]; //!< low corner of the bounding box around the linked cells (including halo)
	double _haloBoundingBoxMax[3]; //!< high corner of the bounding box around the linked cells (including halo)

	int _cellsInCutoff; //!< Minimal number of cells within the cutoff radius
	int _cellsPerDimension[3]; //!< Number of Cells in each spacial dimension (including halo)
	int _haloWidthInNumCells[3]; //!< Halo width (in cells) in each dimension
	int _boxWidthInNumCells[3]; //!< Box width (in cells) in each dimension
	double _haloLength[3]; //!< width of the halo strip (in size units)
	double _cellLength[3]; //!< length of the cell (for each dimension)
	double _cutoffRadius; //!< RDF/electrostatics cutoff radius
    double _LJCutoffRadius;
	int _localInsertionsMinusDeletions; //!< balance of the grand canonical ensemble

	//! @brief True if all Particles are in the right cell
	//!
	//! The particles themselves are not stored in cells, but in one large
	//! list. The cells only contain pointers to particles. But at some
	//! points during the simulation, Those pointers are invalid. This happens
	//! e.g. after the integrator has changed the positions of particles.
	//! If this happens, no method must be able to access the particles via
	//! the cells. Therefore, whenever a piece of code causes the cells to
	//! become possibly invalid, _cellsValid has to been set to false. Methods
	//! accessing cells have to check whether _cellsValid is true (and e.g.
	//! abort the program if not). After the cells are updated, _cellsValid
	//! should be set to true.
	bool _cellsValid;

};

#endif /* LINKEDCELLS_H_ */
