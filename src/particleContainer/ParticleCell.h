#ifndef PARTICLE_CELL_H_
#define PARTICLE_CELL_H_

#include "Cell.h"

#include <vector>

class Molecule;

//! @brief ParticleCell data structure.
//! @author Martin Buchholz
//!
//! A ParticleCell represents a small cuboid area of the domain and stores a list of 
//! pointers to the molecules in that area. Depending on the actual position
//! of the cell, it belongs to one of four different regions: \n
//! - completele outside (more than the cutoffradius away from all cells that
//!                       belong directly to the MoleculeContainer)
//!                       Such a cell shouldn't exist!!!
//! - halo region (not belonging directly to the MoleculeContainer, but within
//!                       the cutoffradius of at least one cell of the MoleculeContainer)
//! - boundary region (belonging directly to the MoleculeContainer, but not more than
//!                           the cutoffradius away from the boundary)
//! - inner region (more than the cutoffradius away from the boundary)
//! 
//! There are three boolean member variables for the last three regions. \n
//! If more than one of them is true, there must be an error in the code \n
//! If none of them is true, the cell wasn't assigned to any region yet.
//! A cell which is completely outside shouldn't exist, as it completely useless.
//! 
//! For cuboid domains, a cell doesn't need much information
//! (you could even avoid a seperate class completely)
//! Most information (position, neighbours,...) can be calculated
//! from the cell's index in the linked cell vector.
//! e.g. for a domain of 10x10x10 cells, the cell with index (counting from zero)
//! 561 has spacial index (counting from zero)  (1,6,5) and the cell indices of
//! the neighbours are e.g. 450, 451, 452, 460, 461, ... 
//! But for either unregular domains or unregular spacial decompositions
//! in the parallelisation (e.g. space filling curves), it's not easy to store that
//! information in the linked cell vector (if it is a vector at all)
//! That's why this class was introduced, so additional information which is
//! specific to a cell should be stored here.


class ParticleCell : public Cell {
public:

	//! removes all elements from the list molecules
	void removeAllParticles();

	//! insert a single molecule into this cell
	void addParticle(Molecule* particle_ptr);

	//! return a reference to the list of molecules (molecule pointers) in this cell
	std::vector<Molecule*>& getParticlePointers();

	bool deleteMolecule(unsigned long molid);

	//! return the number of molecules contained in this cell
	int getMoleculeCount() const;
	
private:
	std::vector<Molecule *> molecules;

};

#endif /* PARTICLE CELL_H_ */
