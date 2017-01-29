#ifndef PARTICLE_CELL_H_
#define PARTICLE_CELL_H_

#include <vector>

#include "Cell.h"
#include "particleContainer/ParticleCellBase.h"
#include "particleContainer/adapter/CellDataSoA.h"

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
/**
 * \details <br>(Johannes Heckl)<br>
 * Also stores data for various CellProcessor%s.<br>
 * If you add a new data member, update the _assign() method with deep copy<br>
 * semantics for the new data member.<br>
 * Uses the default copy constructor and the default assignment operator despite<br>
 * having pointer data members. This is because these data members are not controlled<br>
 * by the ParticleCell itself, but by the various CellProcessor%s so ParticleCell can not<br>
 * know the proper copy semantics. This should not cause any problems because no copy<br>
 * actions should be executed during CellProcessor applications.
 */
class ParticleCell : public ParticleCellBase {
/*private:
	ParticleCell(const ParticleCell& that);*/
public:
	/**
	 * \brief Initialize data pointers to 0.
	 * \author Johannes Heckl
	 */
	ParticleCell() ;
	/**
	 * \brief Destructor.
	 * \author Johannes Heckl
	 */
	~ParticleCell() ;

	//! removes and deallocates all elements
	void deallocateAllParticles();

	//! insert a single molecule into this cell
	bool addParticle(Molecule& particle, bool checkWhetherDuplicate = false);

	Molecule& moleculesAt(size_t i) {
		return _molecules.at(i);
	}

	const Molecule& moleculesAtConst(size_t i) const {
		return _molecules.at(i);
	}

	bool isEmpty() const;

	bool deleteMoleculeByIndex(size_t index);

	//! return the number of molecules contained in this cell
	int getMoleculeCount() const;
	
	/**
	 * \brief Get the structure of arrays for VectorizedCellProcessor.
	 * \author Johannes Heckl
	 */
	CellDataSoA& getCellDataSoA() {
		return _cellDataSoA;
	}

	bool testInBox(const Molecule& particle) const{
		return particle.inBox(_boxMin, _boxMax);
	}
	/**
	 * filter molecules which have left the box
	 * @return field vector containing leaving molecules
	 */
	//std::vector<Molecule> & filterLeavingMolecules();

	void preUpdateLeavingMolecules();

	void updateLeavingMoleculesBase(ParticleCellBase& otherCell);

	void updateLeavingMolecules(ParticleCell& otherCell);

	void postUpdateLeavingMolecules();

	void getRegion(double lowCorner[3], double highCorner[3], std::vector<Molecule*> &particlePtrs, bool removeFromContainer = false);

    void buildSoACaches();

    void reserveMoleculeStorage(size_t numMols);


private:
	/**
	 * \brief A vector of pointers to the Molecules in this cell.
	 */
	std::vector<Molecule> _molecules;

	/**
	 * \brief A vector of molecules, which have left this cell.
	 */
	std::vector<Molecule> _leavingMolecules;

	/**
	 * \brief Structure of arrays for VectorizedCellProcessor.
	 * \author Johannes Heckl
	 */
	CellDataSoA _cellDataSoA;
};

#endif /* PARTICLE CELL_H_ */
