/***********************************************************************************//**
 *
 * \file ParticleIterator.h
 *
 * \brief Iterator class for cell based particle containers
 *
 * \author Wolfgang Hölzl (hoelzlw), hoelzlw AT in.tum.de
 *
 * \details This class implements a forward iterator for the cell based particle containers.
 * It does a coarse iteration over the cells and, within each cell, a fine iteration over the particles.
 *
 **************************************************************************************/

#ifndef  ParticleIterator_INC
#define  ParticleIterator_INC

#include <vector>
#include <stdexcept>
#include "utils/mardyn_assert.h"
#include "ParticleCellBase.h"

class ParticleContainer;

class ParticleIterator {
public:
	enum Type {
		ALL_CELLS=0, /* iterates every cell */
		ONLY_INNER_AND_BOUNDARY=1, /* iterates only inner and boundary cells, i.e. no halo cells */
	};

	typedef ParticleContainer CellContainer_T;
	typedef CellContainer_T* CellContainer_T_ptr;
	typedef size_t CellIndex_T;
	typedef size_t MolIndex_T;

	ParticleIterator ();
	ParticleIterator (Type t_arg, CellContainer_T_ptr cells_arg, const CellIndex_T offset_arg, const CellIndex_T stride_arg, const bool initialize=true);
	ParticleIterator& operator=(const ParticleIterator& other);

	virtual ~ParticleIterator(){}

	void operator ++ ();

	bool operator == (const ParticleIterator& other) const;
	bool operator != (const ParticleIterator& other) const;

	Molecule& operator *  () const;
	Molecule* operator -> () const;

	static ParticleIterator invalid ();

	void deleteCurrentParticle();

	CellIndex_T getCellIndex(){return _cell_index;}

protected:
	Type _type;
	CellContainer_T_ptr _cells;

	CellIndex_T _cell_index;
	MolIndex_T _mol_index;

	bool _currentParticleDeleted;

	const CellIndex_T _stride;

	virtual void make_invalid ();
	virtual void next_non_empty_cell();
};


#endif /* #ifndef ParticleIterator_INC */
