/*
 * CellDataSoAWR.h
 *
 * @Date: 11.01.2017
 * @Author: tchipevn
 */

#ifndef CELLDATASOA_WR_H_
#define CELLDATASOA_WR_H_

#include "utils/AlignedArrayTriplet.h"
#include "vectorization/SIMD_TYPES.h"
#include "molecules/Molecule.h"
#include "CellDataSoABase.h"
#include <cstdint>

/**
 * \brief Structure of Arrays for single-center lennard-Jones molecules for
 * the WR run.
 * \author Nikola Tchipev
 */
class CellDataSoA_WR : public CellDataSoABase {
public:
	CellDataSoA_WR(size_t mol_arg) {
		resize(mol_arg);
	}

	size_t _mol_num;

	// entries per molecule
	AlignedArrayTriplet<vcp_real_calc> _mol_r;
	AlignedArrayTriplet<vcp_real_calc> _mol_v;
	AlignedArray<uint64_t> _mol_uid;

	vcp_inline vcp_real_calc* r_xBegin()   const { return _mol_r.xBegin();}
	vcp_inline vcp_real_calc* r_yBegin()   const { return _mol_r.yBegin();}
	vcp_inline vcp_real_calc* r_zBegin()   const { return _mol_r.zBegin();}
	vcp_inline vcp_real_calc* v_xBegin()   const { return _mol_v.xBegin();}
	vcp_inline vcp_real_calc* v_yBegin()   const { return _mol_v.yBegin();}
	vcp_inline vcp_real_calc* v_zBegin()   const { return _mol_v.zBegin();}

	void resize(size_t molecules_arg) {
		const bool allow_shrink = false; // TODO shrink at some point in the future

		_mol_num = molecules_arg;

		// entries per molecule
		_mol_r.resize_zero_shrink(_mol_num);
		_mol_v.resize_zero_shrink(_mol_num);
		_mol_uid.resize(_mol_num);
	}

	size_t getDynamicSize() const {
		size_t total = 0;

		total += _mol_r.get_dynamic_memory();
		total += _mol_v.get_dynamic_memory();
		total += _mol_uid.get_dynamic_memory();

		return total;
	}

	void appendMolecule(Molecule_WR& m) {
		_mol_r.appendValueTriplet(m.r(0), m.r(1), m.r(2), _mol_num);
		_mol_v.appendValueTriplet(m.v(0), m.v(1), m.v(2), _mol_num);
		_mol_uid.appendValue(m.id(), _mol_num);
		++_mol_num;
	}

	void readImmutableMolecule(size_t index, Molecule_WR& m) const {
		// changes in AOS storage will not be saved
		m.setStorageState(Molecule_WR::STORAGE_AOS);
		m.setr(0, _mol_r.x(index));
		m.setr(1, _mol_r.y(index));
		m.setr(2, _mol_r.z(index));
		m.setv(0, _mol_v.x(index));
		m.setv(1, _mol_v.y(index));
		m.setv(2, _mol_v.z(index));
		m.setid(_mol_uid[index]);
	}

	void readMutableMolecule(size_t index, Molecule_WR& m) {
		// changes in SOA storage will be saved
		m.setStorageState(Molecule_WR::STORAGE_SOA);
		m.setSoA(this);
		m.setStartIndexSoA_LJ(index);
	}
};

#endif /* CELLDATASOA_WR_H_ */
