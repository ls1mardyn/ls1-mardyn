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
#include "molecules/Molecule_WR.h"
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

	vcp_inline vcp_real_calc* r_xBegin() { return _mol_r.xBegin();}
	vcp_inline vcp_real_calc* r_yBegin() { return _mol_r.yBegin();}
	vcp_inline vcp_real_calc* r_zBegin() { return _mol_r.zBegin();}
	vcp_inline vcp_real_calc* v_xBegin() { return _mol_v.xBegin();}
	vcp_inline vcp_real_calc* v_yBegin() { return _mol_v.yBegin();}
	vcp_inline vcp_real_calc* v_zBegin() { return _mol_v.zBegin();}

	const vcp_inline vcp_real_calc* r_xBegin()   const { return _mol_r.xBegin();}
	const vcp_inline vcp_real_calc* r_yBegin()   const { return _mol_r.yBegin();}
	const vcp_inline vcp_real_calc* r_zBegin()   const { return _mol_r.zBegin();}
	const vcp_inline vcp_real_calc* v_xBegin()   const { return _mol_v.xBegin();}
	const vcp_inline vcp_real_calc* v_yBegin()   const { return _mol_v.yBegin();}
	const vcp_inline vcp_real_calc* v_zBegin()   const { return _mol_v.zBegin();}

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

	void appendMolecule(MoleculeInterface& m) {
		Molecule_WR& m_wr = downcastMoleculeReferenceWR(m);

		_mol_r.appendValueTriplet(m_wr.r(0), m_wr.r(1), m_wr.r(2), _mol_num);
		_mol_v.appendValueTriplet(m_wr.v(0), m_wr.v(1), m_wr.v(2), _mol_num);
		_mol_uid.appendValue(m_wr.id(), _mol_num);
		++_mol_num;
	}

	void increaseStorage(size_t additionalMolecules) {
		_mol_r.increaseStorage(_mol_num, additionalMolecules);
		_mol_v.increaseStorage(_mol_num, additionalMolecules);
		_mol_uid.increaseStorage(_mol_num, additionalMolecules);
	}

	void readImmutableMolecule(size_t index, MoleculeInterface& m) const {
		Molecule_WR& m_wr = downcastMoleculeReferenceWR(m);

		// changes in AOS storage will not be saved
		m_wr.setStorageState(Molecule_WR::STORAGE_AOS);
		m_wr.setr(0, _mol_r.x(index));
		m_wr.setr(1, _mol_r.y(index));
		m_wr.setr(2, _mol_r.z(index));
		m_wr.setv(0, _mol_v.x(index));
		m_wr.setv(1, _mol_v.y(index));
		m_wr.setv(2, _mol_v.z(index));
		m_wr.setid(_mol_uid[index]);
	}

	void readMutableMolecule(size_t index, MoleculeInterface& m) {
		Molecule_WR& m_wr = downcastMoleculeReferenceWR(m);

		// changes in SOA storage will be saved
		m_wr.setStorageState(Molecule_WR::STORAGE_SOA);
		m_wr.setSoA(this);
		m_wr.setStartIndexSoA_LJ(index);
	}

	void writeMolecule(size_t i, const MoleculeInterface& m) {
//		Molecule_WR& m_wr = downcastReferenceWR(m);

		_mol_r.x(i) = m.r(0);
		_mol_r.y(i) = m.r(1);
		_mol_r.z(i) = m.r(2);
		_mol_v.x(i) = m.v(0);
		_mol_v.y(i) = m.v(1);
		_mol_v.z(i) = m.v(2);
		_mol_uid[i] = m.id();
	}

	void deleteMolecule(size_t index) {
		mardyn_assert(index < static_cast<size_t>(_mol_num));
		if(_mol_num > 1 and index < _mol_num - 1) {
			_mol_r.x(index) = _mol_r.x(_mol_num-1);
			_mol_r.y(index) = _mol_r.y(_mol_num-1);
			_mol_r.z(index) = _mol_r.z(_mol_num-1);
			_mol_v.x(index) = _mol_v.x(_mol_num-1);
			_mol_v.y(index) = _mol_v.y(_mol_num-1);
			_mol_v.z(index) = _mol_v.z(_mol_num-1);
			_mol_uid[index] = _mol_uid[_mol_num-1];
		}
		--_mol_num;
	}

	void prefetch() const {
		_mol_r.prefetch();
		_mol_v.prefetch();
	}

	size_t getMolNum() const {
		return _mol_num;
	}

	void setMolNum(size_t molNum) {
		_mol_num = molNum;
	}

	vcp_real_calc getMolR(unsigned short d, size_t index) const {
		mardyn_assert(d < 3);
		switch(d) {
		case 0:
			return _mol_r.x(index);
		case 1:
			return _mol_r.y(index);
		case 2:
			return _mol_r.z(index);
		default:
			mardyn_assert(false);
			return static_cast<vcp_real_calc>(0.0);
		}
	}

	void setMolR(unsigned short d, size_t index, vcp_real_calc molR) {
		mardyn_assert(d < 3);
		switch(d) {
		case 0:
			_mol_r.x(index) = molR;
			break;
		case 1:
			_mol_r.y(index) = molR;
			break;
		case 2:
			_mol_r.z(index) = molR;
			break;
		default:
			mardyn_assert(false);
		}
	}

	vcp_real_calc getMolV(unsigned short d, size_t index) const {
		mardyn_assert(d < 3);
		switch(d) {
		case 0:
			return _mol_v.x(index);
		case 1:
			return _mol_v.y(index);
		case 2:
			return _mol_v.z(index);
		default:
			mardyn_assert(false);
			return static_cast<vcp_real_calc>(0.0);
		}
	}

	void setMolV(unsigned short d, size_t index, vcp_real_calc molV) {
		mardyn_assert(d < 3);
		switch(d) {
		case 0:
			_mol_v.x(index) = molV;
			break;
		case 1:
			_mol_v.y(index) = molV;
			break;
		case 2:
			_mol_v.z(index) = molV;
			break;
		default:
			mardyn_assert(false);
		}
	}

	uint64_t getMolUid(size_t index) const {
		return _mol_uid[index];
	}

	void setMolUid(size_t index, unsigned long molUid) {
		_mol_uid[index] = molUid;
	}

private:
	size_t _mol_num;

	// entries per molecule
	AlignedArrayTriplet<vcp_real_calc> _mol_r;
	AlignedArrayTriplet<vcp_real_calc> _mol_v;
	AlignedArray<uint64_t> _mol_uid;
};

#endif /* CELLDATASOA_WR_H_ */
