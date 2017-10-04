#ifndef CELLDATASOARMM_H_
#define CELLDATASOARMM_H_

#include "utils/ConcatenatedAlignedArrayRMM.h"
#include "vectorization/SIMD_TYPES.h"
#include "molecules/Molecule.h"
#include "molecules/MoleculeRMM.h"
#include "CellDataSoABase.h"
#include <cstdint>

/**
 * \brief Structure of Arrays for single-center lennard-Jones molecules for
 * the RMM run.
 * \author Nikola Tchipev
 */
class CellDataSoARMM : public CellDataSoABase {
	typedef ConcatenatedAlignedArrayRMM<vcp_real_calc, uint64_t>::Quantity_t Quantity_t;
public:
	CellDataSoARMM(size_t mol_arg) {
		resize(mol_arg);
	}

	vcp_inline vcp_real_calc* r_xBegin() { return _data.begin_real(Quantity_t::RX);}
	vcp_inline vcp_real_calc* r_yBegin() { return _data.begin_real(Quantity_t::RY);}
	vcp_inline vcp_real_calc* r_zBegin() { return _data.begin_real(Quantity_t::RZ);}
	vcp_inline vcp_real_calc* v_xBegin() { return _data.begin_real(Quantity_t::VX);}
	vcp_inline vcp_real_calc* v_yBegin() { return _data.begin_real(Quantity_t::VY);}
	vcp_inline vcp_real_calc* v_zBegin() { return _data.begin_real(Quantity_t::VZ);}

	const vcp_inline vcp_real_calc* r_xBegin() const { return _data.begin_real(Quantity_t::RX);}
	const vcp_inline vcp_real_calc* r_yBegin() const { return _data.begin_real(Quantity_t::RY);}
	const vcp_inline vcp_real_calc* r_zBegin() const { return _data.begin_real(Quantity_t::RZ);}
	const vcp_inline vcp_real_calc* v_xBegin() const { return _data.begin_real(Quantity_t::VX);}
	const vcp_inline vcp_real_calc* v_yBegin() const { return _data.begin_real(Quantity_t::VY);}
	const vcp_inline vcp_real_calc* v_zBegin() const { return _data.begin_real(Quantity_t::VZ);}

	void resize(size_t molecules_arg) {
		const bool allow_shrink = false; // TODO shrink at some point in the future

		_mol_num = molecules_arg;

		// entries per molecule
		_data.resize(_mol_num);
	}

	size_t getDynamicSize() const {
		return _data.get_dynamic_memory();
	}

	void appendMolecule(MoleculeInterface& m) {
		MoleculeRMM& m_RMM = downcastMoleculeReferenceRMM(m);
		std::array<vcp_real_calc,6> vals = {
			static_cast<vcp_real_calc>(m_RMM.r(0)),
			static_cast<vcp_real_calc>(m_RMM.r(1)),
			static_cast<vcp_real_calc>(m_RMM.r(2)),
			static_cast<vcp_real_calc>(m_RMM.v(0)),
			static_cast<vcp_real_calc>(m_RMM.v(1)),
			static_cast<vcp_real_calc>(m_RMM.v(2))
		};

		_data.appendValues(vals, m_RMM.id(), _mol_num);
		++_mol_num;
	}

	void increaseStorage(size_t additionalMolecules) {
		_data.increaseStorage(_mol_num, additionalMolecules);
	}

	void readImmutableMolecule(size_t index, MoleculeInterface& m) const {
		MoleculeRMM& m_RMM = downcastMoleculeReferenceRMM(m);

		// changes in AOS storage will not be saved
		m_RMM.setStorageState(MoleculeRMM::STORAGE_AOS);
		m_RMM.setr(0, getMolR(0,index));
		m_RMM.setr(1, getMolR(1,index));
		m_RMM.setr(2, getMolR(2,index));
		m_RMM.setv(0, getMolV(0,index));
		m_RMM.setv(1, getMolV(1,index));
		m_RMM.setv(2, getMolV(2,index));
		m_RMM.setid(getMolUid(index));
	}

	void readMutableMolecule(size_t index, MoleculeInterface& m) {
		MoleculeRMM& m_RMM = downcastMoleculeReferenceRMM(m);

		// changes in SOA storage will be saved
		m_RMM.setStorageState(MoleculeRMM::STORAGE_SOA);
		m_RMM.setSoA(this);
		m_RMM.setStartIndexSoA_LJ(index);
	}

	void writeMolecule(size_t i, const MoleculeInterface& m) {
		setMolR(0, i, static_cast<vcp_real_calc>(m.r(0)));
		setMolR(1, i, static_cast<vcp_real_calc>(m.r(1)));
		setMolR(2, i, static_cast<vcp_real_calc>(m.r(2)));
		setMolV(0, i, static_cast<vcp_real_calc>(m.v(0)));
		setMolV(1, i, static_cast<vcp_real_calc>(m.v(1)));
		setMolV(2, i, static_cast<vcp_real_calc>(m.v(2)));
		setMolUid(i, m.id());
	}

	void deleteMolecule(size_t index) {
		mardyn_assert(index < _mol_num);
		if(_mol_num > 1 and index < _mol_num - 1) {
			setMolR(0, index, getMolR(0,_mol_num-1));
			setMolR(1, index, getMolR(1,_mol_num-1));
			setMolR(2, index, getMolR(2,_mol_num-1));
			setMolV(0, index, getMolV(0,_mol_num-1));
			setMolV(1, index, getMolV(1,_mol_num-1));
			setMolV(2, index, getMolV(2,_mol_num-1));
			setMolUid(index, getMolUid(_mol_num-1));
		}
		--_mol_num;
	}

	void prefetchForForce() const {
		_data.prefetchForForce();
	}

	size_t getMolNum() const {
		return _mol_num;
	}

	void setMolNum(size_t molNum) {
		_mol_num = molNum;
	}

	vcp_real_calc getMolR(unsigned short d, size_t index) const {
		mardyn_assert(d < 3);
		Quantity_t q = static_cast<Quantity_t>(d);
		return _data.get_real(q,index);
	}

	void setMolR(unsigned short d, size_t index, vcp_real_calc molR) {
		mardyn_assert(d < 3);
		Quantity_t q = static_cast<Quantity_t>(d);
		_data.get_real(q,index) = molR;
	}

	vcp_real_calc getMolV(unsigned short d, size_t index) const {
		mardyn_assert(d < 3);
		Quantity_t q = static_cast<Quantity_t>(d + 3); // +3 to convert to VX, VY, VZ
		return _data.get_real(q,index);
	}

	void setMolV(unsigned short d, size_t index, vcp_real_calc molV) {
		mardyn_assert(d < 3);
		Quantity_t q = static_cast<Quantity_t>(d + 3); // +3 to convert to VX, VY, VZ
		_data.get_real(q,index) = molV;
	}

	uint64_t getMolUid(size_t index) const {
		Quantity_t q = Quantity_t::UID;
		return _data.get_uid(q, index);
	}

	void setMolUid(size_t index, unsigned long molUid) {
		Quantity_t q = Quantity_t::UID;
		_data.get_uid(q, index) = molUid;
	}

private:
	size_t _mol_num;

	// entries per molecule
	ConcatenatedAlignedArrayRMM<vcp_real_calc, uint64_t> _data;
};

#endif /* CELLDATASOARMM_H_ */

