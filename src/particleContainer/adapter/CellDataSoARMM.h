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
	typedef ConcatenatedAlignedArrayRMM<vcp_real_calc, vcp_real_accum, uint64_t>::Quantity_t Quantity_t;
public:
	CellDataSoARMM(size_t mol_arg) {
		resize(mol_arg);
	}

	vcp_inline vcp_real_calc* r_xBegin() { return _data.begin_calc(Quantity_t::RX);}
	vcp_inline vcp_real_calc* r_yBegin() { return _data.begin_calc(Quantity_t::RY);}
	vcp_inline vcp_real_calc* r_zBegin() { return _data.begin_calc(Quantity_t::RZ);}
	vcp_inline vcp_real_accum* v_xBegin() { return _data.begin_accum(Quantity_t::VX);}
	vcp_inline vcp_real_accum* v_yBegin() { return _data.begin_accum(Quantity_t::VY);}
	vcp_inline vcp_real_accum* v_zBegin() { return _data.begin_accum(Quantity_t::VZ);}

	const vcp_inline vcp_real_calc* r_xBegin() const { return _data.begin_calc(Quantity_t::RX);}
	const vcp_inline vcp_real_calc* r_yBegin() const { return _data.begin_calc(Quantity_t::RY);}
	const vcp_inline vcp_real_calc* r_zBegin() const { return _data.begin_calc(Quantity_t::RZ);}
	const vcp_inline vcp_real_accum* v_xBegin() const { return _data.begin_accum(Quantity_t::VX);}
	const vcp_inline vcp_real_accum* v_yBegin() const { return _data.begin_accum(Quantity_t::VY);}
	const vcp_inline vcp_real_accum* v_zBegin() const { return _data.begin_accum(Quantity_t::VZ);}

	void resize(size_t molecules_arg) {
		const bool allow_shrink = false; // TODO shrink at some point in the future

		setMolNum(molecules_arg);

		// entries per molecule
		_data.resize(getMolNum());
	}

	size_t getDynamicSize() const {
		return _data.get_dynamic_memory();
	}

	void appendMolecule(MoleculeInterface& m) {
		MoleculeRMM& m_RMM = downcastMoleculeReferenceRMM(m);
		std::array<vcp_real_calc, 3> calcs = {
			static_cast<vcp_real_calc>(m_RMM.r(0)),
			static_cast<vcp_real_calc>(m_RMM.r(1)),
			static_cast<vcp_real_calc>(m_RMM.r(2))
		};
		std::array<vcp_real_accum, 3> accums = {
			static_cast<vcp_real_accum>(m_RMM.v(0)),
			static_cast<vcp_real_accum>(m_RMM.v(1)),
			static_cast<vcp_real_accum>(m_RMM.v(2))
		};

		_data.appendValues(calcs, accums, m_RMM.id(), getMolNum());
		incrementMolNum();
	}

	void increaseStorage(size_t additionalMolecules) {
		_data.increaseStorage(getMolNum(), additionalMolecules);
	}

	Molecule buildAoSMolecule(size_t index) const {
		return Molecule (
			getMolUid(index), nullptr,
			getMolR(0,index), getMolR(1,index), getMolR(2,index),
			getMolV(0,index), getMolV(1,index), getMolV(2,index));
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
		mardyn_assert(index < getMolNum());
		if(getMolNum() > 1 and index < getMolNum() - 1) {
			setMolR(0, index, getMolR(0,getMolNum()-1));
			setMolR(1, index, getMolR(1,getMolNum()-1));
			setMolR(2, index, getMolR(2,getMolNum()-1));
			setMolV(0, index, getMolV(0,getMolNum()-1));
			setMolV(1, index, getMolV(1,getMolNum()-1));
			setMolV(2, index, getMolV(2,getMolNum()-1));
			setMolUid(index, getMolUid(getMolNum()-1));
		}
		decrementMolNum();
	}

	void prefetchForForce() const {
		_data.prefetchForForce();
	}

	vcp_real_calc getMolR(unsigned short d, size_t index) const {
		mardyn_assert(d < 3);
		Quantity_t q = static_cast<Quantity_t>(d);
		return _data.get_calc(q,index);
	}

	void setMolR(unsigned short d, size_t index, vcp_real_calc molR) {
		mardyn_assert(d < 3);
		Quantity_t q = static_cast<Quantity_t>(d);
		_data.get_calc(q,index) = molR;
	}

	vcp_real_accum getMolV(unsigned short d, size_t index) const {
		mardyn_assert(d < 3);
		Quantity_t q = static_cast<Quantity_t>(d + 3); // +3 to convert to VX, VY, VZ
		return _data.get_accum(q,index);
	}

	void setMolV(unsigned short d, size_t index, vcp_real_accum molV) {
		mardyn_assert(d < 3);
		Quantity_t q = static_cast<Quantity_t>(d + 3); // +3 to convert to VX, VY, VZ
		_data.get_accum(q,index) = molV;
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
	// entries per molecule
	ConcatenatedAlignedArrayRMM<vcp_real_calc, vcp_real_accum, uint64_t> _data;
};

#endif /* CELLDATASOARMM_H_ */

