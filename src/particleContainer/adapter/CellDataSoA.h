/*
 * CellDataSoA.h
 *
 * @Date: 25.03.2013
 * @Author: eckhardw
 */

#ifndef CELLDATASOA_H_
#define CELLDATASOA_H_

#include "CellDataSoABase.h"
#include "utils/AlignedArrayTriplet.h"
#include "utils/ConcatenatedSites.h"
#include "vectorization/SIMD_TYPES.h"
#include <cstdint>
#include <array>

/**
 * \brief Structure of Arrays for vectorized force calculation.
 * \author Johannes Heckl, Wolfgang Eckhardt, Uwe Ehmann
 */
class CellDataSoA : public CellDataSoABase {

	//for better readability:
	typedef ConcatenatedSites<vcp_real_calc>::SiteType 			SiteType;
	typedef ConcatenatedSites<vcp_real_calc>::CoordinateType	CoordinateType;


public:
	CellDataSoA(size_t mol_arg, size_t ljc_arg, size_t charges_arg, size_t dipoles_arg, size_t quadrupoles_arg) {
		resize(mol_arg, ljc_arg, charges_arg, dipoles_arg, quadrupoles_arg);
	}

	/**
	 * \brief
	 */
	enum class QuantityType {
		MOL_POSITION, CENTER_POSITION, FORCE, VIRIAL
	};

	size_t getMolNum() const { return _mol_num; }

	size_t _mol_num;
	size_t _ljc_num;
	size_t _charges_num;
	size_t _dipoles_num;
	size_t _quadrupoles_num;
	size_t _centers_num;

	// entries per molecule
	AlignedArrayTriplet<vcp_real_calc> _mol_pos;
	AlignedArray<int> _mol_ljc_num;
	AlignedArray<int> _mol_charges_num;
	AlignedArray<int> _mol_dipoles_num;
	AlignedArray<int> _mol_quadrupoles_num;

	// entries per center
	ConcatenatedSites<vcp_real_calc> _centers_m_r;
	ConcatenatedSites<vcp_real_calc> _centers_r;
	ConcatenatedSites<vcp_real_calc> _centers_f;
	ConcatenatedSites<vcp_real_calc> _centers_V;

	// entries per lj center
	AlignedArray<vcp_ljc_id_t> _ljc_id;

	// entries per charge
	AlignedArray<vcp_real_calc> _charges_q;

	// entries per dipole
	AlignedArray<vcp_real_calc> _dipoles_p; // dipole moment
	AlignedArrayTriplet<vcp_real_calc> _dipoles_e; // orientation vector of dipole moment
	AlignedArrayTriplet<vcp_real_calc> _dipoles_M; // torque vector

	// entries per quadrupole
	AlignedArray<vcp_real_calc> _quadrupoles_m; // quadrupole moment
	AlignedArrayTriplet<vcp_real_calc> _quadrupoles_e; // orientation vector of quadrupole moment
	AlignedArrayTriplet<vcp_real_calc> _quadrupoles_M; // torque vector


	/**
	 * \brief	Get Pointer to the beginning of the specified data
	 * \details
	 * \tparam	qt	The Quantity you want to access, as defined in enum QuantityType
	 * \tparam	st	The site one wants to access in the quantity (LJC, CHARGE, DIPOLE or QUADRUPOLE)
	 * \tparam	coord	Choose the coordinate of the site you want
	 * \return	Pointer to the first element of the data you requested
	 */
	vcp_inline vcp_real_calc* getBegin(QuantityType qt, SiteType st, CoordinateType coord) {
		ConcatenatedSites<vcp_real_calc> * thisQuantity = resolveQuantity(qt);
		return thisQuantity->getBeginPointer(st, coord);
	}

	vcp_inline const vcp_real_calc* getBegin(QuantityType qt, SiteType st, CoordinateType coord) const {
		const ConcatenatedSites<vcp_real_calc> * thisQuantity = resolveQuantity(qt);
		return thisQuantity->getBeginPointer(st, coord);
	}

	/**
	 * \brief	Get a triplet of data from a ConcatenatedSites at specific index
	 */
	vcp_inline std::array<vcp_real_calc, 3> getTriplet(QuantityType qt, SiteType st, size_t index) const {
		const ConcatenatedSites<vcp_real_calc>* thisQuantity = resolveQuantity(qt);
		return thisQuantity->getTriplet(st, index);
	}

	/**
	 * \brief	Set a triplet of data in a ConcatenatedSites to specified values
	 */
	vcp_inline void setTriplet(std::array<vcp_real_calc, 3> t, QuantityType qt, SiteType st, size_t index) {
		ConcatenatedSites<vcp_real_calc>* thisQuantity = resolveQuantity(qt);
		thisQuantity->setTriplet(t, st, index);
	}

	/**
	 * \brief	Add a set of LJC-data at position index
	 */
	void pushBackLJC(const size_t index, std::array<vcp_real_calc,3> moleculePos, std::array<vcp_real_calc,3> centerPos, vcp_ljc_id_t lookUpIndex) {
		setTriplet(moleculePos, QuantityType::MOL_POSITION, SiteType::LJC, index);
		setTriplet(centerPos, QuantityType::CENTER_POSITION, SiteType::LJC, index);
		_ljc_id[index] = lookUpIndex;
	}

	/**
	 * \brief	Add a set of charge-data at position index
	 */
	void pushBackCharge(const size_t index, std::array<vcp_real_calc,3> moleculePos, std::array<vcp_real_calc,3> centerPos, vcp_real_calc charge) {
		setTriplet(moleculePos, QuantityType::MOL_POSITION, SiteType::CHARGE, index);
		setTriplet(centerPos, QuantityType::CENTER_POSITION, SiteType::CHARGE, index);
		_charges_q[index] = charge;
	}

	/**
	 * \brief	Add a set of dipole-data at position index
	 */
	void pushBackDipole(const size_t index, std::array<vcp_real_calc,3> moleculePos, std::array<vcp_real_calc,3> centerPos,
			vcp_real_calc dipoleMoment, std::array<vcp_real_calc,3> orientation) {
		setTriplet(moleculePos, QuantityType::MOL_POSITION, SiteType::DIPOLE, index);
		setTriplet(centerPos, QuantityType::CENTER_POSITION, SiteType::DIPOLE, index);
		_dipoles_p[index] = dipoleMoment;
		_dipoles_e.x(index) = orientation[0];
		_dipoles_e.y(index) = orientation[1];
		_dipoles_e.z(index) = orientation[2];
	}

	/**
	 * \brief	Add a set of quadrupole-data at position index
	 */
	void pushBackQuadrupole(const size_t index, std::array<vcp_real_calc,3> moleculePos, std::array<vcp_real_calc,3> centerPos,
			vcp_real_calc quadrupoleMoment, std::array<vcp_real_calc,3> orientation) {
		setTriplet(moleculePos, QuantityType::MOL_POSITION, SiteType::QUADRUPOLE, index);
		setTriplet(centerPos, QuantityType::CENTER_POSITION, SiteType::QUADRUPOLE, index);
		_quadrupoles_m[index] = quadrupoleMoment;
		_quadrupoles_e.x(index) = orientation[0];
		_quadrupoles_e.y(index) = orientation[1];
		_quadrupoles_e.z(index) = orientation[2];
	}

	void vcp_inline initDistLookupPointers(
			AlignedArray<vcp_lookupOrMask_single>& centers_dist_lookup,
			vcp_lookupOrMask_single*& ljc_dist_lookup,
			vcp_lookupOrMask_single*& charges_dist_lookup,
			vcp_lookupOrMask_single*& dipoles_dist_lookup,
			vcp_lookupOrMask_single*& quadrupoles_dist_lookup) const {

		size_t ljc_size 	= AlignedArray<vcp_real_calc>::_round_up(_ljc_num);
		size_t charges_size = AlignedArray<vcp_real_calc>::_round_up(_charges_num);
		size_t dipoles_size = AlignedArray<vcp_real_calc>::_round_up(_dipoles_num);
		size_t quadrupoles_size = AlignedArray<vcp_real_calc>::_round_up(_quadrupoles_num);
		size_t centers_size = ljc_size + charges_size + dipoles_size + quadrupoles_size;

		centers_dist_lookup.resize_zero_shrink(centers_size);
		setPaddingToZero(centers_dist_lookup);

		ljc_dist_lookup = centers_dist_lookup;
		charges_dist_lookup = ljc_dist_lookup + (ljc_size + VCP_INDICES_PER_LOOKUP_SINGLE_M1)/VCP_INDICES_PER_LOOKUP_SINGLE;
		dipoles_dist_lookup = charges_dist_lookup + (charges_size + VCP_INDICES_PER_LOOKUP_SINGLE_M1)/VCP_INDICES_PER_LOOKUP_SINGLE;
		quadrupoles_dist_lookup = dipoles_dist_lookup + (dipoles_size + VCP_INDICES_PER_LOOKUP_SINGLE_M1)/VCP_INDICES_PER_LOOKUP_SINGLE;
	}

	void vcp_inline initDistLookupPointersSingle(
			AlignedArray<vcp_lookupOrMask_single>& centers_dist_lookup,
			vcp_lookupOrMask_single*& sites_dist_lookup,
			size_t sites_num) const {

		centers_dist_lookup.resize_zero_shrink(sites_num, true, false);
		sites_dist_lookup = centers_dist_lookup;
	}

	void resize(size_t molecules_arg, size_t ljcenters_arg, size_t charges_arg, size_t dipoles_arg, size_t quadrupoles_arg) {
//		const bool allow_shrink = false; // TODO shrink at some point in the future

		_mol_num = molecules_arg;
		_ljc_num = ljcenters_arg;
		_charges_num = charges_arg;
		_dipoles_num = dipoles_arg;
		_quadrupoles_num = quadrupoles_arg;

		// entries per molecule
		_mol_pos			.resize_zero_shrink(_mol_num);
		_mol_ljc_num		.resize_zero_shrink(_mol_num);
		_mol_charges_num	.resize_zero_shrink(_mol_num);
		_mol_dipoles_num	.resize_zero_shrink(_mol_num);
		_mol_quadrupoles_num.resize_zero_shrink(_mol_num);

		_centers_m_r.resize(ljcenters_arg, charges_arg, dipoles_arg, quadrupoles_arg);
		_centers_r	.resize(ljcenters_arg, charges_arg, dipoles_arg, quadrupoles_arg);
		_centers_f	.resize(ljcenters_arg, charges_arg, dipoles_arg, quadrupoles_arg);
		_centers_V	.resize(ljcenters_arg, charges_arg, dipoles_arg, quadrupoles_arg);

		// entries per lj center
		_ljc_id.resize_zero_shrink(_ljc_num, true);

		// entries per charge
		_charges_q.resize_zero_shrink(_charges_num);

		// entries per dipole
		_dipoles_p.resize_zero_shrink(_dipoles_num);
		_dipoles_e.resize_zero_shrink(_dipoles_num);
		_dipoles_M.resize_zero_shrink(_dipoles_num);

		// entries per quadrupole
		_quadrupoles_m.resize_zero_shrink(_quadrupoles_num);
		_quadrupoles_e.resize_zero_shrink(_quadrupoles_num);
		_quadrupoles_M.resize_zero_shrink(_quadrupoles_num);

	}

	size_t getDynamicSize() const {
		size_t total = 0;

		total += _mol_pos.get_dynamic_memory();
		total += _mol_ljc_num.get_dynamic_memory();
		total += _mol_charges_num.get_dynamic_memory();
		total += _mol_dipoles_num.get_dynamic_memory();
		total += _mol_quadrupoles_num.get_dynamic_memory();

		total += _centers_m_r.get_dynamic_memory();
		total += _centers_r.get_dynamic_memory();
		total += _centers_f.get_dynamic_memory();
		total += _centers_V.get_dynamic_memory();

		total += _dipoles_p.get_dynamic_memory();
		total += _dipoles_e.get_dynamic_memory();
		total += _dipoles_M.get_dynamic_memory();

		total += _quadrupoles_m.get_dynamic_memory();
		total += _quadrupoles_e.get_dynamic_memory();
		total += _quadrupoles_M.get_dynamic_memory();

		return total;
	}

private:

	template <typename T>
	vcp_inline void setPaddingToZero(AlignedArray<T>& t) const {
		size_t ljc_size = t._round_up(_ljc_num);
		size_t charges_size = t._round_up(_charges_num);
		size_t dipoles_size = t._round_up(_dipoles_num);

		t.zero(_ljc_num);
		t.zero(ljc_size + _charges_num);
		t.zero(ljc_size + charges_size + _dipoles_num);
		t.zero(ljc_size + charges_size + dipoles_size + _quadrupoles_num);
	}

	/**
	 * \brief Matches the given QuantityType and returns a pointer to the associated ConcatenatedSites
	 */
	vcp_inline ConcatenatedSites<vcp_real_calc>* resolveQuantity(QuantityType qt) {
		ConcatenatedSites<vcp_real_calc>* returnQuantity;
		switch(qt) {
		case QuantityType::MOL_POSITION:
			returnQuantity = &_centers_m_r;
			break;
		case QuantityType::CENTER_POSITION:
			returnQuantity = &_centers_r;
			break;
		case QuantityType::FORCE:
			returnQuantity = &_centers_f;
			break;
		case QuantityType::VIRIAL:
			returnQuantity = &_centers_V;
			break;
		default:
			returnQuantity = nullptr;
			break;
		}

		mardyn_assert(returnQuantity != nullptr);
		return returnQuantity;
	}

	vcp_inline const ConcatenatedSites<vcp_real_calc>* resolveQuantity(QuantityType qt) const {
		const ConcatenatedSites<vcp_real_calc>* returnQuantity;
		switch(qt) {
		case QuantityType::MOL_POSITION:
			returnQuantity = &_centers_m_r;
			break;
		case QuantityType::CENTER_POSITION:
			returnQuantity = &_centers_r;
			break;
		case QuantityType::FORCE:
			returnQuantity = &_centers_f;
			break;
		case QuantityType::VIRIAL:
			returnQuantity = &_centers_V;
			break;
		default:
			returnQuantity = nullptr;
			break;
		}

		mardyn_assert(returnQuantity != nullptr);
		return returnQuantity;
	}

};

#endif /* CELLDATASOA_H_ */
