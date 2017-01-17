/*
 * CellDataSoA.h
 *
 * @Date: 25.03.2013
 * @Author: eckhardw
 */

#ifndef CELLDATASOA_H_
#define CELLDATASOA_H_

#include "utils/AlignedArrayTriplet.h"
#include "vectorization/SIMD_TYPES.h"
#include <cstdint>

/**
 * \brief Structure of Arrays for vectorized force calculation.
 * \author Johannes Heckl, Wolfgang Eckhardt, Uwe Ehmann
 */
class CellDataSoA {
public:
	CellDataSoA(size_t mol_arg, size_t ljc_arg, size_t charges_arg, size_t dipoles_arg, size_t quadrupoles_arg) {
		resize(mol_arg, ljc_arg, charges_arg, dipoles_arg, quadrupoles_arg);
	}

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
	AlignedArrayTriplet<vcp_real_calc> _centers_m_r;
	AlignedArrayTriplet<vcp_real_calc> _centers_r;
	AlignedArrayTriplet<vcp_real_calc> _centers_f;
	AlignedArrayTriplet<vcp_real_calc> _centers_V;

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

	vcp_inline vcp_real_calc* ljc_m_r_xBegin() const { return _centers_m_r.xBegin();}
	vcp_inline vcp_real_calc* ljc_m_r_yBegin() const { return _centers_m_r.yBegin();}
	vcp_inline vcp_real_calc* ljc_m_r_zBegin() const { return _centers_m_r.zBegin();}
	vcp_inline vcp_real_calc* ljc_r_xBegin()   const { return _centers_r  .xBegin();}
	vcp_inline vcp_real_calc* ljc_r_yBegin()   const { return _centers_r  .yBegin();}
	vcp_inline vcp_real_calc* ljc_r_zBegin()   const { return _centers_r  .zBegin();}
	vcp_inline vcp_real_calc* ljc_f_xBegin()   const { return _centers_f  .xBegin();}
	vcp_inline vcp_real_calc* ljc_f_yBegin()   const { return _centers_f  .yBegin();}
	vcp_inline vcp_real_calc* ljc_f_zBegin()   const { return _centers_f  .zBegin();}
	vcp_inline vcp_real_calc* ljc_V_xBegin()   const { return _centers_V  .xBegin();}
	vcp_inline vcp_real_calc* ljc_V_yBegin()   const { return _centers_V  .yBegin();}
	vcp_inline vcp_real_calc* ljc_V_zBegin()   const { return _centers_V  .zBegin();}

	vcp_inline vcp_real_calc* charges_m_r_xBegin() const { return ljc_m_r_xBegin() + _centers_m_r._round_up(_ljc_num);}
	vcp_inline vcp_real_calc* charges_m_r_yBegin() const { return ljc_m_r_yBegin() + _centers_m_r._round_up(_ljc_num);}
	vcp_inline vcp_real_calc* charges_m_r_zBegin() const { return ljc_m_r_zBegin() + _centers_m_r._round_up(_ljc_num);}
	vcp_inline vcp_real_calc* charges_r_xBegin()   const { return ljc_r_xBegin()   + _centers_r  ._round_up(_ljc_num);}
	vcp_inline vcp_real_calc* charges_r_yBegin()   const { return ljc_r_yBegin()   + _centers_r  ._round_up(_ljc_num);}
	vcp_inline vcp_real_calc* charges_r_zBegin()   const { return ljc_r_zBegin()   + _centers_r  ._round_up(_ljc_num);}
	vcp_inline vcp_real_calc* charges_f_xBegin()   const { return ljc_f_xBegin()   + _centers_f  ._round_up(_ljc_num);}
	vcp_inline vcp_real_calc* charges_f_yBegin()   const { return ljc_f_yBegin()   + _centers_f  ._round_up(_ljc_num);}
	vcp_inline vcp_real_calc* charges_f_zBegin()   const { return ljc_f_zBegin()   + _centers_f  ._round_up(_ljc_num);}
	vcp_inline vcp_real_calc* charges_V_xBegin()   const { return ljc_V_xBegin()   + _centers_V  ._round_up(_ljc_num);}
	vcp_inline vcp_real_calc* charges_V_yBegin()   const { return ljc_V_yBegin()   + _centers_V  ._round_up(_ljc_num);}
	vcp_inline vcp_real_calc* charges_V_zBegin()   const { return ljc_V_zBegin()   + _centers_V  ._round_up(_ljc_num);}

	vcp_inline vcp_real_calc* dipoles_m_r_xBegin() const { return charges_m_r_xBegin() + _centers_m_r._round_up(_charges_num);}
	vcp_inline vcp_real_calc* dipoles_m_r_yBegin() const { return charges_m_r_yBegin() + _centers_m_r._round_up(_charges_num);}
	vcp_inline vcp_real_calc* dipoles_m_r_zBegin() const { return charges_m_r_zBegin() + _centers_m_r._round_up(_charges_num);}
	vcp_inline vcp_real_calc* dipoles_r_xBegin()   const { return charges_r_xBegin()   + _centers_r  ._round_up(_charges_num);}
	vcp_inline vcp_real_calc* dipoles_r_yBegin()   const { return charges_r_yBegin()   + _centers_r  ._round_up(_charges_num);}
	vcp_inline vcp_real_calc* dipoles_r_zBegin()   const { return charges_r_zBegin()   + _centers_r  ._round_up(_charges_num);}
	vcp_inline vcp_real_calc* dipoles_f_xBegin()   const { return charges_f_xBegin()   + _centers_f  ._round_up(_charges_num);}
	vcp_inline vcp_real_calc* dipoles_f_yBegin()   const { return charges_f_yBegin()   + _centers_f  ._round_up(_charges_num);}
	vcp_inline vcp_real_calc* dipoles_f_zBegin()   const { return charges_f_zBegin()   + _centers_f  ._round_up(_charges_num);}
	vcp_inline vcp_real_calc* dipoles_V_xBegin()   const { return charges_V_xBegin()   + _centers_V  ._round_up(_charges_num);}
	vcp_inline vcp_real_calc* dipoles_V_yBegin()   const { return charges_V_yBegin()   + _centers_V  ._round_up(_charges_num);}
	vcp_inline vcp_real_calc* dipoles_V_zBegin()   const { return charges_V_zBegin()   + _centers_V  ._round_up(_charges_num);}

	vcp_inline vcp_real_calc* quadrupoles_m_r_xBegin() const { return dipoles_m_r_xBegin() + _centers_m_r._round_up(_dipoles_num);}
	vcp_inline vcp_real_calc* quadrupoles_m_r_yBegin() const { return dipoles_m_r_yBegin() + _centers_m_r._round_up(_dipoles_num);}
	vcp_inline vcp_real_calc* quadrupoles_m_r_zBegin() const { return dipoles_m_r_zBegin() + _centers_m_r._round_up(_dipoles_num);}
	vcp_inline vcp_real_calc* quadrupoles_r_xBegin()   const { return dipoles_r_xBegin()   + _centers_r  ._round_up(_dipoles_num);}
	vcp_inline vcp_real_calc* quadrupoles_r_yBegin()   const { return dipoles_r_yBegin()   + _centers_r  ._round_up(_dipoles_num);}
	vcp_inline vcp_real_calc* quadrupoles_r_zBegin()   const { return dipoles_r_zBegin()   + _centers_r  ._round_up(_dipoles_num);}
	vcp_inline vcp_real_calc* quadrupoles_f_xBegin()   const { return dipoles_f_xBegin()   + _centers_f  ._round_up(_dipoles_num);}
	vcp_inline vcp_real_calc* quadrupoles_f_yBegin()   const { return dipoles_f_yBegin()   + _centers_f  ._round_up(_dipoles_num);}
	vcp_inline vcp_real_calc* quadrupoles_f_zBegin()   const { return dipoles_f_zBegin()   + _centers_f  ._round_up(_dipoles_num);}
	vcp_inline vcp_real_calc* quadrupoles_V_xBegin()   const { return dipoles_V_xBegin()   + _centers_V  ._round_up(_dipoles_num);}
	vcp_inline vcp_real_calc* quadrupoles_V_yBegin()   const { return dipoles_V_yBegin()   + _centers_V  ._round_up(_dipoles_num);}
	vcp_inline vcp_real_calc* quadrupoles_V_zBegin()   const { return dipoles_V_zBegin()   + _centers_V  ._round_up(_dipoles_num);}

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

	template<class T>
	vcp_inline
	void setPaddingToZero(AlignedArray<T>& t) const {
		size_t ljc_size = t._round_up(_ljc_num);
		size_t charges_size = t._round_up(_charges_num);
		size_t dipoles_size = t._round_up(_dipoles_num);

		t.zero(_ljc_num);
		t.zero(ljc_size + _charges_num);
		t.zero(ljc_size + charges_size + _dipoles_num);
		t.zero(ljc_size + charges_size + dipoles_size + _quadrupoles_num);
	}

	void resize(size_t molecules_arg, size_t ljcenters_arg, size_t charges_arg, size_t dipoles_arg, size_t quadrupoles_arg) {
		const bool allow_shrink = false; // TODO shrink at some point in the future

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

		// entries per center
		size_t num_centers_comp =
				AlignedArray<vcp_real_calc>::_round_up(_ljc_num) +
				AlignedArray<vcp_real_calc>::_round_up(_charges_num) +
				AlignedArray<vcp_real_calc>::_round_up(_dipoles_num) +
				AlignedArray<vcp_real_calc>::_round_up(_quadrupoles_num);

		size_t num_centers_acc =
				AlignedArray<vcp_real_calc>::_round_up(_ljc_num) +
				AlignedArray<vcp_real_calc>::_round_up(_charges_num) +
				AlignedArray<vcp_real_calc>::_round_up(_dipoles_num) +
				AlignedArray<vcp_real_calc>::_round_up(_quadrupoles_num);

		_centers_m_r.resize_zero_shrink(num_centers_comp);
		_centers_r	.resize_zero_shrink(num_centers_comp);
		_centers_f	.resize_zero_shrink(num_centers_acc );
		_centers_V	.resize_zero_shrink(num_centers_acc );

		// set padding to zero
		setPaddingToZero(_centers_m_r);
		setPaddingToZero(_centers_r);
		setPaddingToZero(_centers_f);
		setPaddingToZero(_centers_V);


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
};

#endif /* CELLDATASOA_H_ */
