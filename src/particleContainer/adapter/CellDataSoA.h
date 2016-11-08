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
	AlignedArrayTriplet<double> _mol_pos;
	AlignedArray<int> _mol_ljc_num;
	AlignedArray<int> _mol_charges_num;
	AlignedArray<int> _mol_dipoles_num;
	AlignedArray<int> _mol_quadrupoles_num;

	// entries per center
	AlignedArrayTriplet<double> _centers_m_r;
	AlignedArrayTriplet<double> _centers_r;
	AlignedArrayTriplet<double> _centers_f;
	AlignedArrayTriplet<double> _centers_V;

	// entries per lj center
	AlignedArray<size_t> _ljc_id;

	// entries per charge
	AlignedArray<double> _charges_q;

	// entries per dipole
	AlignedArray<double> _dipoles_p; // dipole moment
	AlignedArrayTriplet<double> _dipoles_e; // orientation vector of dipole moment
	AlignedArrayTriplet<double> _dipoles_M; // torque vector

	// entries per quadrupole
	AlignedArray<double> _quadrupoles_m; // quadrupole moment
	AlignedArrayTriplet<double> _quadrupoles_e; // orientation vector of quadrupole moment
	AlignedArrayTriplet<double> _quadrupoles_M; // torque vector

	vcp_inline double* ljc_m_r_xBegin() const { return _centers_m_r.xBegin();}
	vcp_inline double* ljc_m_r_yBegin() const { return _centers_m_r.yBegin();}
	vcp_inline double* ljc_m_r_zBegin() const { return _centers_m_r.zBegin();}
	vcp_inline double* ljc_r_xBegin()   const { return _centers_r  .xBegin();}
	vcp_inline double* ljc_r_yBegin()   const { return _centers_r  .yBegin();}
	vcp_inline double* ljc_r_zBegin()   const { return _centers_r  .zBegin();}
	vcp_inline double* ljc_f_xBegin()   const { return _centers_f  .xBegin();}
	vcp_inline double* ljc_f_yBegin()   const { return _centers_f  .yBegin();}
	vcp_inline double* ljc_f_zBegin()   const { return _centers_f  .zBegin();}
	vcp_inline double* ljc_V_xBegin()   const { return _centers_V  .xBegin();}
	vcp_inline double* ljc_V_yBegin()   const { return _centers_V  .yBegin();}
	vcp_inline double* ljc_V_zBegin()   const { return _centers_V  .zBegin();}

	vcp_inline double* charges_m_r_xBegin() const { return ljc_m_r_xBegin() + _centers_m_r._round_up(_ljc_num);}
	vcp_inline double* charges_m_r_yBegin() const { return ljc_m_r_yBegin() + _centers_m_r._round_up(_ljc_num);}
	vcp_inline double* charges_m_r_zBegin() const { return ljc_m_r_zBegin() + _centers_m_r._round_up(_ljc_num);}
	vcp_inline double* charges_r_xBegin()   const { return ljc_r_xBegin()   + _centers_r  ._round_up(_ljc_num);}
	vcp_inline double* charges_r_yBegin()   const { return ljc_r_yBegin()   + _centers_r  ._round_up(_ljc_num);}
	vcp_inline double* charges_r_zBegin()   const { return ljc_r_zBegin()   + _centers_r  ._round_up(_ljc_num);}
	vcp_inline double* charges_f_xBegin()   const { return ljc_f_xBegin()   + _centers_f  ._round_up(_ljc_num);}
	vcp_inline double* charges_f_yBegin()   const { return ljc_f_yBegin()   + _centers_f  ._round_up(_ljc_num);}
	vcp_inline double* charges_f_zBegin()   const { return ljc_f_zBegin()   + _centers_f  ._round_up(_ljc_num);}
	vcp_inline double* charges_V_xBegin()   const { return ljc_V_xBegin()   + _centers_V  ._round_up(_ljc_num);}
	vcp_inline double* charges_V_yBegin()   const { return ljc_V_yBegin()   + _centers_V  ._round_up(_ljc_num);}
	vcp_inline double* charges_V_zBegin()   const { return ljc_V_zBegin()   + _centers_V  ._round_up(_ljc_num);}

	vcp_inline double* dipoles_m_r_xBegin() const { return charges_m_r_xBegin() + _centers_m_r._round_up(_charges_num);}
	vcp_inline double* dipoles_m_r_yBegin() const { return charges_m_r_yBegin() + _centers_m_r._round_up(_charges_num);}
	vcp_inline double* dipoles_m_r_zBegin() const { return charges_m_r_zBegin() + _centers_m_r._round_up(_charges_num);}
	vcp_inline double* dipoles_r_xBegin()   const { return charges_r_xBegin()   + _centers_r  ._round_up(_charges_num);}
	vcp_inline double* dipoles_r_yBegin()   const { return charges_r_yBegin()   + _centers_r  ._round_up(_charges_num);}
	vcp_inline double* dipoles_r_zBegin()   const { return charges_r_zBegin()   + _centers_r  ._round_up(_charges_num);}
	vcp_inline double* dipoles_f_xBegin()   const { return charges_f_xBegin()   + _centers_f  ._round_up(_charges_num);}
	vcp_inline double* dipoles_f_yBegin()   const { return charges_f_yBegin()   + _centers_f  ._round_up(_charges_num);}
	vcp_inline double* dipoles_f_zBegin()   const { return charges_f_zBegin()   + _centers_f  ._round_up(_charges_num);}
	vcp_inline double* dipoles_V_xBegin()   const { return charges_V_xBegin()   + _centers_V  ._round_up(_charges_num);}
	vcp_inline double* dipoles_V_yBegin()   const { return charges_V_yBegin()   + _centers_V  ._round_up(_charges_num);}
	vcp_inline double* dipoles_V_zBegin()   const { return charges_V_zBegin()   + _centers_V  ._round_up(_charges_num);}

	vcp_inline double* quadrupoles_m_r_xBegin() const { return dipoles_m_r_xBegin() + _centers_m_r._round_up(_dipoles_num);}
	vcp_inline double* quadrupoles_m_r_yBegin() const { return dipoles_m_r_yBegin() + _centers_m_r._round_up(_dipoles_num);}
	vcp_inline double* quadrupoles_m_r_zBegin() const { return dipoles_m_r_zBegin() + _centers_m_r._round_up(_dipoles_num);}
	vcp_inline double* quadrupoles_r_xBegin()   const { return dipoles_r_xBegin()   + _centers_r  ._round_up(_dipoles_num);}
	vcp_inline double* quadrupoles_r_yBegin()   const { return dipoles_r_yBegin()   + _centers_r  ._round_up(_dipoles_num);}
	vcp_inline double* quadrupoles_r_zBegin()   const { return dipoles_r_zBegin()   + _centers_r  ._round_up(_dipoles_num);}
	vcp_inline double* quadrupoles_f_xBegin()   const { return dipoles_f_xBegin()   + _centers_f  ._round_up(_dipoles_num);}
	vcp_inline double* quadrupoles_f_yBegin()   const { return dipoles_f_yBegin()   + _centers_f  ._round_up(_dipoles_num);}
	vcp_inline double* quadrupoles_f_zBegin()   const { return dipoles_f_zBegin()   + _centers_f  ._round_up(_dipoles_num);}
	vcp_inline double* quadrupoles_V_xBegin()   const { return dipoles_V_xBegin()   + _centers_V  ._round_up(_dipoles_num);}
	vcp_inline double* quadrupoles_V_yBegin()   const { return dipoles_V_yBegin()   + _centers_V  ._round_up(_dipoles_num);}
	vcp_inline double* quadrupoles_V_zBegin()   const { return dipoles_V_zBegin()   + _centers_V  ._round_up(_dipoles_num);}

	void vcp_inline initDistLookupPointers(
			AlignedArray<vcp_lookupOrMask_single>& centers_dist_lookup,
			vcp_lookupOrMask_single*& ljc_dist_lookup,
			vcp_lookupOrMask_single*& charges_dist_lookup,
			vcp_lookupOrMask_single*& dipoles_dist_lookup,
			vcp_lookupOrMask_single*& quadrupoles_dist_lookup) const {

		size_t ljc_size 	= AlignedArray<double>::_round_up(_ljc_num);
		size_t charges_size = AlignedArray<double>::_round_up(_charges_num);
		size_t dipoles_size = AlignedArray<double>::_round_up(_dipoles_num);
		size_t quadrupoles_size = AlignedArray<double>::_round_up(_quadrupoles_num);
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
				AlignedArray<double>::_round_up(_ljc_num) +
				AlignedArray<double>::_round_up(_charges_num) +
				AlignedArray<double>::_round_up(_dipoles_num) +
				AlignedArray<double>::_round_up(_quadrupoles_num);

		size_t num_centers_acc =
				AlignedArray<double>::_round_up(_ljc_num) +
				AlignedArray<double>::_round_up(_charges_num) +
				AlignedArray<double>::_round_up(_dipoles_num) +
				AlignedArray<double>::_round_up(_quadrupoles_num);

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

#if 0

	template<class T>
	static vcp_inline
	void resizeLastZero(AlignedArray<T>& array, const size_t& size, const size_t& startZero){
		array.resize(size);
	}

	template<class T>
	vcp_inline
	void setPaddingToZero(T* ptr) const {
		size_t ljc_size = AlignedArray<double>::_round_up(_ljc_num);

		//memset(array, 0, size * sizeof(T));//sets all to zero
		memset(ptr + _ljc_num, 0, (ljc_size - _ljc_num) * sizeof(T)); //ljc
		ptr += ljc_size;
		memset(ptr + _charges_num, 0, (_charges_size - _charges_num) * sizeof(T)); //charges
		ptr += _charges_size;
		memset(ptr + _dipoles_num, 0, (_dipoles_size - _dipoles_num) * sizeof(T)); //dipoles
		ptr += _dipoles_size;
		memset(ptr + _quadrupoles_num, 0, (_quadrupoles_size - _quadrupoles_num) * sizeof(T)); //quadrupoles
	}

	/**
	 * resizes an array for all the centers and ensures, that the additionally allocated space is at least set once (valgrind error prevention reasons)
	 * @param array
	 * @param size
	 */
	template<class T>
	vcp_inline
	void resizeCentersZero(AlignedArray<T>& array, const size_t& size) const{
		array.resize(size);
		T* ptr = array;
		setPaddingToZero(ptr);
	}

	/**
	 * resizes an array for all the centers and ensures, that the additionally allocated space is at least set once (valgrind error prevention reasons)
	 * @param array
	 * @param size
	 */
	template<class T>
	vcp_inline
	void resizeCentersZero(AlignedArrayTriplet<T>& triplet, const size_t& size) const{

		triplet.resize(size);

		setPaddingToZero(triplet.xBegin());
		setPaddingToZero(triplet.yBegin());
		setPaddingToZero(triplet.zBegin());
	}

	void resize(size_t molecules_arg, size_t ljcenters_arg, size_t charges_arg, size_t dipoles_arg, size_t quadrupoles_arg) {
		_mol_num = molecules_arg;
		_ljc_num = ljcenters_arg;
		_charges_num = charges_arg;
		_dipoles_num = dipoles_arg;
		_quadrupoles_num = quadrupoles_arg;

		if (_ljc_num > _ljc_size ||
			_charges_num > _charges_size ||
			_dipoles_num > _dipoles_size ||
			_quadrupoles_num > _quadrupoles_size ) {

			if (_ljc_num > _ljc_size) {
				_ljc_size = AlignedArray<double>::_round_up(_ljc_num);
				_ljc_id.resize(_ljc_size);//set0 later on...
			}


			if (_charges_num > _charges_size) {
				_charges_size = AlignedArray<double>::_round_up(_charges_num);
				resizeLastZero(_charges_q,_charges_size,_charges_num);
			}

			if (_dipoles_num > _dipoles_size) {
				_dipoles_size = AlignedArray<double>::_round_up(_dipoles_num);
				resizeLastZero(_dipoles_p,_dipoles_size, _dipoles_num);
				_dipoles_e.resize(_dipoles_size);
				_dipoles_M.resize(_dipoles_size);
			}

			if (_quadrupoles_num > _quadrupoles_size) {
				_quadrupoles_size = AlignedArray<double>::_round_up(_quadrupoles_num);
				resizeLastZero(_quadrupoles_m,_quadrupoles_size, _quadrupoles_num);
				_quadrupoles_e.resize(_quadrupoles_size);
				_quadrupoles_M.resize(_quadrupoles_size);
			}

			if (_centers_size < _ljc_size + _charges_size + _dipoles_size + _quadrupoles_size)
			{
				_centers_size = _ljc_size + _charges_size + _dipoles_size + _quadrupoles_size;//divisible by 8, since all others are
				_centers_num = _ljc_num + _charges_num + _dipoles_num + _quadrupoles_num;
				resizeCentersZero(_centers_m_r, _centers_size);
				resizeCentersZero(_centers_r, _centers_size);
				resizeCentersZero(_centers_f, _centers_size);
				resizeCentersZero(_centers_V, _centers_size);
			}
		}

		if (_mol_num > _mol_size) {
			_mol_size = AlignedArray<double>::_round_up(_mol_num);
			_mol_pos.resize(_mol_size);
			resizeLastZero(_mol_ljc_num,_mol_size, _mol_num);
			resizeLastZero(_mol_charges_num,_mol_size, _mol_num);
			resizeLastZero(_mol_dipoles_num,_mol_size, _mol_num);
			resizeLastZero(_mol_quadrupoles_num,_mol_size, _mol_num);
		}
		memset(_ljc_id + _ljc_num, 0, (_ljc_size - _ljc_num) * sizeof(size_t));//set the remaining values to zero.
			//This is needed to allow vectorization even of the last elements, their count does not necessarily divide by VCP_VEC_SIZE.
			//The array size is however long enough to vectorize over the last few entries.
			//This sets the entries, that do not make sense in that vectorization to zero. In this case this is needed to allow indirect access using this vector.
	}
#endif /* 0 */

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
