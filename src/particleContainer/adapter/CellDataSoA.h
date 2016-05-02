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
	typedef AlignedArray<size_t> IndexArray;
	typedef AlignedArray<double> DoubleArray;

	CellDataSoA(size_t molecules_arg, size_t lj_centers_arg, size_t charges_arg, size_t dipoles_arg, size_t quadrupoles_arg) :
		_mol_num(molecules_arg),
		_ljc_num(lj_centers_arg),
		_charges_num(charges_arg),
		_dipoles_num(dipoles_arg),
		_quadrupoles_num(quadrupoles_arg),
		_mol_size( molecules_arg + (molecules_arg & 1)),
		_ljc_size(lj_centers_arg + (lj_centers_arg & 1)),
		_charges_size(charges_arg + (charges_arg & 1)),
		_dipoles_size(dipoles_arg + (dipoles_arg & 1)),
		_quadrupoles_size(quadrupoles_arg + (quadrupoles_arg & 1)),
		_centers_size(_ljc_size + _charges_size + _dipoles_size + _quadrupoles_size),
		_mol_pos(_mol_size),
		_mol_ljc_num(_mol_size),
		_mol_charges_num(_mol_size),
		_mol_dipoles_num(_mol_size),
		_mol_quadrupoles_num(_mol_size),
		_centers_m_r(_centers_size),
		_centers_r(_centers_size),
		_centers_f(_centers_size),
		_centers_V(_centers_size),
		_ljc_id(_ljc_size),
		_charges_q(_charges_size),
		_dipoles_p(_dipoles_size),
		_dipoles_e(_dipoles_size),
		_dipoles_M(_dipoles_size),
		_quadrupoles_m(_quadrupoles_size),
		_quadrupoles_e(_quadrupoles_size),
		_quadrupoles_M(_quadrupoles_size)
		{}

	size_t _mol_num;
	size_t _ljc_num;
	size_t _charges_num;
	size_t _dipoles_num;
	size_t _quadrupoles_num;
	size_t _centers_num;

	size_t _mol_size;
	size_t _ljc_size;
	size_t _charges_size;
	size_t _dipoles_size;
	size_t _quadrupoles_size;
	size_t _centers_size;

	// entries per molecule
	AlignedArrayTriplet _mol_pos;
	AlignedArray<int> _mol_ljc_num;
	AlignedArray<int> _mol_charges_num;
	AlignedArray<int> _mol_dipoles_num;
	AlignedArray<int> _mol_quadrupoles_num;

	// entries per center
	AlignedArrayTriplet _centers_m_r;
	AlignedArrayTriplet _centers_r;
	AlignedArrayTriplet _centers_f;
	AlignedArrayTriplet _centers_V;

	double* _ljc_m_r_x;
	double* _ljc_m_r_y;
	double* _ljc_m_r_z;
	double* _ljc_r_x;
	double* _ljc_r_y;
	double* _ljc_r_z;
	double* _ljc_f_x;
	double* _ljc_f_y;
	double* _ljc_f_z;
	double* _ljc_V_x;
	double* _ljc_V_y;
	double* _ljc_V_z;

	double* _charges_m_r_x;
	double* _charges_m_r_y;
	double* _charges_m_r_z;
	double* _charges_r_x;
	double* _charges_r_y;
	double* _charges_r_z;
	double* _charges_f_x;
	double* _charges_f_y;
	double* _charges_f_z;
	double* _charges_V_x;
	double* _charges_V_y;
	double* _charges_V_z;

	double* _dipoles_m_r_x;
	double* _dipoles_m_r_y;
	double* _dipoles_m_r_z;
	double* _dipoles_r_x;
	double* _dipoles_r_y;
	double* _dipoles_r_z;
	double* _dipoles_f_x;
	double* _dipoles_f_y;
	double* _dipoles_f_z;
	double* _dipoles_V_x;
	double* _dipoles_V_y;
	double* _dipoles_V_z;

	double* _quadrupoles_m_r_x;
	double* _quadrupoles_m_r_y;
	double* _quadrupoles_m_r_z;
	double* _quadrupoles_r_x;
	double* _quadrupoles_r_y;
	double* _quadrupoles_r_z;
	double* _quadrupoles_f_x;
	double* _quadrupoles_f_y;
	double* _quadrupoles_f_z;
	double* _quadrupoles_V_x;
	double* _quadrupoles_V_y;
	double* _quadrupoles_V_z;

	// entries per lj center
	IndexArray _ljc_id;

	// entries per charge
	DoubleArray _charges_q;

	// entries per dipole
	DoubleArray _dipoles_p; // dipole moment
	AlignedArrayTriplet _dipoles_e; // orientation vector of dipole moment
	AlignedArrayTriplet _dipoles_M; // torque vector

	// entries per quadrupole
	DoubleArray _quadrupoles_m; // quadrupole moment
	AlignedArrayTriplet _quadrupoles_e; // orientation vector of quadrupole moment
	AlignedArrayTriplet _quadrupoles_M; // torque vector


	void vcp_inline initCenterPointers()
	{
		_ljc_m_r_x = _centers_m_r.xBegin();
		_ljc_m_r_y = _centers_m_r.yBegin();
		_ljc_m_r_z = _centers_m_r.zBegin();
		_ljc_r_x = _centers_r.xBegin();
		_ljc_r_y = _centers_r.yBegin();
		_ljc_r_z = _centers_r.zBegin();
		_ljc_f_x = _centers_f.xBegin();
		_ljc_f_y = _centers_f.yBegin();
		_ljc_f_z = _centers_f.zBegin();
		_ljc_V_x = _centers_V.xBegin();
		_ljc_V_y = _centers_V.yBegin();
		_ljc_V_z = _centers_V.zBegin();

		_charges_m_r_x = _ljc_m_r_x + _ljc_size;
		_charges_m_r_y = _ljc_m_r_y + _ljc_size;
		_charges_m_r_z = _ljc_m_r_z + _ljc_size;
		_charges_r_x = _ljc_r_x + _ljc_size;
		_charges_r_y = _ljc_r_y + _ljc_size;
		_charges_r_z = _ljc_r_z + _ljc_size;
		_charges_f_x = _ljc_f_x + _ljc_size;
		_charges_f_y = _ljc_f_y + _ljc_size;
		_charges_f_z = _ljc_f_z + _ljc_size;
		_charges_V_x = _ljc_V_x + _ljc_size;
		_charges_V_y = _ljc_V_y + _ljc_size;
		_charges_V_z = _ljc_V_z + _ljc_size;

		_dipoles_m_r_x = _charges_m_r_x + _charges_size;
		_dipoles_m_r_y = _charges_m_r_y + _charges_size;
		_dipoles_m_r_z = _charges_m_r_z + _charges_size;
		_dipoles_r_x = _charges_r_x + _charges_size;
		_dipoles_r_y = _charges_r_y + _charges_size;
		_dipoles_r_z = _charges_r_z + _charges_size;
		_dipoles_f_x = _charges_f_x + _charges_size;
		_dipoles_f_y = _charges_f_y + _charges_size;
		_dipoles_f_z = _charges_f_z + _charges_size;
		_dipoles_V_x = _charges_V_x + _charges_size;
		_dipoles_V_y = _charges_V_y + _charges_size;
		_dipoles_V_z = _charges_V_z + _charges_size;

		_quadrupoles_m_r_x = _dipoles_m_r_x + _dipoles_size;
		_quadrupoles_m_r_y = _dipoles_m_r_y + _dipoles_size;
		_quadrupoles_m_r_z = _dipoles_m_r_z + _dipoles_size;
		_quadrupoles_r_x = _dipoles_r_x + _dipoles_size;
		_quadrupoles_r_y = _dipoles_r_y + _dipoles_size;
		_quadrupoles_r_z = _dipoles_r_z + _dipoles_size;
		_quadrupoles_f_x = _dipoles_f_x + _dipoles_size;
		_quadrupoles_f_y = _dipoles_f_y + _dipoles_size;
		_quadrupoles_f_z = _dipoles_f_z + _dipoles_size;
		_quadrupoles_V_x = _dipoles_V_x + _dipoles_size;
		_quadrupoles_V_y = _dipoles_V_y + _dipoles_size;
		_quadrupoles_V_z = _dipoles_V_z + _dipoles_size;
	}

	void vcp_inline initDistLookupPointers(const AlignedArray<vcp_lookupOrMask_single>& centers_dist_lookup, vcp_lookupOrMask_single*& ljc_dist_lookup,
			vcp_lookupOrMask_single*& charges_dist_lookup, vcp_lookupOrMask_single*& dipoles_dist_lookup, vcp_lookupOrMask_single*& quadrupoles_dist_lookup) const{
		ljc_dist_lookup = centers_dist_lookup;
		charges_dist_lookup = ljc_dist_lookup + (_ljc_size + VCP_INDICES_PER_LOOKUP_SINGLE_M1)/VCP_INDICES_PER_LOOKUP_SINGLE;
		dipoles_dist_lookup = charges_dist_lookup + (_charges_size + VCP_INDICES_PER_LOOKUP_SINGLE_M1)/VCP_INDICES_PER_LOOKUP_SINGLE;
		quadrupoles_dist_lookup = dipoles_dist_lookup + (_dipoles_size + VCP_INDICES_PER_LOOKUP_SINGLE_M1)/VCP_INDICES_PER_LOOKUP_SINGLE;
	}

	template<class T>
	static vcp_inline
	void resizeLastZero(AlignedArray<T>& array, const size_t& size, const size_t& startZero){
		array.resize(size, startZero);
	}

	template<class T>
	vcp_inline
	void setPaddingToZero(T* ptr) const {
		//memset(array, 0, size * sizeof(T));//sets all to zero
		memset(ptr + _ljc_num, 0, (_ljc_size - _ljc_num) * sizeof(T)); //ljc
		ptr += _ljc_size;
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
		array.resize(size, size);
		T* ptr = array;
		setPaddingToZero(ptr);
	}

	/**
	 * resizes an array for all the centers and ensures, that the additionally allocated space is at least set once (valgrind error prevention reasons)
	 * @param array
	 * @param size
	 */
	vcp_inline
	void resizeCentersZero(AlignedArrayTriplet& triplet, const size_t& size) const{

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

			//TODO: "/ 8) * 8" things to be redone with AlignedArray::_round_up()

			if (_ljc_num > _ljc_size) {
				_ljc_size = ceil( (double)_ljc_num / 8) * 8;
				_ljc_id.resize(_ljc_size,_ljc_size);//set0 later on...
			}


			if (_charges_num > _charges_size) {
				_charges_size = ceil( (double)_charges_num / 8) * 8;
				resizeLastZero(_charges_q,_charges_size,_charges_num);
			}

			if (_dipoles_num > _dipoles_size) {
				_dipoles_size = ceil( (double)_dipoles_num / 8) * 8;
				resizeLastZero(_dipoles_p,_dipoles_size, _dipoles_num);
				_dipoles_e.resize(_dipoles_size);
				_dipoles_M.resize(_dipoles_size);
			}

			if (_quadrupoles_num > _quadrupoles_size) {
				_quadrupoles_size = ceil( (double)_quadrupoles_num / 8) * 8;
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
			_mol_size = ceil( (double)molecules_arg / 8) * 8;
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

		initCenterPointers();
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
