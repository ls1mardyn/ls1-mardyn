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

	CellDataSoA(size_t molecules_arg, size_t lj_centers_arg, size_t charges_arg, size_t dipoles_arg, size_t quadrupoles_arg) :
		_mol_num(molecules_arg),
		_ljc_num(lj_centers_arg),
		_charges_num(charges_arg),
		_dipoles_num(dipoles_arg),
		_quadrupoles_num(quadrupoles_arg),
		_centers_num(lj_centers_arg + charges_arg + dipoles_arg + quadrupoles_arg),
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
	IndexArray _ljc_id;

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
	vcp_inline double* ljc_r_xBegin()   const { return _centers_r.xBegin()  ;}
	vcp_inline double* ljc_r_yBegin()   const { return _centers_r.yBegin()  ;}
	vcp_inline double* ljc_r_zBegin()   const { return _centers_r.zBegin()  ;}
	vcp_inline double* ljc_f_xBegin()   const { return _centers_f.xBegin()  ;}
	vcp_inline double* ljc_f_yBegin()   const { return _centers_f.yBegin()  ;}
	vcp_inline double* ljc_f_zBegin()   const { return _centers_f.zBegin()  ;}
	vcp_inline double* ljc_V_xBegin()   const { return _centers_V.xBegin()  ;}
	vcp_inline double* ljc_V_yBegin()   const { return _centers_V.yBegin()  ;}
	vcp_inline double* ljc_V_zBegin()   const { return _centers_V.zBegin()  ;}

	vcp_inline double* charges_m_r_xBegin() const { return ljc_m_r_xBegin() + _ljc_size;}
	vcp_inline double* charges_m_r_yBegin() const { return ljc_m_r_yBegin() + _ljc_size;}
	vcp_inline double* charges_m_r_zBegin() const { return ljc_m_r_zBegin() + _ljc_size;}
	vcp_inline double* charges_r_xBegin()   const { return ljc_r_xBegin()   + _ljc_size;}
	vcp_inline double* charges_r_yBegin()   const { return ljc_r_yBegin()   + _ljc_size;}
	vcp_inline double* charges_r_zBegin()   const { return ljc_r_zBegin()   + _ljc_size;}
	vcp_inline double* charges_f_xBegin()   const { return ljc_f_xBegin()   + _ljc_size;}
	vcp_inline double* charges_f_yBegin()   const { return ljc_f_yBegin()   + _ljc_size;}
	vcp_inline double* charges_f_zBegin()   const { return ljc_f_zBegin()   + _ljc_size;}
	vcp_inline double* charges_V_xBegin()   const { return ljc_V_xBegin()   + _ljc_size;}
	vcp_inline double* charges_V_yBegin()   const { return ljc_V_yBegin()   + _ljc_size;}
	vcp_inline double* charges_V_zBegin()   const { return ljc_V_zBegin()   + _ljc_size;}

	vcp_inline double* dipoles_m_r_xBegin() const { return charges_m_r_xBegin() + _charges_size;}
	vcp_inline double* dipoles_m_r_yBegin() const { return charges_m_r_yBegin() + _charges_size;}
	vcp_inline double* dipoles_m_r_zBegin() const { return charges_m_r_zBegin() + _charges_size;}
	vcp_inline double* dipoles_r_xBegin()   const { return charges_r_xBegin()   + _charges_size;}
	vcp_inline double* dipoles_r_yBegin()   const { return charges_r_yBegin()   + _charges_size;}
	vcp_inline double* dipoles_r_zBegin()   const { return charges_r_zBegin()   + _charges_size;}
	vcp_inline double* dipoles_f_xBegin()   const { return charges_f_xBegin()   + _charges_size;}
	vcp_inline double* dipoles_f_yBegin()   const { return charges_f_yBegin()   + _charges_size;}
	vcp_inline double* dipoles_f_zBegin()   const { return charges_f_zBegin()   + _charges_size;}
	vcp_inline double* dipoles_V_xBegin()   const { return charges_V_xBegin()   + _charges_size;}
	vcp_inline double* dipoles_V_yBegin()   const { return charges_V_yBegin()   + _charges_size;}
	vcp_inline double* dipoles_V_zBegin()   const { return charges_V_zBegin()   + _charges_size;}

	vcp_inline double* quadrupoles_m_r_xBegin() const { return dipoles_m_r_xBegin() + _dipoles_size;}
	vcp_inline double* quadrupoles_m_r_yBegin() const { return dipoles_m_r_yBegin() + _dipoles_size;}
	vcp_inline double* quadrupoles_m_r_zBegin() const { return dipoles_m_r_zBegin() + _dipoles_size;}
	vcp_inline double* quadrupoles_r_xBegin()   const { return dipoles_r_xBegin()   + _dipoles_size;}
	vcp_inline double* quadrupoles_r_yBegin()   const { return dipoles_r_yBegin()   + _dipoles_size;}
	vcp_inline double* quadrupoles_r_zBegin()   const { return dipoles_r_zBegin()   + _dipoles_size;}
	vcp_inline double* quadrupoles_f_xBegin()   const { return dipoles_f_xBegin()   + _dipoles_size;}
	vcp_inline double* quadrupoles_f_yBegin()   const { return dipoles_f_yBegin()   + _dipoles_size;}
	vcp_inline double* quadrupoles_f_zBegin()   const { return dipoles_f_zBegin()   + _dipoles_size;}
	vcp_inline double* quadrupoles_V_xBegin()   const { return dipoles_V_xBegin()   + _dipoles_size;}
	vcp_inline double* quadrupoles_V_yBegin()   const { return dipoles_V_yBegin()   + _dipoles_size;}
	vcp_inline double* quadrupoles_V_zBegin()   const { return dipoles_V_zBegin()   + _dipoles_size;}

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
		array.resize(size);
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

			//TODO: "/ 8) * 8" things to be redone with AlignedArray::_round_up()

			if (_ljc_num > _ljc_size) {
				_ljc_size = ceil( (double)_ljc_num / 8) * 8;
				_ljc_id.resize(_ljc_size);//set0 later on...
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
