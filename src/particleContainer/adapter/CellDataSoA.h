/*
 * CellDataSoA.h
 *
 * @Date: 25.03.2013
 * @Author: eckhardw
 */

#ifndef CELLDATASOA_H_
#define CELLDATASOA_H_

#include "utils/AlignedArray.h"

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
		_mol_pos_x(_mol_size), _mol_pos_y(_mol_size), _mol_pos_z(_mol_size),
		_mol_ljc_num(_mol_size),
		_mol_charges_num(_mol_size),
		_mol_dipoles_num(_mol_size),
		_mol_quadrupoles_num(_mol_size),
		_centers_m_r_x(_centers_size), _centers_m_r_y(_centers_size), _centers_m_r_z(_centers_size),
		_centers_r_x(_centers_size), _centers_r_y(_centers_size), _centers_r_z(_centers_size),
		_centers_f_x(_centers_size), _centers_f_y(_centers_size), _centers_f_z(_centers_size),
		_ljc_id(_ljc_size),
		_charges_q(_charges_size),
		_dipoles_p(_dipoles_size),
		_dipoles_e_x(_dipoles_size), _dipoles_e_y(_dipoles_size), _dipoles_e_z(_dipoles_size),
		_dipoles_M_x(_dipoles_size), _dipoles_M_y(_dipoles_size), _dipoles_M_z(_dipoles_size),
		_quadrupoles_m(_quadrupoles_size),
		_quadrupoles_e_x(_quadrupoles_size), _quadrupoles_e_y(_quadrupoles_size), _quadrupoles_e_z(_quadrupoles_size),
		_quadrupoles_M_x(_quadrupoles_size), _quadrupoles_M_y(_quadrupoles_size), _quadrupoles_M_z(_quadrupoles_size)
		{}

	size_t _mol_num;
	size_t _ljc_num;
	size_t _charges_num;
	size_t _dipoles_num;
	size_t _quadrupoles_num;

	size_t _mol_size;
	size_t _ljc_size;
	size_t _charges_size;
	size_t _dipoles_size;
	size_t _quadrupoles_size;
	size_t _centers_size;

	// entries per molecule
	DoubleArray _mol_pos_x;
	DoubleArray _mol_pos_y;
	DoubleArray _mol_pos_z;
	AlignedArray<int> _mol_ljc_num;
	AlignedArray<int> _mol_charges_num;
	AlignedArray<int> _mol_dipoles_num;
	AlignedArray<int> _mol_quadrupoles_num;

	// entries per center
	DoubleArray _centers_m_r_x;
	DoubleArray _centers_m_r_y;
	DoubleArray _centers_m_r_z;
	DoubleArray _centers_r_x;
	DoubleArray _centers_r_y;
	DoubleArray _centers_r_z;
	DoubleArray _centers_f_x;
	DoubleArray _centers_f_y;
	DoubleArray _centers_f_z;
	DoubleArray _centers_dist_lookup;

	double* _ljc_m_r_x;
	double* _ljc_m_r_y;
	double* _ljc_m_r_z;
	double* _ljc_r_x;
	double* _ljc_r_y;
	double* _ljc_r_z;
	double* _ljc_f_x;
	double* _ljc_f_y;
	double* _ljc_f_z;
	double* _ljc_dist_lookup;

	double* _charges_m_r_x;
	double* _charges_m_r_y;
	double* _charges_m_r_z;
	double* _charges_r_x;
	double* _charges_r_y;
	double* _charges_r_z;
	double* _charges_f_x;
	double* _charges_f_y;
	double* _charges_f_z;
	double* _charges_dist_lookup;

	double* _dipoles_m_r_x;
	double* _dipoles_m_r_y;
	double* _dipoles_m_r_z;
	double* _dipoles_r_x;
	double* _dipoles_r_y;
	double* _dipoles_r_z;
	double* _dipoles_f_x;
	double* _dipoles_f_y;
	double* _dipoles_f_z;
	double* _dipoles_dist_lookup;

	double* _quadrupoles_m_r_x;
	double* _quadrupoles_m_r_y;
	double* _quadrupoles_m_r_z;
	double* _quadrupoles_r_x;
	double* _quadrupoles_r_y;
	double* _quadrupoles_r_z;
	double* _quadrupoles_f_x;
	double* _quadrupoles_f_y;
	double* _quadrupoles_f_z;
	double* _quadrupoles_dist_lookup;

	// entries per lj center
	IndexArray _ljc_id;

	// entries per charge
	DoubleArray _charges_q;

	// entries per dipole
	DoubleArray _dipoles_p; // dipole moment
	DoubleArray _dipoles_e_x; // orientation vector of dipole moment
	DoubleArray _dipoles_e_y;
	DoubleArray _dipoles_e_z;
	DoubleArray _dipoles_M_x; // torque vector
	DoubleArray _dipoles_M_y;
	DoubleArray _dipoles_M_z;

	// entries per quadrupole
	DoubleArray _quadrupoles_m; // quadrupole moment
	DoubleArray _quadrupoles_e_x; // orientation vector of quadrupole moment
	DoubleArray _quadrupoles_e_y;
	DoubleArray _quadrupoles_e_z;
	DoubleArray _quadrupoles_M_x; // torque vector
	DoubleArray _quadrupoles_M_y;
	DoubleArray _quadrupoles_M_z;


	void initCenterPointers()
	{
		_ljc_m_r_x = _centers_m_r_x;
		_ljc_m_r_y = _centers_m_r_y;
		_ljc_m_r_z = _centers_m_r_z;
		_ljc_r_x = _centers_r_x;
		_ljc_r_y = _centers_r_y;
		_ljc_r_z = _centers_r_z;
		_ljc_f_x = _centers_f_x;
		_ljc_f_y = _centers_f_y;
		_ljc_f_z = _centers_f_z;
		_ljc_dist_lookup = _centers_dist_lookup;

		_charges_m_r_x = _ljc_m_r_x + _ljc_size;
		_charges_m_r_y = _ljc_m_r_y + _ljc_size;
		_charges_m_r_z = _ljc_m_r_z + _ljc_size;
		_charges_r_x = _ljc_r_x + _ljc_size;
		_charges_r_y = _ljc_r_y + _ljc_size;
		_charges_r_z = _ljc_r_z + _ljc_size;
		_charges_f_x = _ljc_f_x + _ljc_size;
		_charges_f_y = _ljc_f_y + _ljc_size;
		_charges_f_z = _ljc_f_z + _ljc_size;
		_charges_dist_lookup = _ljc_dist_lookup + _ljc_size;

		_dipoles_m_r_x = _charges_m_r_x + _charges_size;
		_dipoles_m_r_y = _charges_m_r_y + _charges_size;
		_dipoles_m_r_z = _charges_m_r_z + _charges_size;
		_dipoles_r_x = _charges_r_x + _charges_size;
		_dipoles_r_y = _charges_r_y + _charges_size;
		_dipoles_r_z = _charges_r_z + _charges_size;
		_dipoles_f_x = _charges_f_x + _charges_size;
		_dipoles_f_y = _charges_f_y + _charges_size;
		_dipoles_f_z = _charges_f_z + _charges_size;
		_dipoles_dist_lookup = _charges_dist_lookup + _charges_size;

		_quadrupoles_m_r_x = _dipoles_m_r_x + _dipoles_size;
		_quadrupoles_m_r_y = _dipoles_m_r_y + _dipoles_size;
		_quadrupoles_m_r_z = _dipoles_m_r_z + _dipoles_size;
		_quadrupoles_r_x = _dipoles_r_x + _dipoles_size;
		_quadrupoles_r_y = _dipoles_r_y + _dipoles_size;
		_quadrupoles_r_z = _dipoles_r_z + _dipoles_size;
		_quadrupoles_f_x = _dipoles_f_x + _dipoles_size;
		_quadrupoles_f_y = _dipoles_f_y + _dipoles_size;
		_quadrupoles_f_z = _dipoles_f_z + _dipoles_size;
		_quadrupoles_dist_lookup = _dipoles_dist_lookup + _dipoles_size;
	}

	void setDistLookup(DoubleArray & centers_dist_lookup_arg)
	{
		_centers_dist_lookup = centers_dist_lookup_arg;
		resizeDistLookup();
	}

	void resizeDistLookup()
	{
		if (_centers_dist_lookup.get_size() < _centers_size)
		{
			_centers_dist_lookup.resize(_centers_size);
		}
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
				_ljc_size = ceil( (double)_ljc_num / 8) * 8;
				_ljc_id.resize(_ljc_size);
				memset(_ljc_id + _ljc_num, 0, (_ljc_size-_ljc_num) * sizeof(double));//set the remaining values to zero.
				//This is needed to allow vectorization even of the last elements, their count does not necessarily divide by VCP_VEC_SIZE.
				//The array size is however long enough to vectorize over the last few entries.
				//This sets the entries, that do not make sense in that vectorization to zero. In this case this is needed to allow indirect access using this vector.
			}

			if (_charges_num > _charges_size) {
				_charges_size = ceil( (double)_charges_num / 8) * 8;
				_charges_q.resize(_charges_size);
			}

			if (_dipoles_num > _dipoles_size) {
				_dipoles_size = ceil( (double)_dipoles_num / 8) * 8;
				_dipoles_p.resize(_dipoles_size);
				_dipoles_e_x.resize(_dipoles_size);
				_dipoles_e_y.resize(_dipoles_size);
				_dipoles_e_z.resize(_dipoles_size);
				_dipoles_M_x.resize(_dipoles_size);
				_dipoles_M_y.resize(_dipoles_size);
				_dipoles_M_z.resize(_dipoles_size);
			}

			if (_quadrupoles_num > _quadrupoles_size) {
				_quadrupoles_size = ceil( (double)_quadrupoles_num / 8) * 8;
				_quadrupoles_m.resize(_quadrupoles_size);
				_quadrupoles_e_x.resize(_quadrupoles_size);
				_quadrupoles_e_y.resize(_quadrupoles_size);
				_quadrupoles_e_z.resize(_quadrupoles_size);
				_quadrupoles_M_x.resize(_quadrupoles_size);
				_quadrupoles_M_y.resize(_quadrupoles_size);
				_quadrupoles_M_z.resize(_quadrupoles_size);
			}

			if (_centers_size < _ljc_size + _charges_size + _dipoles_size + _quadrupoles_size)
			{
				_centers_size = _ljc_size + _charges_size + _dipoles_size + _quadrupoles_size;

				_centers_m_r_x.resize(_centers_size);
				_centers_m_r_y.resize(_centers_size);
				_centers_m_r_z.resize(_centers_size);
				_centers_r_x.resize(_centers_size);
				_centers_r_y.resize(_centers_size);
				_centers_r_z.resize(_centers_size);
				_centers_f_x.resize(_centers_size);
				_centers_f_y.resize(_centers_size);
				_centers_f_z.resize(_centers_size);
			}
		}

		if (_mol_num > _mol_size) {
			_mol_size = ceil( (double)molecules_arg / 8) * 8;
			_mol_pos_x.resize(_mol_size);
			_mol_pos_y.resize(_mol_size);
			_mol_pos_z.resize(_mol_size);
			_mol_ljc_num.resize(_mol_size);
			_mol_charges_num.resize(_mol_size);
			_mol_dipoles_num.resize(_mol_size);
			_mol_quadrupoles_num.resize(_mol_size);
		}

		resizeDistLookup();
		initCenterPointers();
	}
};

#endif /* CELLDATASOA_H_ */
