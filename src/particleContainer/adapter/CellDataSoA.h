/*
 * CellDataSoA.h
 *
 * @Date: 25.03.2013
 * @Author: eckhardw
 */

#ifndef LENNARDJONESSOA_H_
#define LENNARDJONESSOA_H_

#include "utils/AlignedArray.h"

/**
 * \brief Structure of Arrays for vectorized Lennard Jones force calculation.
 * \author Johannes Heckl, Wolfgang Eckhardt
 */
class CellDataSoA {
public:
	typedef AlignedArray<size_t> IndexArray;
	typedef AlignedArray<double> DoubleArray;

	CellDataSoA(size_t molecules_arg, size_t centers_arg) :
		_num_molecules(molecules_arg),
		_num_ljcenters(centers_arg),
		_molecules_size( molecules_arg + (molecules_arg & 1)),
		_ljcenters_size(centers_arg + (centers_arg & 1)),
		_mol_pos_x(_molecules_size), _mol_pos_y(_molecules_size), _mol_pos_z(_molecules_size), _mol_num_ljc(_molecules_size),
		_m_r_x(_ljcenters_size), _m_r_y(_ljcenters_size), _m_r_z(
		_ljcenters_size), _ljc_r_x(_ljcenters_size), _ljc_r_y(_ljcenters_size), _ljc_r_z(_ljcenters_size), _ljc_f_x(
		_ljcenters_size), _ljc_f_y(_ljcenters_size), _ljc_f_z(_ljcenters_size), _ljc_id(_ljcenters_size) {

	}

	size_t _num_molecules;
	size_t _num_ljcenters;
	size_t _molecules_size;
	size_t _ljcenters_size;
	// entries per molecule
	DoubleArray _mol_pos_x;
	DoubleArray _mol_pos_y;
	DoubleArray _mol_pos_z;
	AlignedArray<int> _mol_num_ljc;

	// entries per center
	DoubleArray _m_r_x;
	DoubleArray _m_r_y;
	DoubleArray _m_r_z;
	DoubleArray _ljc_r_x;
	DoubleArray _ljc_r_y;
	DoubleArray _ljc_r_z;
	DoubleArray _ljc_f_x;
	DoubleArray _ljc_f_y;
	DoubleArray _ljc_f_z;
	IndexArray _ljc_id;

	void resize(size_t molecules_arg, size_t centers_arg) {
			if (centers_arg > _ljcenters_size) {
				_ljcenters_size = ceil( (double)centers_arg / 4) * 4;
				_m_r_x.resize(_ljcenters_size);
				_m_r_y.resize(_ljcenters_size);
				_m_r_z.resize(_ljcenters_size);
				_ljc_r_x.resize(_ljcenters_size);
				_ljc_r_y.resize(_ljcenters_size);
				_ljc_r_z.resize(_ljcenters_size);
				_ljc_f_x.resize(_ljcenters_size);
				_ljc_f_y.resize(_ljcenters_size);
				_ljc_f_z.resize(_ljcenters_size);
				_ljc_id.resize(_ljcenters_size);
			}

			if (molecules_arg > _molecules_size) {
				_molecules_size = ceil( (double)molecules_arg / 4) * 4;
				_mol_pos_x.resize(_molecules_size);
				_mol_pos_y.resize(_molecules_size);
				_mol_pos_z.resize(_molecules_size);
				_mol_num_ljc.resize(_molecules_size);
			}
		}
};

#endif /* LENNARDJONESSOA_H_ */
