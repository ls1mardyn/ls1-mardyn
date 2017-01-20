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
#include <cstdint>

/**
 * \brief Structure of Arrays for single-center lennard-Jones molecules for
 * the WR run.
 * \author Nikola Tchipev
 */
class CellDataSoA_WR {
public:
	CellDataSoA_WR(size_t mol_arg) {
		resize(mol_arg);
	}

	size_t _mol_num;

	// entries per molecule
	AlignedArrayTriplet<vcp_real_calc> _mol_r;
	AlignedArrayTriplet<vcp_real_calc> _mol_v;
	AlignedArray<uint64_t> _mol_uid;

	vcp_inline vcp_real_calc* r_xBegin()   const { return _mol_r.xBegin();}
	vcp_inline vcp_real_calc* r_yBegin()   const { return _mol_r.yBegin();}
	vcp_inline vcp_real_calc* r_zBegin()   const { return _mol_r.zBegin();}
	vcp_inline vcp_real_calc* v_xBegin()   const { return _mol_v.xBegin();}
	vcp_inline vcp_real_calc* v_yBegin()   const { return _mol_v.yBegin();}
	vcp_inline vcp_real_calc* v_zBegin()   const { return _mol_v.zBegin();}

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
};

#endif /* CELLDATASOA_WR_H_ */
