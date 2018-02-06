/*
 * PositionCellProcessorRMM.h
 *
 *  Created on: 5 Oct 2017
 *      Author: tchipevn
 */

#ifndef SRC_INTEGRATORS_POSITIONCELLPROCESSORRMM_H_
#define SRC_INTEGRATORS_POSITIONCELLPROCESSORRMM_H_

#include "particleContainer/ParticleCell.h"
#include "particleContainer/adapter/CellProcessor.h"

// explicit Leapfrog position update rule in a vectorized fashion for the RMM

class PositionCellProcessorRMM : public CellProcessor {
public:
	PositionCellProcessorRMM(double timeStep) : CellProcessor(0.0, 0.0), _timeStep(static_cast<vcp_real_calc>(timeStep)) {}
	void initTraversal() {}

	void preprocessCell(ParticleCell& cell) {}

	void processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll = false) {} // does this need bool?

	double processSingleMolecule(Molecule* m1, ParticleCell& cell2) { return 0.0; }

	void postprocessCell(ParticleCell& cell) {}

	void endTraversal() {}

	void processCell(ParticleCell& cell) {
		ParticleCellRMM & c = downcastCellReferenceRMM(cell);
		CellDataSoARMM & soa = c.getCellDataSoA();
		const size_t end_i = vcp_floor_to_vec_size(soa.getMolNum());

		      vcp_real_calc * const soa_r_x = soa.r_xBegin();
		      vcp_real_calc * const soa_r_y = soa.r_yBegin();
		      vcp_real_calc * const soa_r_z = soa.r_zBegin();
		const vcp_real_calc * const soa_v_x = soa.v_xBegin();
		const vcp_real_calc * const soa_v_y = soa.v_yBegin();
		const vcp_real_calc * const soa_v_z = soa.v_zBegin();

		const RealCalcVec dt = RealCalcVec::set1(_timeStep);

		size_t i = 0;
		for (; i < end_i; i += VCP_VEC_SIZE) {
			RealCalcVec r_x = RealCalcVec::aligned_load(soa_r_x + i);
			RealCalcVec r_y = RealCalcVec::aligned_load(soa_r_y + i);
			RealCalcVec r_z = RealCalcVec::aligned_load(soa_r_z + i);

			const RealCalcVec v_x = RealCalcVec::aligned_load(soa_v_x + i);
			const RealCalcVec v_y = RealCalcVec::aligned_load(soa_v_y + i);
			const RealCalcVec v_z = RealCalcVec::aligned_load(soa_v_z + i);

			r_x = RealCalcVec::fmadd(dt, v_x, r_x);
			r_y = RealCalcVec::fmadd(dt, v_y, r_y);
			r_z = RealCalcVec::fmadd(dt, v_z, r_z);

			r_x.aligned_store(soa_r_x + i);
			r_y.aligned_store(soa_r_y + i);
			r_z.aligned_store(soa_r_z + i);

		}
		const MaskCalcVec remainderMask = vcp_simd_getRemainderMask(soa.getMolNum());
		if (remainderMask.movemask()) {
			RealCalcVec r_x = RealCalcVec::aligned_load_mask(soa_r_x + i, remainderMask);
			RealCalcVec r_y = RealCalcVec::aligned_load_mask(soa_r_y + i, remainderMask);
			RealCalcVec r_z = RealCalcVec::aligned_load_mask(soa_r_z + i, remainderMask);

			const RealCalcVec v_x = RealCalcVec::aligned_load_mask(soa_v_x + i, remainderMask);
			const RealCalcVec v_y = RealCalcVec::aligned_load_mask(soa_v_y + i, remainderMask);
			const RealCalcVec v_z = RealCalcVec::aligned_load_mask(soa_v_z + i, remainderMask);

			r_x = RealCalcVec::fmadd(dt, v_x, r_x);
			r_y = RealCalcVec::fmadd(dt, v_y, r_y);
			r_z = RealCalcVec::fmadd(dt, v_z, r_z);

            // TODO: handle aligned masked store properly in intrinsics.
			// For now just store, not caring about mask. Due to internal padding, this will work without problems and will not overwrite stuff.
			r_x.aligned_store(soa_r_x + i);
			r_y.aligned_store(soa_r_y + i);
			r_z.aligned_store(soa_r_z + i);
		}
	}

private:
	vcp_real_calc _timeStep;
};

#endif /* SRC_INTEGRATORS_POSITIONCELLPROCESSORRMM_H_ */
