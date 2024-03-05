/*
 * VelocityCellProcessorRMM.h
 *
 *  Created on: 5 Oct 2017
 *      Author: tchipevn
 */

#ifndef SRC_INTEGRATORS_VELOCITYCELLPROCESSORRMM_H_
#define SRC_INTEGRATORS_VELOCITYCELLPROCESSORRMM_H_

#include "WrapOpenMP.h"
#include "particleContainer/ParticleCell.h"
#include "particleContainer/adapter/CellProcessor.h"
#include "particleContainer/adapter/vectorization/SIMD_VectorizedCellProcessorHelpers.h"

// just compute summv2 in a vectorized fashion in RMM

class VelocityCellProcessorRMM: public CellProcessor {
public:
	VelocityCellProcessorRMM() :
			CellProcessor(0.0, 0.0), _N(0), _summv2(0.0) {

		_threadData.resize(mardyn_get_max_threads());

		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			ThreadData * myown = new ThreadData();
			const int myid = mardyn_get_thread_num();
			_threadData[myid] = myown;
		} // end pragma omp parallel
	}
	~VelocityCellProcessorRMM() {
		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			const int myid = mardyn_get_thread_num();
			delete _threadData[myid];
		}
	}

	void preprocessCell(ParticleCell& cell) {}

	void processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll = false) {}

	double processSingleMolecule(Molecule* m1, ParticleCell& cell2) { return 0.0; }

	void postprocessCell(ParticleCell& cell) {}

	void initTraversal() {
		#if defined(_OPENMP)
		#pragma omp master
		#endif
		{
			_N = 0;
			_summv2 = 0.0;
		} // end pragma omp master

		Log::global_log->debug() << "VelocityCellProcessorRMM::initTraversal()." << std::endl;
	}

	void endTraversal() {
		vcp_real_accum glob_summv2 = 0.0;
		unsigned long glob_N = 0;

		#if defined(_OPENMP)
		#pragma omp parallel reduction(+:glob_summv2, glob_N)
		#endif
		{
			const int tid = mardyn_get_thread_num();

			// reduce vectors and clear local variable
			vcp_real_accum thread_summv2 = 0.0;

			load_hSum_Store_Clear(&thread_summv2, _threadData[tid]->_thread_summv2V);

			// add to global sum
			glob_summv2 += thread_summv2;
			glob_N += _threadData[tid]->_thread_N;
			_threadData[tid]->_thread_N = 0;
		} // end pragma omp parallel reduction

		_summv2 = glob_summv2;
		_N = glob_N;
	}

	void processCell(ParticleCell& cell) {
		// just compute the velocity sums

		ParticleCellRMM & c = downcastCellReferenceRMM(cell);
		CellDataSoARMM & soa = c.getCellDataSoA();

		const size_t molNum = soa.getMolNum();
		const size_t end_i = vcp_floor_to_vec_size(molNum);

		const int tid = mardyn_get_thread_num();
		ThreadData &my_threadData = *_threadData[tid];
		my_threadData._thread_N += static_cast<unsigned long>(molNum);

		const vcp_real_accum * const soa_v_x = soa.v_xBegin();
		const vcp_real_accum * const soa_v_y = soa.v_yBegin();
		const vcp_real_accum * const soa_v_z = soa.v_zBegin();

		RealAccumVec sum_summv2 = RealAccumVec::zero();

		size_t i = 0;
		for (; i < end_i; i += VCP_VEC_SIZE) {
			const RealAccumVec v_x = RealAccumVec::aligned_load(soa_v_x + i);
			const RealAccumVec v_y = RealAccumVec::aligned_load(soa_v_y + i);
			const RealAccumVec v_z = RealAccumVec::aligned_load(soa_v_z + i);

			const RealAccumVec v2 = RealAccumVec::scal_prod(v_x, v_y, v_z, v_x, v_y, v_z);
			sum_summv2 = sum_summv2 + v2;

		}
		const MaskCalcVec remainderMask = vcp_simd_getRemainderMask(soa.getMolNum());
		if (remainderMask.movemask()) {
			const RealAccumVec v_x = RealAccumVec::aligned_load_mask(soa_v_x + i, remainderMask);
			const RealAccumVec v_y = RealAccumVec::aligned_load_mask(soa_v_y + i, remainderMask);
			const RealAccumVec v_z = RealAccumVec::aligned_load_mask(soa_v_z + i, remainderMask);

			const RealAccumVec v2 = RealAccumVec::scal_prod(v_x, v_y, v_z, v_x, v_y, v_z);
			sum_summv2 = sum_summv2 + v2;
		}
		sum_summv2.aligned_load_add_store(&(my_threadData._thread_summv2V[0]));
	}

	class ThreadData {
	public:
		ThreadData() {
			_thread_summv2V.resize(VCP_VEC_SIZE);

			for (size_t j = 0; j < VCP_VEC_SIZE; ++j) {
				_thread_summv2V[j] = 0.0;
			}
			_thread_N = 0;
		}

		AlignedArray<vcp_real_accum> _thread_summv2V;
		unsigned long _thread_N;
	};

	unsigned long getN() const {
		return _N;
	}

	double getSummv2() const {
		return _summv2;
	}
private:
	unsigned long _N;
	double _summv2;
	std::vector<ThreadData *> _threadData;
};

#endif /* SRC_INTEGRATORS_VELOCITYCELLPROCESSORRMM_H_ */
