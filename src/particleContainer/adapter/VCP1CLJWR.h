/*
 * VCP1CLJWR.h
 *
 *  Created on: 30 Jan 2017
 *      Author: tchipevn
 */

#ifndef SRC_PARTICLECONTAINER_ADAPTER_VCP1CLJWR_H_
#define SRC_PARTICLECONTAINER_ADAPTER_VCP1CLJWR_H_

#include "CellProcessor.h"

#include "particleContainer/adapter/vectorization/SIMD_TYPES.h"
#include "particleContainer/adapter/vectorization/SIMD_VectorizedCellProcessorHelpers.h"
#include "WrapOpenMP.h"

#include "molecules/MoleculeForwardDeclaration.h"

#include <vector>

class Component;
class Domain;
class Comp2Param;
class CellDataSoA_WR;

class VCP1CLJ_WR: public CellProcessor {
public:
	VCP1CLJ_WR(Domain & domain, double cutoffRadius, double LJcutoffRadius);
	~VCP1CLJ_WR();


	/**
	 * \brief Reset macroscopic values to 0.0.
	 */
	void initTraversal();
	/**
	 * \brief Load the CellDataSoA for cell.
	 */
	void preprocessCell(ParticleCell& /*cell*/) {}
	/**
	 * \brief Calculate forces between pairs of Molecules in cell1 and cell2.
	 */
	void processCellPair(ParticleCell& cell1, ParticleCell& cell2);

	double processSingleMolecule(Molecule* /*m1*/, ParticleCell& /*cell2*/) {
		return 0.0;
	}

	/**
	 * \brief Calculate forces between pairs of Molecules in cell.
	 */
	void processCell(ParticleCell& cell);
	/**
	 * \brief Free the LennardJonesSoA for cell.
	 */
	void postprocessCell(ParticleCell& /*cell*/) {}
	/**
	 * \brief Store macroscopic values in the Domain.
	 */
	void endTraversal();

	void setDtInvm(double dtInvm) {
		_dtInvm = static_cast<vcp_real_calc>(dtInvm);
	}

	double getDtInvm() const {
		return _dtInvm;
	}

private:
	/**
	 * \brief The Domain where macroscopic values will be stored.
	 */
	Domain & _domain;


	vcp_real_calc _eps24, _sig2, _shift6, _dtInvm;

	/**
	 * \brief Sum of all LJ potentials.
	 * \details Multiplied by 6.0 for performance reasons.
	 */
	double _upot6lj;

	/**
	 * \brief The virial.
	 */
	double _virial;

	struct VCP1CLJWRThreadData {
	public:
		VCP1CLJWRThreadData(): _ljc_dist_lookup(nullptr){
			_upot6ljV.resize(_numVectorElements);
			_virialV.resize(_numVectorElements);

			for (size_t j = 0; j < _numVectorElements; ++j) {
				_upot6ljV[j] = 0.0;
				_virialV[j] = 0.0;
			}
		}

		/**
		 * \brief array, that stores the dist_lookup.
		 * For all vectorization methods, that utilize masking, this stores masks.
		 * To utilize the gather operations of the KNC architecture, the dist_lookup is able to store the indices of the required particles.
		 */
		AlignedArray<vcp_lookupOrMask_single> _centers_dist_lookup;

		/**
		 * \brief pointer to the starting point of the dist_lookup of the lennard jones particles.
		 */
		vcp_lookupOrMask_single* _ljc_dist_lookup;

		AlignedArray<vcp_real_calc> _upot6ljV, _virialV;
	};

	std::vector<VCP1CLJWRThreadData *> _threadData;

	static const size_t _numVectorElements = VCP_VEC_SIZE;
	size_t _numThreads;

	template<bool calculateMacroscopic>
	inline
	void _loopBodyLJ(
		const RealCalcVec& c_dx, const RealCalcVec& c_dy, const RealCalcVec& c_dz, const RealCalcVec& c_r2,
		RealCalcVec& f_x, RealCalcVec& f_y, RealCalcVec& f_z,
		RealCalcVec& sum_upot6lj, RealCalcVec& sum_virial,
		const MaskVec& forceMask,
		const RealCalcVec& eps_24, const RealCalcVec& sig2,
		const RealCalcVec& shift6);

	/**
	 * \brief Force calculation with abstraction of cell pairs.
	 * \details The differences between single cell and cell pair calculation<br>
	 * have been moved into two policy class templates.<br>
	 * <br>
	 * The ForcePolicy class must provide the following methods:<br>
	 * <br>
	 * static size_t InitJ(size_t i);<br>
	 * Returns the value which j is to be initialized as in the inner loop<br>
	 * depending on the vectorization method.<br>
	 * <br>
	 * static bool Condition(double m_r2, double rc2);<br>
	 * Returns whether to calculate the force for a non-vectorized pair.<br>
	 * <br>
	 * If the code is to be vectorized:<br>
	 * static DoubleVec GetForceMask(DoubleVec m_r2, DoubleVec rc2);<br>
	 * Returns the mask indicating which pairs to calculate in the vectorized code.<br>
	 * <br>
	 * The boolean CalculateMacroscopic should specify, whether macroscopic values are to be calculated or not.
	 * <br>
	 * The class MaskGatherChooser is a class, that specifies the used loading,storing and masking routines.
	 */
	template<class ForcePolicy, bool CalculateMacroscopic, class MaskGatherChooser>
	void _calculatePairs(CellDataSoA_WR & soa1, CellDataSoA_WR & soa2);

};

#endif /* SRC_PARTICLECONTAINER_ADAPTER_VCP1CLJWR_H_ */
