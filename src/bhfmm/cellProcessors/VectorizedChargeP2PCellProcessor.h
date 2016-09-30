/**
 * \file
 * \brief VectorizedChargeP2PCellProcessor.h
 * \author Johannes Heckl, Wolfgang Eckhardt, Uwe Ehmann, Steffen Seckler
 */

#pragma once

#include <stdlib.h>
#include "bhfmm/containers/ParticleCellPointers.h"
#include "utils/AlignedArray.h"
#include "utils/Timer.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "particleContainer/adapter/vectorization/SIMD_TYPES.h"
#include "particleContainer/adapter/vectorization/SIMD_VectorizedCellProcessorHelpers.h"
#include "particleContainer/adapter/CellDataSoA.h"
#include "WrapOpenMP.h"

class Component;
class Domain;
class Comp2Param;
class Molecule;
//class CellDataSoA;
namespace bhfmm {
/**
 * \brief Vectorized calculation of the force.
 * \author Johannes Heckl
 */
class VectorizedChargeP2PCellProcessor {
public:
	typedef std::vector<Component> ComponentList;

	/**
	 * \brief Construct and set up the internal parameter table.
	 * \details Components and parameters should be finalized before this call.
	 */
	VectorizedChargeP2PCellProcessor(Domain & domain, double cutoffRadius=0, double LJcutoffRadius=0);

	~VectorizedChargeP2PCellProcessor();

	/**
	 * \brief Reset macroscopic values to 0.0.
	 */
	void initTraversal();
	/**
	 * \brief Load the CellDataSoA for cell.
	 */
	void preprocessCell(ParticleCellPointers& cell);
	/**
	 * \brief Calculate forces between pairs of Molecules in cell1 and cell2.
	 */
	void processCellPair(ParticleCellPointers& cell1, ParticleCellPointers& cell2);
	/**
	 * \brief Calculate forces between pairs of Molecules in cell.
	 */
	void processCell(ParticleCellPointers& cell);
	/**
	 * \brief Free the LennardJonesSoA for cell.
	 */
	void postprocessCell(ParticleCellPointers& cell);
	/**
	 * \brief Store macroscopic values in the Domain.
	 */
	void endTraversal();

	void printTimers();

private:
	double _cutoffRadiusSquare;

	Timer _timer;

	/**
	 * \brief An aligned array of doubles.
	 */
	typedef AlignedArray<double> DoubleArray;

	/**
	 * \brief a vector of Molecule pointers.
	 */
	typedef std::vector<Molecule *> MoleculeList;

	/**
	 * \brief The Domain where macroscopic values will be stored.
	 */
	Domain & _domain;

	/**
	 * \brief Sum of all Xpole potentials.
	 */
	double _upotXpoles;

	/**
	 * \brief The virial.
	 */
	double _virial;

	struct VCP2PCPThreadData {
	public:
		VCP2PCPThreadData(): _charges_dist_lookup(nullptr){
			_upotXpolesV.resize(_numVectorElements);
			_virialV.resize(_numVectorElements);

			for (size_t j = 0; j < _numVectorElements; ++j) {
				_upotXpolesV[j] = 0.0;
				_virialV[j] = 0.0;
			}
		}

		/**
		 * \brief array, that stores the dist_lookup.
		 * For all vectorization methods, that utilize masking, this stores masks.
		 * To utilize the gather operations of the MIC architecture, the dist_lookup is able to store the indices of the required particles.
		 */
		AlignedArray<vcp_lookupOrMask_single> _centers_dist_lookup;

		/**
		 * \brief pointer to the starting point of the dist_lookup of the charge particles.
		 */
		vcp_lookupOrMask_single* _charges_dist_lookup;

		AlignedArray<double> _upotXpolesV, _virialV;
	};

	std::vector<VCP2PCPThreadData *> _threadData;

	static const size_t _numVectorElements = VCP_VEC_SIZE;
	size_t _numThreads;

	template<bool calculateMacroscopic>
	inline void _loopBodyCharge(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& qii,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& qjj,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& V_x, vcp_double_vec& V_y, vcp_double_vec& V_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_mask_vec& forceMask);

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
	 * static vcp_double_vec GetForceMask(vcp_double_vec m_r2, vcp_double_vec rc2);<br>
	 * Returns the mask indicating which pairs to calculate in the vectorized code.<br>
	 * <br>
	 * The boolean CalculateMacroscopic should specify, whether macroscopic values are to be calculated or not.
	 * <br>
	 * The class MaskGatherChooser is a class, that specifies the used loading,storing and masking routines.
	 */
	template<class ForcePolicy, bool CalculateMacroscopic, class MaskGatherChooser>
	void _calculatePairs(const CellDataSoA & soa1, const CellDataSoA & soa2);

}; /* end of class VectorizedChargeP2PCellProcessor */

} // namespace bhfmm

