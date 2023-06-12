/**
 * \file
 * \brief VectorizedLJP2PCellProcessor.h
 * \author Johannes Heckl, Wolfgang Eckhardt, Uwe Ehmann, Steffen Seckler
 */

#pragma once

#include <stdlib.h>
#include "particleContainer/adapter/CellProcessor.h"
#include "utils/AlignedArray.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "particleContainer/adapter/vectorization/SIMD_TYPES.h"
#include "particleContainer/adapter/vectorization/SIMD_VectorizedCellProcessorHelpers.h"
#include "WrapOpenMP.h"

#include "molecules/MoleculeForwardDeclaration.h"

class Component;
class Domain;
class Comp2Param;
class CellDataSoA;

namespace bhfmm {
/**
 * \brief Vectorized calculation of the force.
 * \author Johannes Heckl
 */
class VectorizedLJP2PCellProcessor : public CellProcessor {
public:
	typedef std::vector<Component> ComponentList;

	/**
	 * \brief Construct and set up the internal parameter table.
	 * \details Components and parameters should be finalized before this call.
	 */
	VectorizedLJP2PCellProcessor(Domain & domain, double cutoffRadius, double LJcutoffRadius);

	~VectorizedLJP2PCellProcessor();

	/**
	 * \brief Reset macroscopic values to 0.0.
	 */
	void initTraversal();
	/**
	 * \brief Load the CellDataSoA for cell.
	 */
	void preprocessCell(ParticleCell& /*cell*/) {}

	double processSingleMolecule(Molecule* /*m1*/, ParticleCell& /*cell2*/) {
		return 0.0;
	}

	/**
	 * \brief Calculate forces between pairs of Molecules in cell.
	 */
	void processCell(ParticleCell& cell);

        void processCellPair(ParticleCell& c1, ParticleCell& c2, bool sumAll = false);
	/**
	 * \brief Free the LennardJonesSoA for cell.
	 */
	void postprocessCell(ParticleCell& /*cell*/) {}
	/**
	 * \brief Store macroscopic values in the Domain.
	 */
	void endTraversal();

	void printTimers();


private:
	/**
	 * \brief a vector of Molecule pointers.
	 */
	typedef std::vector<Molecule *> MoleculeList;

	/**
	 * \brief The Domain where macroscopic values will be stored.
	 */
	Domain & _domain;

	/**
	 * \brief The squared cutoff radius.
	 */
	//const double _cutoffRadiusSquare;

	/**
	 * \brief The squared LJ cutoff radius.
	 */
	//const double _LJcutoffRadiusSquare;

	/**
	 * \brief Parameter for the reaction field method (see description in Domain.h and Comp2Param.cpp).
	 */
	const double _epsRFInvrc3;

	/**
	 * \brief One LJ center enumeration start index for each component.
	 * \details All the LJ centers of all components are enumerated.<br>
	 * Comp1 gets indices 0 through n1 - 1, Comp2 n1 through n2 - 1 and so on.<br>
	 * This is necessary for finding the respective parameters for each interaction<br>
	 * between two centers.
	 */

	/**
	 * \brief Epsilon and sigma for pairs of LJcenters.
	 * \details Each DoubleArray contains parameters for one center combined with all centers.<br>
	 * Each set of parameters is a pair (epsilon*24.0, sigma^2).
	 */
	std::vector<AlignedArray<vcp_real_calc> > _eps_sig;

	/**
	 * \brief Shift for pairs of LJcenters.
	 * \details Each DoubleArray contains the LJ shift*6.0 for one center combined<br>
	 * with all centers.
	 */
	std::vector<AlignedArray<vcp_real_calc> > _shift6;

	/**
	 * \brief Sum of all LJ potentials.
	 * \details Multiplied by 6.0 for performance reasons.
	 */
	double _upot6lj;

	/**
	 * \brief The virial.
	 */
	double _virial;

	struct VLJP2PCPThreadData {
	public:
		VLJP2PCPThreadData(): _ljc_dist_lookup(nullptr){
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
		 * To utilize the gather operations of the KNL architecture, the dist_lookup is able to store the indices of the required particles.
		 */
		AlignedArray<vcp_lookupOrMask_single> _centers_dist_lookup;

		/**
		 * \brief pointer to the starting point of the dist_lookup of the lennard jones particles.
		 */
		vcp_lookupOrMask_single* _ljc_dist_lookup;

		AlignedArray<vcp_real_accum> _upot6ljV, _virialV;
	};

	std::vector<VLJP2PCPThreadData *> _threadData;

	static const size_t _numVectorElements = VCP_VEC_SIZE;
	size_t _numThreads;

	template<bool calculateMacroscopic>
	inline
	void _loopBodyLJ(
			const RealCalcVec& m1_r_x, const RealCalcVec& m1_r_y, const RealCalcVec& m1_r_z,
			const RealCalcVec& r1_x, const RealCalcVec& r1_y, const RealCalcVec& r1_z,
			const RealCalcVec& m2_r_x, const RealCalcVec& m2_r_y, const RealCalcVec& m2_r_z,
			const RealCalcVec& r2_x, const RealCalcVec& r2_y, const RealCalcVec& r2_z,
			RealCalcVec& f_x, RealCalcVec& f_y, RealCalcVec& f_z,
			RealAccumVec& V_x, RealAccumVec& V_y, RealAccumVec& V_z,
			RealAccumVec& sum_upot6lj, RealAccumVec& sum_virial,
			const MaskCalcVec& forceMask,
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
	void _calculatePairs(CellDataSoA & soa1, CellDataSoA & soa2);

}; /* end of class VectorizedLJP2PCellProcessor */

} // namespace bhfmm
