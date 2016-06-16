/**
 * \file
 * \brief VectorizedLJP2PCellProcessor.h
 * \author Johannes Heckl, Wolfgang Eckhardt, Uwe Ehmann, Steffen Seckler
 */

#pragma once


#include "particleContainer/adapter/CellProcessor.h"
#include "utils/AlignedArray.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "particleContainer/adapter/vectorization/SIMD_TYPES.h"
#include "particleContainer/adapter/vectorization/SIMD_VectorizedCellProcessorHelpers.h"
#include "utils/Timer.h"

class Component;
class Domain;
class Comp2Param;
class Molecule;
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
	void initTraversal(const size_t numCells);
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
	int countNeighbours(Molecule* /*m1*/, ParticleCell& /*cell2*/, double /*RR*/) {
		exit(0);
		return 0;
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

	void printTimers();

private:
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
	std::vector<DoubleArray> _eps_sig;
	/**
	 * \brief Shift for pairs of LJcenters.
	 * \details Each DoubleArray contains the LJ shift*6.0 for one center combined<br>
	 * with all centers.
	 */
	std::vector<DoubleArray> _shift6;
	/**
	 * \brief Sum of all LJ potentials.
	 * \details Multiplied by 6.0 for performance reasons.
	 */
	double _upot6lj;

	/**
	 * \brief The virial.
	 */
	double _virial;


	/**
	 * \brief array, that stores the dist_lookup.
	 * For all vectorization methods, that utilize masking, this stores masks.
	 * To utilize the gather operations of the MIC architecture, the dist_lookup is able to store the indices of the required particles.
	 */
	AlignedArray<vcp_lookupOrMask_single> _centers_dist_lookup;

	/**
	 * \brief pointer to the starting point of the dist_lookup of the lennard jones particles.
	 */
	vcp_lookupOrMask_single* _ljc_dist_lookup;

	template<bool calculateMacroscopic>
	inline
	void _loopBodyLJ(
			const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
			const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
			const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
			const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
			vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
			vcp_double_vec& V_x, vcp_double_vec& V_y, vcp_double_vec& V_z,
			vcp_double_vec& sum_upot6lj, vcp_double_vec& sum_virial,
			const vcp_mask_vec& forceMask,
			const vcp_double_vec& eps_24, const vcp_double_vec& sig2,
			const vcp_double_vec& shift6);


	/**
	 * \brief The dist lookup for a molecule and all centers of a type
	 * \author Robert Hajda
	 */
	template<class ForcePolicy>
	vcp_mask_vec
	calcDistLookup (const CellDataSoA & soa1, const size_t & i, const size_t & i_center_idx, const size_t & soa2_num_centers, const double & cutoffRadiusSquare,
			vcp_lookupOrMask_single* const soa2_center_dist_lookup, const double* const soa2_m_r_x, const double* const soa2_m_r_y, const double* const soa2_m_r_z,
			const vcp_double_vec & cutoffRadiusSquareD, size_t end_j, const vcp_double_vec m1_r_x, const vcp_double_vec m1_r_y, const vcp_double_vec m1_r_z,
			countertype32& counter
			);



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

	/**
	 * \brief Policy class for single cell force calculation.
	 */
	class SingleCellPolicy_ {
	public:

		inline static bool Condition(double m_r2, double rc2)
		{
			// If m_r2 == 0, it has to be 2 LJ centers in the same molecule
			// (or 2 molecules are at the same location). These are ignored.
			return (m_r2 < rc2) && (m_r2 != 0.0);
		}

		inline static bool DetectSingleCell ()
		{
			return true;
		}

		inline static size_t InitJ (const size_t i)//needed for alignment. (guarantees, that one simd_load always accesses the same cache line.
		{//i+1: only calculate j>i
			// however we do a floor for alignment purposes. ->  we have to mark some of the indices to not be computed (this is handled using the InitJ_Mask)
			return vcp_floor_to_vec_size(i+1); // this is i+1 if i+1 is divisible by VCP_VEC_SIZE otherwise the next smaller multiple of VCP_VEC_SIZE
		}

		inline static size_t InitJ2 (const size_t i)//needed for alignment. (guarantees, that one simd_load always accesses the same cache line.
		{
#if VCP_VEC_TYPE!=VCP_VEC_MIC_GATHER
			return InitJ(i);
#else
			return 0;
#endif
		}

#if VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_VEC_AVX2 or VCP_VEC_TYPE==VCP_VEC_MIC or VCP_VEC_TYPE==VCP_VEC_MIC_GATHER
		inline static vcp_mask_vec GetForceMask (const vcp_double_vec& m_r2, const vcp_double_vec& rc2, vcp_mask_vec& j_mask)
		{
			vcp_mask_vec result = vcp_simd_and( vcp_simd_and(vcp_simd_lt(m_r2, rc2), vcp_simd_neq(m_r2, VCP_SIMD_ZEROV) ), j_mask);
			j_mask = VCP_SIMD_ONESVM;
			return result;
		}
	#if VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_VEC_AVX2
		inline static vcp_mask_vec InitJ_Mask (const size_t i)//calculations only for i+1 onwards.
		{
			switch (i & static_cast<size_t>(VCP_VEC_SIZE_M1)) {
				case 0: return _mm256_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, 0, 0);
				case 1: return _mm256_set_epi32(~0, ~0, ~0, ~0, 0, 0, 0, 0);
				case 2: return _mm256_set_epi32(~0, ~0, 0, 0, 0, 0, 0, 0);
				default: return VCP_SIMD_ONESVM;
			}
		}
	#else //mic
		inline static vcp_mask_vec InitJ_Mask (const size_t i)//calculations only for i+1 onwards.
		{
			static const vcp_mask_vec possibleInitJMasks[VCP_VEC_SIZE] = { 0xFE, 0xFC, 0xF8, 0xF0, 0xE0, 0xC0, 0x80, 0xFF };
			return possibleInitJMasks[i & static_cast<size_t>(VCP_VEC_SIZE_M1)];
		}
	#endif
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
		// Erstellen der Bitmaske, analog zu Condition oben
		inline static vcp_mask_vec GetForceMask(vcp_double_vec m_r2, vcp_double_vec rc2)
		{
			return vcp_simd_and(vcp_simd_lt(m_r2, rc2), vcp_simd_neq(m_r2, VCP_SIMD_ZEROV) );
		}
#elif VCP_VEC_TYPE==VCP_NOVEC
		inline static vcp_mask_vec GetForceMask(vcp_double_vec m_r2, vcp_double_vec rc2)
		{
			return vcp_simd_lt(m_r2, rc2);
		}
#endif /* definition of InitJ_Mask & GetForceMask */
	}; /* end of class SingleCellPolicy_ */

	/**
	 * \brief Policy class for cell pair force calculation.
	 */
	class CellPairPolicy_ {
	public:

		inline static bool Condition(double m_r2, double rc2)
		{
			// Because we have 2 different cells, no 2 centers on the same
			// molecule can form a pair.
			return m_r2 < rc2;
		}

		inline static size_t InitJ (const size_t /*i*/)
		{
			return 0;
		}
		inline static size_t InitJ2 (const size_t i)//needed for alignment. (guarantees, that one simd_load always accesses the same cache line.
		{
			return InitJ(i);
		}


		inline static bool DetectSingleCell ()
		{
			return false;
		}

		inline static vcp_mask_vec GetForceMask (const vcp_double_vec& m_r2, const vcp_double_vec& rc2
#if VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_VEC_AVX2 or VCP_VEC_TYPE==VCP_VEC_MIC or VCP_VEC_TYPE==VCP_VEC_MIC_GATHER
				, vcp_mask_vec& /*j_mask*/
#endif
				)
		{
			// Provide a mask with the same logic as used in
			// bool Condition(double m_r2, double rc2)
			return vcp_simd_lt(m_r2, rc2);
		}


#if VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_VEC_AVX2 or VCP_VEC_TYPE==VCP_VEC_MIC or VCP_VEC_TYPE==VCP_VEC_MIC_GATHER
		inline static vcp_mask_vec InitJ_Mask (const size_t /*i*/)
		{
			return VCP_SIMD_ZEROVM;//totally unimportant, since not used...
		}
#endif
 /* definition of GetForceMask */
	}; /* end of class CellPairPolicy_ */

}; /* end of class VectorizedLJP2PCellProcessor */

} // namespace bhfmm
