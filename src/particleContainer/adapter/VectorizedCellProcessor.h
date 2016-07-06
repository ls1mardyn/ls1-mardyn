/**
 * \file
 * \brief VectorizedCellProcessor.h
 * \author Johannes Heckl, Wolfgang Eckhardt, Uwe Ehmann, Steffen Seckler
 */

#ifndef VECTORIZEDCELLPROCESSOR_H_
#define VECTORIZEDCELLPROCESSOR_H_


#include "CellProcessor.h"
#include "utils/AlignedArray.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "vectorization/SIMD_TYPES.h"
#include "vectorization/SIMD_VectorizedCellProcessorHelpers.h"


class Component;
class Domain;
class Comp2Param;
class Molecule;
class CellDataSoA;

/**
 * \brief Vectorized calculation of the force.
 * \author Johannes Heckl
 */
class VectorizedCellProcessor : public CellProcessor {
public:
	typedef std::vector<Component> ComponentList;

	VectorizedCellProcessor& operator=(const VectorizedCellProcessor&) = delete;

	/**
	 * \brief Construct and set up the internal parameter table.
	 * \details Components and parameters should be finalized before this call.
	 */
	VectorizedCellProcessor(Domain & domain, double cutoffRadius, double LJcutoffRadius);

	~VectorizedCellProcessor();

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

	double processSingleMolecule(Molecule* /*m1*/, ParticleCell& /*cell2*/) { return 0.0; }

	// provisionally, the code from the legacy cell processor is used here
	//
        int countNeighbours(Molecule* m1, ParticleCell& cell2, double RR);

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
private:
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
	 * \brief Sum of all Xpole potentials.
	 */
	double _upotXpoles;

	/**
	 * \brief The virial.
	 */
	double _virial;

	/**
	 * \brief MyRF contribution of all pairs
	 */
	double _myRF;

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

	/**
	 * \brief pointer to the starting point of the dist_lookup of the charge particles.
	 */
	vcp_lookupOrMask_single* _charges_dist_lookup;

	/**
	 * \brief pointer to the starting point of the dist_lookup of the dipole particles.
	 */
	vcp_lookupOrMask_single* _dipoles_dist_lookup;

	/**
	 * \brief pointer to the starting point of the dist_lookup of the quadrupole particles.
	 */
	vcp_lookupOrMask_single* _quadrupoles_dist_lookup;

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

	template<bool calculateMacroscopic>
	inline void _loopBodyChargeDipole(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& q,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& e_x, const vcp_double_vec& e_y, const vcp_double_vec& e_z,
		const vcp_double_vec& p,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& V_x, vcp_double_vec& V_y, vcp_double_vec& V_z,
		vcp_double_vec& M_x, vcp_double_vec& M_y, vcp_double_vec& M_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_mask_vec& forceMask);

	template<bool calculateMacroscopic>
	inline void _loopBodyDipole(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& eii_x, const vcp_double_vec& eii_y, const vcp_double_vec& eii_z,
		const vcp_double_vec& pii,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& ejj_x, const vcp_double_vec& ejj_y, const vcp_double_vec& ejj_z,
		const vcp_double_vec& pjj,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& V_x, vcp_double_vec& V_y, vcp_double_vec& V_z,
		vcp_double_vec& M1_x, vcp_double_vec& M1_y, vcp_double_vec& M1_z,
		vcp_double_vec& M2_x, vcp_double_vec& M2_y, vcp_double_vec& M2_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial, vcp_double_vec& sum_myRF,
		const vcp_mask_vec& forceMask,
		const vcp_double_vec& epsRFInvrc3);

	template<bool calculateMacroscopic>
	inline void _loopBodyChargeQuadrupole(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& q,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& ejj_x, const vcp_double_vec& ejj_y, const vcp_double_vec& ejj_z,
		const vcp_double_vec& m,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& V_x, vcp_double_vec& V_y, vcp_double_vec& V_z,
		vcp_double_vec& M_x, vcp_double_vec& M_y, vcp_double_vec& M_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_mask_vec& forceMask);

	template<bool calculateMacroscopic>
	inline void _loopBodyDipoleQuadrupole(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& eii_x, const vcp_double_vec& eii_y, const vcp_double_vec& eii_z,
		const vcp_double_vec& p,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& ejj_x, const vcp_double_vec& ejj_y, const vcp_double_vec& ejj_z,
		const vcp_double_vec& m,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& V_x, vcp_double_vec& V_y, vcp_double_vec& V_z,
		vcp_double_vec& M1_x, vcp_double_vec& M1_y, vcp_double_vec& M1_z,
		vcp_double_vec& M2_x, vcp_double_vec& M2_y, vcp_double_vec& M2_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_mask_vec& forceMask);

	template<bool calculateMacroscopic>
	inline void _loopBodyQuadrupole(
		const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
		const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
		const vcp_double_vec& eii_x, const vcp_double_vec& eii_y, const vcp_double_vec& eii_z,
		const vcp_double_vec& mii,
		const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
		const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
		const vcp_double_vec& ejj_x, const vcp_double_vec& ejj_y, const vcp_double_vec& ejj_z,
		const vcp_double_vec& mjj,
		vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
		vcp_double_vec& V_x, vcp_double_vec& V_y, vcp_double_vec& V_z,
		vcp_double_vec& Mii_x, vcp_double_vec& Mii_y, vcp_double_vec& Mii_z,
		vcp_double_vec& Mjj_x, vcp_double_vec& Mjj_y, vcp_double_vec& Mjj_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_mask_vec& forceMask);


	/**
	 * \brief The dist lookup for a molecule and all centers of a type
	 * \author Robert Hajda
	 */
	template<class ForcePolicy>
	vcp_mask_vec
	calcDistLookup (const size_t & i_center_idx, const size_t & soa2_num_centers, const double & cutoffRadiusSquare,
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

		inline static size_t InitJ (const size_t i)//needed for alignment. (guarantees, that one simd_load always accesses the same cache line.
		{//i: only calculate j>=i
			// however we do a floor for alignment purposes. ->  we have to mark some of the indices to not be computed (this is handled using the InitJ_Mask)
			return vcp_floor_to_vec_size(i); // this is i if i is divisible by VCP_VEC_SIZE otherwise the next smaller multiple of VCP_VEC_SIZE
		}

		inline static size_t InitJ2 (const size_t i __attribute__((unused)))//needed for alignment. (guarantees, that one simd_load always accesses the same cache line.
		{
#if VCP_VEC_TYPE!=VCP_VEC_MIC_GATHER
			return InitJ(i);
#else
			return 0;
#endif
		}


		inline static vcp_mask_vec GetForceMask (const vcp_double_vec& m_r2, const vcp_double_vec& rc2, vcp_mask_vec& j_mask)
		{
			vcp_mask_vec result = vcp_simd_and( vcp_simd_and(vcp_simd_lt(m_r2, rc2), vcp_simd_neq(m_r2, VCP_SIMD_ZEROV) ), j_mask);
			j_mask = VCP_SIMD_ONESVM;
			return result;
		}
	#if VCP_VEC_TYPE==VCP_VEC_SSE3
		inline static vcp_mask_vec InitJ_Mask (const size_t i)//calculations only for i onwards.
		{
			switch (i & static_cast<size_t>(VCP_VEC_SIZE_M1)) {
				case 0: return _mm_set_epi32(~0, ~0, ~0, ~0);
				default: return _mm_set_epi32(~0, ~0, 0, 0);
			}
		}
	#elif VCP_VEC_TYPE==VCP_VEC_AVX or VCP_VEC_TYPE==VCP_VEC_AVX2
		inline static vcp_mask_vec InitJ_Mask (const size_t i)//calculations only for i onwards.
		{
			switch (i & static_cast<size_t>(VCP_VEC_SIZE_M1)) {
				case 0: return _mm256_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0);
				case 1: return _mm256_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, 0, 0);
				case 2: return _mm256_set_epi32(~0, ~0, ~0, ~0, 0, 0, 0, 0);
				default: return _mm256_set_epi32(~0, ~0, 0, 0, 0, 0, 0, 0);
			}
		}
	#elif VCP_VEC_TYPE==VCP_NOVEC
		inline static vcp_mask_vec InitJ_Mask (const size_t /*i*/)//calculations only for i onwards.
			{
				return true; // calculate everything
			}
	#else //mic
		inline static vcp_mask_vec InitJ_Mask (const size_t i)//calculations only for i onwards.
		{
			static const vcp_mask_vec possibleInitJMasks[VCP_VEC_SIZE] = { 0xFF, 0xFE, 0xFC, 0xF8, 0xF0, 0xE0, 0xC0, 0x80 };
			return possibleInitJMasks[i & static_cast<size_t>(VCP_VEC_SIZE_M1)];
		}
	#endif

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

		inline static vcp_mask_vec GetForceMask (const vcp_double_vec& m_r2, const vcp_double_vec& rc2, vcp_mask_vec& /*j_mask*/)
		{
			// Provide a mask with the same logic as used in
			// bool Condition(double m_r2, double rc2)
			return vcp_simd_lt(m_r2, rc2);
		}

		inline static vcp_mask_vec InitJ_Mask (const size_t /*i*/)
		{
			return VCP_SIMD_ONESVM;//totally unimportant, since not used...
		}

	}; /* end of class CellPairPolicy_ */

}; /* end of class VectorizedCellProcessor */

#endif /* VECTORIZEDCELLPROCESSOR_H_ */
