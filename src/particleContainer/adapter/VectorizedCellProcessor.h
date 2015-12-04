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

	/**
	 * \brief Construct and set up the internal parameter table.
	 * \details Components and parameters should be finalized before this call.
	 */
	VectorizedCellProcessor(Domain & domain, double cutoffRadius, double LJcutoffRadius);

	~VectorizedCellProcessor();

	/**
	 * \brief Reset macroscopic values to 0.0.
	 */
	void initTraversal(const size_t numCells);
	/**
	 * \brief Load the CellDataSoA for cell.
	 */
	void preprocessCell(ParticleCell& cell);
	/**
	 * \brief Calculate forces between pairs of Molecules in cell1 and cell2.
	 */
	void processCellPair(ParticleCell& cell1, ParticleCell& cell2);

	double processSingleMolecule(Molecule* m1, ParticleCell& cell2) { return 0.0; }

	/**
	 * \brief Calculate forces between pairs of Molecules in cell.
	 */
	void processCell(ParticleCell& cell);
	/**
	 * \brief Free the LennardJonesSoA for cell.
	 */
	void postprocessCell(ParticleCell& cell);
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
	 * \brief One LJ center enumeration start index for each component.
	 * \details All the LJ centers of all components are enumerated.<br>
	 * Comp1 gets indices 0 through n1 - 1, Comp2 n1 through n2 - 1 and so on.<br>
	 * This is necessary for finding the respective parameters for each interaction<br>
	 * between two centers.
	 */
	std::vector<size_t> _compIDs;
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

	// vector holding pointers to the above objects, used as a stack for
	// managing free objects
	std::vector<CellDataSoA*> _particleCellDataVector;

	template<bool calculateMacroscopic>
	inline
	void _loopBodyLJ(
			const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
			const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
			const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
			const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
			vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
			vcp_double_vec& sum_upot6lj, vcp_double_vec& sum_virial,
			const vcp_double_vec& forceMask,
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
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_double_vec& forceMask);

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
		vcp_double_vec& M_x, vcp_double_vec& M_y, vcp_double_vec& M_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_double_vec& forceMask, const vcp_double_vec& switched);

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
		vcp_double_vec& M1_x, vcp_double_vec& M1_y, vcp_double_vec& M1_z,
		vcp_double_vec& M2_x, vcp_double_vec& M2_y, vcp_double_vec& M2_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial, vcp_double_vec& sum_myRF,
		const vcp_double_vec& forceMask,
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
		vcp_double_vec& M_x, vcp_double_vec& M_y, vcp_double_vec& M_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_double_vec& forceMask, const vcp_double_vec& switched);

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
		vcp_double_vec& M1_x, vcp_double_vec& M1_y, vcp_double_vec& M1_z,
		vcp_double_vec& M2_x, vcp_double_vec& M2_y, vcp_double_vec& M2_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_double_vec& forceMask, const vcp_double_vec& switched);

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
		vcp_double_vec& Mii_x, vcp_double_vec& Mii_y, vcp_double_vec& Mii_z,
		vcp_double_vec& Mjj_x, vcp_double_vec& Mjj_y, vcp_double_vec& Mjj_z,
		vcp_double_vec& sum_upotXpoles, vcp_double_vec& sum_virial,
		const vcp_double_vec& forceMask);


	/**
	 * \brief The dist lookup for a molecule and all centers of a type
	 * \author Robert Hajda
	 */
	template<class ForcePolicy>
	vcp_doublesizedmask_vec
	calcDistLookup (const CellDataSoA & soa1, const size_t & i, const size_t & i_center_idx, const size_t & soa2_num_centers, const double & cutoffRadiusSquare,
			double* const soa2_center_dist_lookup, const double* const soa2_m_r_x, const double* const soa2_m_r_y, const double* const soa2_m_r_z
			, const vcp_double_vec & cutoffRadiusSquareD, size_t end_j, const vcp_double_vec m1_r_x, const vcp_double_vec m1_r_y, const vcp_double_vec m1_r_z
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
	 *
	 */
	template<class ForcePolicy, bool CalculateMacroscopic>
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
		{
			return vcp_floor_to_vec_size(i+1); // this is i+1 if i+1 is divisible by VCP_VEC_SIZE otherwise the next smaller multiple of VCP_VEC_SIZE
		}

#if VCP_VEC_TYPE==VCP_VEC_AVX
		inline static vcp_double_vec GetForceMask (const vcp_double_vec& m_r2, const vcp_double_vec& rc2, vcp_double_vec& j_mask)
		{
			static vcp_double_vec ones = vcp_simd_ones();
			vcp_double_vec result = vcp_simd_applymask( vcp_simd_and(vcp_simd_lt(m_r2, rc2), vcp_simd_neq(m_r2, vcp_simd_zerov()) ), j_mask);
			j_mask = ones;
			return result;
		}

		inline static vcp_double_vec InitJ_Mask (const size_t i)//calculations only for i+1 onwards.
		{
			switch (i & static_cast<size_t>(VCP_VEC_SIZE_M1)) {
				case 0: return _mm256_castsi256_pd(_mm256_set_epi32(~0, ~0, ~0, ~0, ~0, ~0, 0, 0));
				case 1: return _mm256_castsi256_pd(_mm256_set_epi32(~0, ~0, ~0, ~0, 0, 0, 0, 0));
				case 2: return _mm256_castsi256_pd(_mm256_set_epi32(~0, ~0, 0, 0, 0, 0, 0, 0));
				default: return vcp_simd_ones();
			}
		}
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
		// Erstellen der Bitmaske, analog zu Condition oben
		inline static vcp_double_vec GetForceMask(vcp_double_vec m_r2, vcp_double_vec rc2)
		{
			return vcp_simd_and(vcp_simd_lt(m_r2, rc2), vcp_simd_neq(m_r2, vcp_simd_zerov()) );
		}
#elif VCP_VEC_TYPE==VCP_NOVEC
		inline static vcp_double_vec GetForceMask(vcp_double_vec m_r2, vcp_double_vec rc2)
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

		inline static size_t InitJ (const size_t i)
		{
			return 0;
		}

		inline static bool DetectSingleCell ()
		{
			return false;
		}

		inline static vcp_double_vec GetForceMask (const vcp_double_vec& m_r2, const vcp_double_vec& rc2
#if VCP_VEC_TYPE==VCP_VEC_AVX
				, vcp_double_vec& sj_mask
#endif
				)
		{
			// Provide a mask with the same logic as used in
			// bool Condition(double m_r2, double rc2)
			return vcp_simd_lt(m_r2, rc2);
		}


#if VCP_VEC_TYPE==VCP_VEC_AVX
		inline static vcp_double_vec InitJ_Mask (const size_t i)
		{
			return vcp_simd_zerov(); //TODO: initj check avx <-> sse3
		}
#endif
 /* definition of GetForceMask */
	}; /* end of class CellPairPolicy_ */

}; /* end of class VectorizedCellProcessor */

#endif /* VECTORIZEDCELLPROCESSOR_H_ */
