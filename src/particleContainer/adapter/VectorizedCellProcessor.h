/**
 * \file
 * \brief VectorizedCellProcessor.h
 * \author Johannes Heckl, Wolfgang Eckhardt, Uwe Ehmann
 */

#ifndef VECTORIZEDCELLPROCESSOR_H_
#define VECTORIZEDCELLPROCESSOR_H_


#include "CellProcessor.h"
#include "utils/AlignedArray.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "vectorization/SIMD_TYPES.h"


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
	const double _cutoffRadiusSquare;

	/**
	 * \brief The squared LJ cutoff radius.
	 */
	const double _LJcutoffRadiusSquare;

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

	// lookup array for the distance molecule-molecule on a molecule-center basis.
	DoubleArray _centers_dist_lookup;

	// vector holding pointers to the above objects, used as a stack for
	// managing free objects
	std::vector<CellDataSoA*> _particleCellDataVector;

	/**
	 * \brief The body of the inner loop of the non-vectorized force calculation between LJ centers.
	 */
	template<class MacroPolicy>
		void _loopBodyNovecLJ(const CellDataSoA & soa1, size_t i, const CellDataSoA & soa2, size_t j, const double *const forceMask);

	/**
	 * \brief The body of the inner loop of the non-vectorized force calculation between charges.
	 * \author Robert Hajda
	 */
	template<class MacroPolicy>
		void _loopBodyNovecCharges(const CellDataSoA & soa1, size_t i, const CellDataSoA & soa2, size_t j, const double *const forceMask);

	/**
	 * \brief The body of the inner loop of the non-vectorized force calculation between charges and dipoles.
	 * \author Robert Hajda
	 */
	template<class MacroPolicy>
		void _loopBodyNovecChargesDipoles(const CellDataSoA & soa1, size_t i, const CellDataSoA & soa2, size_t j, const double *const forceMask, const bool& switched);

	/**
	 * \brief The body of the inner loop of the non-vectorized force calculation between dipoles.
	 * \author Robert Hajda
	 */
	template<class MacroPolicy>
		void _loopBodyNovecDipoles(const CellDataSoA & soa1, size_t i, const CellDataSoA & soa2, size_t j, const double *const forceMask);

	/**
	 * \brief Inner loop body of the non-vectorized force calculation between charges and quadrupoles.
	 * \author Uwe Ehmann
	 */
	template<class MacroPolicy>
		void _loopBodyNovecChargesQuadrupoles (const CellDataSoA& soa1, size_t i, const CellDataSoA& soa2, size_t j, const double *const forceMask, const bool& switched);

	/**
	 * \brief Inner loop body of the non-vectorized force calculation between dipoles and quadrupoles.
	 * \author Uwe Ehmann
	 */
	template<class MacroPolicy>
		void _loopBodyNovecDipolesQuadrupoles (const CellDataSoA& soa1, size_t i, const CellDataSoA& soa2, size_t j, const double *const forceMask, const bool& switched);

	/**
	 * \brief Inner loop body of the non-vectorized force calculation between quadrupoles.
	 * \author Uwe Ehmann
	 */
	template<class MacroPolicy>
		void _loopBodyNovecQuadrupoles (const CellDataSoA& soa1, size_t i, const CellDataSoA& soa2, size_t j, const double *const forceMask);
//TODO: try to merge this somehow
#if VCP_VEC_TYPE==VCP_VEC_SSE3
	template<class MacroPolicy>
	inline
	void _loopBodyLJ(
			const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
			const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
			const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
			const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
			vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
			vcp_double_vec& sum_upot6lj, vcp_double_vec& sum_virial,
			const vcp_double_vec& forceMask,
			const vcp_double_vec& e1s1, const vcp_double_vec& e2s2,
			const size_t& id_j0, const size_t& id_j1, const size_t& id_i);
#elif VCP_VEC_TYPE==VCP_VEC_AVX
	template<class MacroPolicy>
	inline
	void _loopBodyLJ(
			const vcp_double_vec& m1_r_x, const vcp_double_vec& m1_r_y, const vcp_double_vec& m1_r_z,
			const vcp_double_vec& r1_x, const vcp_double_vec& r1_y, const vcp_double_vec& r1_z,
			const vcp_double_vec& m2_r_x, const vcp_double_vec& m2_r_y, const vcp_double_vec& m2_r_z,
			const vcp_double_vec& r2_x, const vcp_double_vec& r2_y, const vcp_double_vec& r2_z,
			vcp_double_vec& f_x, vcp_double_vec& f_y, vcp_double_vec& f_z,
			vcp_double_vec& sum_upot6lj, vcp_double_vec& sum_virial,
			const vcp_double_vec& forceMask, const vcp_double_vec& e0s0,
			const vcp_double_vec& e1s1, const vcp_double_vec& e2s2, const vcp_double_vec& e3s3,
			const size_t& id_j0, const size_t& id_j1,  const size_t& id_j2,  const size_t& id_j3,
			const size_t& id_i);
#endif /* _loopBodyLJ SSE3 */

	/**
	 * \brief The dist lookup for a molecule and all centers of a type
	 * \author Robert Hajda
	 */
	template<class ForcePolicy>
#if VCP_VEC_TYPE==VCP_NOVEC
	unsigned long
#else
	vcp_double_vec
#endif
	calcDistLookup (const CellDataSoA & soa1, const size_t & i, const size_t & i_center_idx, const size_t & soa2_num_centers, const double & cutoffRadiusSquare,
			double* const soa2_center_dist_lookup, const double* const soa2_m_r_x, const double* const soa2_m_r_y, const double* const soa2_m_r_z
	#if VCP_VEC_TYPE!=VCP_VEC_NOVEC
			, const vcp_double_vec & cutoffRadiusSquareD, size_t end_j, const vcp_double_vec m1_r_x, const vcp_double_vec m1_r_y, const vcp_double_vec m1_r_z
	#endif
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
	 * static vcp_double_vec GetMask(vcp_double_vec m_r2, vcp_double_vec rc2);<br>
	 * Returns the mask indicating which pairs to calculate in the vectorized code.<br>
	 * <br>
	 * The MacroPolicy class must provide the following methods:<br>
	 * static bool MacroscopicValueCondition(double m_dx, double m_dy, double m_dz);<br>
	 * Returns whether to store macroscopic values for a non-vectorized pair.<br>
	 * <br>
	 * If the code is to be vectorized:<br>
	 * static vcp_double_vec GetMacroMask(vcp_double_vec forceMask, vcp_double_vec m_dx, vcp_double_vec m_dy, vcp_double_vec m_dz);
	 * <br> Returns the mask indicating for which pairs to store macroscopic values in<br>
	 * the vectorized code.
	 *
	 */
	template<class ForcePolicy, class MacroPolicy>
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

#if VCP_VEC_TYPE==VCP_VEC_AVX
		inline static vcp_double_vec GetForceMask (const vcp_double_vec& m_r2, const vcp_double_vec& rc2, vcp_double_vec& j_mask)
		{
			static vcp_double_vec ones = vcp_simd_ones();
			vcp_double_vec result = vcp_simd_and(
										vcp_simd_and(
											vcp_simd_lt(m_r2, rc2),
											vcp_simd_neq(m_r2, vcp_simd_zerov())),
									j_mask);
			j_mask = ones;
			return result;
		}
		inline static size_t InitJ (const size_t i)
		{
			return (i+1) & ~static_cast<size_t>(3); //TODO: initj check avx <-> sse3
		}
		inline static vcp_double_vec InitJ_Mask (const size_t i)
		{
			switch (i & static_cast<size_t>(3)) {
				case 0: return _mm256_castsi256_pd(_mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0, 0));
				case 1: return _mm256_castsi256_pd(_mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,0, 0, 0, 0));
				case 2: return _mm256_castsi256_pd(_mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0, 0, 0, 0, 0, 0));
				default: return vcp_simd_ones();
			}
		}
#elif VCP_VEC_TYPE==VCP_VEC_SSE3
		inline static size_t InitJ (const size_t i)
		{
			return i + (i & static_cast<size_t>(1)); //TODO: initj check avx <-> sse3
		}

		// Erstellen der Bitmaske, analog zu Condition oben
		inline static vcp_double_vec GetForceMask(vcp_double_vec m_r2, vcp_double_vec rc2)
		{
			return vcp_simd_and(vcp_simd_lt(m_r2, rc2), vcp_simd_neq(m_r2, vcp_simd_zerov()));
		}
#else //novec
		inline static size_t InitJ (const size_t i)
		{
			return i + 1;
		}
#endif /* definition of InitJ & GetForceMask */
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

#if VCP_VEC_TYPE!=VCP_NOVEC
		inline static vcp_double_vec GetForceMask (const vcp_double_vec& m_r2, const vcp_double_vec& rc2
			#if VCP_VEC_TYPE!=VCP_VEC_SSE3
				, vcp_double_vec& sj_mask
			#endif
				)
		{
			// Provide a mask with the same logic as used in
			// bool Condition(double m_r2, double rc2)
			return vcp_simd_lt(m_r2, rc2);//TODO: merge with sse3s
		}
	#if VCP_VEC_TYPE==VCP_VEC_AVX
		inline static vcp_double_vec InitJ_Mask (const size_t i)
		{
			return vcp_simd_zerov(); //TODO: initj check avx <-> sse3
		}
	#endif
#endif /* definition of GetForceMask */
	}; /* end of class CellPairPolicy_ */

	/**
	 * \brief A MacroPolicy for adding up all macroscopic values.
	 * \details This is used for single cell calculations and for<br>
	 * double cell calculations with no halo cells.
	 */
	class AllMacroPolicy_ {
	public:
		inline static bool MacroscopicValueCondition(double, double, double)
		{
			// We want all macroscopic values to be calculated.
			return true;
		}

		inline static bool MacroscopicValueConditionSwitched(double, double, double, bool)
		{
			return true;
		}

#if VCP_VEC_TYPE!=VCP_NOVEC
		inline static vcp_double_vec GetMacroMask(vcp_double_vec forceMask, vcp_double_vec, vcp_double_vec, vcp_double_vec)
		{
			// We want all macroscopic values to be calculated, but not those
			// for pairs which we ignore because of cutoff or other reasons.
			return forceMask;
		}

		inline static vcp_double_vec GetMacroMaskSwitched(vcp_double_vec forceMask, vcp_double_vec, vcp_double_vec, vcp_double_vec, vcp_double_vec)
		{
			return forceMask;
		}
#endif /* definition of GetMacroMask and GetMacroMaskSwitched */
	}; /* end of class AllMacroPolicy_ */

	/**
	 * \brief A MacroPolicy.
	 * \details Only adds up the macroscopic values for pairs where<br>
	 * mol1 < mol2 in the 3D ordering. This is used for double cell calculations<br>
	 * with one halo Cell.
	 */
	class SomeMacroPolicy_ {
	public:
		inline static bool MacroscopicValueCondition(double m_dx, double m_dy, double m_dz)
		{
			// Only calculate macroscopic values for pairs where molecule 1
			// "IsLessThan" molecule 2.
			return (m_dz < 0.0) || (
					(m_dz == 0.0) && (
						(m_dy < 0.0) || (
							(m_dy == 0.0) && (m_dx < 0.0)
						)
					)
				);
		}

		inline static bool MacroscopicValueConditionSwitched(double m_dx, double m_dy, double m_dz, bool switched)
		{
			// Only calculate macroscopic values for pairs where molecule 1
			// "IsLessThan" molecule 2.
			return MacroscopicValueCondition(m_dx, m_dy, m_dz) ^ switched;
		}

#if VCP_VEC_TYPE!=VCP_NOVEC
		// Only calculate macroscopic values for pairs where molecule 1
		// "IsLessThan" molecule 2.
		inline static vcp_double_vec GetMacroMask(const vcp_double_vec& forceMask, const vcp_double_vec& m_dx, const vcp_double_vec& m_dy, const vcp_double_vec& m_dz)
		{
			const vcp_double_vec zero = vcp_simd_zerov();

			const vcp_double_vec x_lt = vcp_simd_lt(m_dx, zero);
			const vcp_double_vec y_eq = vcp_simd_eq(m_dy, zero);
			const vcp_double_vec t1 = vcp_simd_and(x_lt, y_eq);

			const vcp_double_vec y_lt = vcp_simd_lt(m_dy, zero);
			const vcp_double_vec t2 = vcp_simd_or(t1, y_lt);

			const vcp_double_vec z_eq = vcp_simd_eq(m_dz, zero);
			const vcp_double_vec t3 = vcp_simd_and(t2, z_eq);

			const vcp_double_vec z_lt = vcp_simd_lt(m_dz, zero);
			const vcp_double_vec t4 = vcp_simd_or(t3, z_lt);

			return vcp_simd_and(t4, forceMask);
		}

		inline static vcp_double_vec GetMacroMaskSwitched(const vcp_double_vec& forceMask, const vcp_double_vec& m_dx, const vcp_double_vec& m_dy, const vcp_double_vec& m_dz, const vcp_double_vec& switched)
		{
			return vcp_simd_xor(GetMacroMask(forceMask, m_dx, m_dy, m_dz), switched);
		}
#endif /* definition of GetMacroMask and GetMacroMaskSwitched */
	}; /* end of class SomeMacroPolicy_ */
}; /* end of class VectorizedCellProcessor */

#endif /* VECTORIZEDCELLPROCESSOR_H_ */
