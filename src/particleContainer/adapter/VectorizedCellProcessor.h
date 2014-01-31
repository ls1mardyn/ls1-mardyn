/**
 * \file
 * \brief VectorizedCellProcessor.h
 * \author Johannes Heckl, Wolfgang Eckhardt
 */

#ifndef VECTORIZEDCELLPROCESSOR_H_
#define VECTORIZEDCELLPROCESSOR_H_

#include "CellProcessor.h"
#include "utils/AlignedArray.h"
#include <iostream>
#include <vector>
#include <cmath>

// The following error should NEVER occur, since it signalizes, that the macros, used by THIS translation unit are defined anywhere else in the program.
#if defined(VLJCP_VEC_TYPE) || defined(VLJCP_NOVEC) || defined(VLJCP_VEC_SSE3)
	#error conflicting macro definitions
#endif

#define VLJCP_NOVEC 0
#define VLJCP_VEC_SSE3 1
#define VLJCP_VEC_AVX 2

#if defined(__AVX__) && not defined(AVX128)
	#define VLJCP_VEC_TYPE VLJCP_VEC_AVX
#elif defined(__AVX__) && defined(AVX128)
	#define VLJCP_VEC_TYPE VLJCP_VEC_SSE3
#elif defined(__SSE3__)
	#define VLJCP_VEC_TYPE VLJCP_VEC_SSE3
#else
	#define VLJCP_VEC_TYPE VLJCP_NOVEC
#endif

#ifdef NOVEC
	#ifdef VLJCP_VEC_TYPE
		#warn Multiple vectorization methods specified. Will not use vectorization at all!
		#undef VLJCP_VEC_TYPE
	#endif
	#define VLJCP_VEC_TYPE VLJCP_NOVEC
#endif

// Include necessary files if we vectorize.
#if VLJCP_VEC_TYPE==VLJCP_VEC_AVX
	#include "immintrin.h"
#elif VLJCP_VEC_TYPE==VLJCP_VEC_SSE3
	#include "pmmintrin.h"
#endif


class Component;
class Domain;
class Comp2Param;
class Molecule;
class CellDataSoA;

/**
 * \brief Vectorized calculation of Lennard Jones force.
 * \author Johannes Heckl
 */
class VectorizedCellProcessor : public CellProcessor {
public:
	typedef std::vector<Component> ComponentList;

	/**
	 * \brief Construct and set up the internal parameter table.
	 * \details Components and parameters should be finalized before this call.
	 */
	VectorizedCellProcessor(Domain & domain, double cutoffRadius);

	~VectorizedCellProcessor();

	/**
	 * \brief Reset macroscopic values to 0.0.
	 */
	void initTraversal(const size_t numCells);
	/**
	 * \brief Load the LennardJonesSoA for cell.
	 */
	void preprocessCell(ParticleCell& cell);
	/**
	 * \brief Calculate forces between pairs of Molecules in cell1 and cell2.
	 */
	void processCellPair(ParticleCell& cell1, ParticleCell& cell2);
	/**
	 * \brief Calculate forces between pairs of Molecules in cell.
	 */
	void processCell(ParticleCell& cell);

	double processSingleMolecule(Molecule* m1, ParticleCell& cell2);

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
	const double _rc2;
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
	 * \brief Sum of all potentials.
	 * \details Multiplied by 6.0 for performance reasons.
	 */
	double _upot6lj;
	/**
	 * \brief The virial.
	 */
	double _virial;

	// vector holding pointers to the above objects, used as a stack for
	// managing free objects
	std::vector<CellDataSoA*> _particleCellDataVector;

	// lookup array for the distance molecule-molecule on a molecule-center basis.
	DoubleArray _center_dist_lookup;

	/**
	 * \brief The body of the inner loop of the non-vectorized force calculation.
	 */
	template<class ForcePolicy, class MacroPolicy>
		void _loopBodyNovec(const CellDataSoA & soa1, size_t i, const CellDataSoA & soa2, size_t j, const double *const forceMask);

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
	 * static __m128d GetMask(__m128d m_r2, __m128d rc2);<br>
	 * Returns the mask indicating which pairs to calculate in the vectorized code.<br>
	 * <br>
	 * The MacroPolicy class must provide the following methods:<br>
	 * static bool MacroscopicValueCondition(double m_dx, double m_dy, double m_dz);<br>
	 * Returns whether to store macroscopic values for a non-vectorized pair.<br>
	 * <br>
	 * If the code is to be vectorized:<br>
	 * static __m128d GetMacroMask(__m128d forceMask, __m128d m_dx, __m128d m_dy, __m128d m_dz);
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

#if VLJCP_VEC_TYPE==VLJCP_VEC_AVX
		inline static __m256d GetForceMask (const __m256d& m_r2, const __m256d& rc2, __m256d& j_mask)
		{
			static __m256d ones = _mm256_castsi256_pd( _mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF) );
			__m256d result = _mm256_and_pd(
									_mm256_and_pd(
											_mm256_cmp_pd(m_r2, rc2, _CMP_LT_OS),
											_mm256_cmp_pd(m_r2, _mm256_setzero_pd(), _CMP_NEQ_OS)),
									j_mask);
			j_mask = ones;
			return result;
		}
		inline static size_t InitJ (const size_t i)
		{
			return (i+1) & ~static_cast<size_t>(3);
		}
		inline static __m256d InitJ_Mask (const size_t i)
		{
			switch (i & static_cast<size_t>(3)) {
				case 0: return _mm256_castsi256_pd(_mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0, 0));
				case 1: return _mm256_castsi256_pd(_mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF,0, 0, 0, 0));
				case 2: return _mm256_castsi256_pd(_mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0, 0, 0, 0, 0, 0));
				default: return _mm256_castsi256_pd(_mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)); 
			}
		}
#else
#if VLJCP_VEC_TYPE==VLJCP_VEC_SSE3
		inline static size_t InitJ (const size_t i)
		{
			return i + (i & static_cast<size_t>(1));
		}
		inline static __m128d GetForceMask(__m128d m_r2, __m128d rc2)
		{
			return _mm_and_pd(_mm_cmplt_pd(m_r2, rc2), _mm_cmpneq_pd(m_r2, _mm_setzero_pd()));
		}
#else
		inline static size_t InitJ (const size_t i)
		{
			return i + 1;
		}
#endif
#endif
	};


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

#if VLJCP_VEC_TYPE==VLJCP_VEC_AVX
		inline static __m256d GetForceMask (const __m256d& m_r2, const __m256d& rc2, __m256d& j_mask)
		{
			return _mm256_cmp_pd(m_r2, rc2, _CMP_LT_OS);
		}
		inline static __m256d InitJ_Mask (const size_t i)
		{
			return _mm256_setzero_pd();
		}
#else
#if VLJCP_VEC_TYPE==VLJCP_VEC_SSE3
		inline static __m128d GetForceMask(__m128d m_r2, __m128d rc2)
		{
			// Provide a mask with the same logic as used in
			// bool Condition(double m_r2, double rc2)
			return _mm_cmplt_pd(m_r2, rc2);
		}
#endif
#endif
	};
	/**
	 * \brief A MacroPolicy for adding up all macroscopic values.
	 * \details This is used for single cell calculations and for<br>
	 * double cell calculations with no halo cells.
	 */
	class AllMacroPolicy_ {
	public:
		static bool MacroscopicValueCondition(double, double, double)
		{
			// We want all macroscopic values to be calculated.
			return true;
		}
#if VLJCP_VEC_TYPE==VLJCP_VEC_SSE3 || VLJCP_VEC_TYPE==VLJCP_VEC_AVX
		static __m128d GetMacroMask(__m128d forceMask, __m128d, __m128d, __m128d)
		{
			// We want all macroscopic values to be calculated, but not those
			// for pairs which we ignore because of cutoff or other reasons.
			return forceMask;
		}
#endif
#if VLJCP_VEC_TYPE==VLJCP_VEC_AVX
		inline static __m256d GetMacroMask(const __m256d& forceMask, const __m256d&, const __m256d&, const __m256d&)
		{
			return forceMask;
		}
#endif
	};
	/**
	 * \brief A MacroPolicy.
	 * \details Only adds up the macroscopic values for pairs where<br>
	 * mol1 < mol2 in the 3D ordering. This is used for double cell calculations<br>
	 * with one halo Cell.
	 */
	class SomeMacroPolicy_ {
	public:
		static bool MacroscopicValueCondition(double m_dx, double m_dy, double m_dz)
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
#if VLJCP_VEC_TYPE==VLJCP_VEC_SSE3 || VLJCP_VEC_TYPE==VLJCP_VEC_AVX
		// Only calculate macroscopic values for pairs where molecule 1
		// "IsLessThan" molecule 2.
		static __m128d GetMacroMask(__m128d forceMask, __m128d m_dx, __m128d m_dy, __m128d m_dz)
		{
			const __m128d zero = _mm_setzero_pd();

			const __m128d x_lt = _mm_cmplt_pd(m_dx, zero);
			const __m128d y_eq = _mm_cmpeq_pd(m_dy, zero);
			const __m128d t1 = _mm_and_pd(y_eq, x_lt);

			const __m128d y_lt = _mm_cmplt_pd(m_dy, zero);
			const __m128d t2 = _mm_or_pd(y_lt, t1);

			const __m128d z_eq = _mm_cmpeq_pd(m_dz, zero);
			const __m128d t3 = _mm_and_pd(z_eq, t2);

			const __m128d z_lt = _mm_cmplt_pd(m_dz, zero);
			const __m128d t4 = _mm_or_pd(z_lt, t3);

			return _mm_and_pd(forceMask, t4);
		}
#endif
#if VLJCP_VEC_TYPE==VLJCP_VEC_AVX
		inline static __m256d GetMacroMask(const __m256d& forceMask, const __m256d& m_dx, const __m256d& m_dy, const __m256d& m_dz)
		{
			const __m256d zero = _mm256_setzero_pd();

			const __m256d x_lt = _mm256_cmp_pd(m_dx, zero, _CMP_LT_OS);
			const __m256d y_eq = _mm256_cmp_pd(m_dy, zero, _CMP_EQ_OS);
			const __m256d t1 = _mm256_and_pd(x_lt, y_eq);

			const __m256d y_lt = _mm256_cmp_pd(m_dy, zero, _CMP_LT_OS);
			const __m256d t2 = _mm256_or_pd(t1, y_lt);

			const __m256d z_eq = _mm256_cmp_pd(m_dz, zero, _CMP_EQ_OS);
			const __m256d t3 = _mm256_and_pd(t2, z_eq);

			const __m256d z_lt = _mm256_cmp_pd(m_dz, zero, _CMP_LT_OS);
			const __m256d t4 = _mm256_or_pd(t3, z_lt);

			return _mm256_and_pd(t4, forceMask);
		}
#endif
	};
};

#endif
