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
#include "WrapOpenMP.h"

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

	double processSingleMolecule(Molecule* /*m1*/, ParticleCell& /*cell2*/) {
		return 0.0;
	}

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

	struct VLJCPThreadData {
	public:
		VLJCPThreadData(): _ljc_dist_lookup(nullptr), _charges_dist_lookup(nullptr), _dipoles_dist_lookup(nullptr), _quadrupoles_dist_lookup(nullptr){
			_upot6ljV.resize(_numVectorElements);
			_upotXpolesV.resize(_numVectorElements);
			_virialV.resize(_numVectorElements);
			_myRFV.resize(_numVectorElements);

			for (size_t j = 0; j < _numVectorElements; ++j) {
				_upot6ljV[j] = 0.0;
				_upotXpolesV[j] = 0.0;
				_virialV[j] = 0.0;
				_myRFV[j] = 0.0;
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

		AlignedArray<double> _upot6ljV, _upotXpolesV, _virialV, _myRFV;
	};

	std::vector<VLJCPThreadData *> _threadData;

	static const size_t _numVectorElements = VCP_VEC_SIZE;
	size_t _numThreads;

	template<bool calculateMacroscopic>
	inline void _loopBodyLJ(
			const DoubleVec& m1_r_x, const DoubleVec& m1_r_y, const DoubleVec& m1_r_z,
			const DoubleVec& r1_x, const DoubleVec& r1_y, const DoubleVec& r1_z,
			const DoubleVec& m2_r_x, const DoubleVec& m2_r_y, const DoubleVec& m2_r_z,
			const DoubleVec& r2_x, const DoubleVec& r2_y, const DoubleVec& r2_z,
			DoubleVec& f_x, DoubleVec& f_y, DoubleVec& f_z,
			DoubleVec& V_x, DoubleVec& V_y, DoubleVec& V_z,
			DoubleVec& sum_upot6lj, DoubleVec& sum_virial,
			const MaskVec& forceMask,
			const DoubleVec& eps_24, const DoubleVec& sig2,
			const DoubleVec& shift6);

	template<bool calculateMacroscopic>
	inline void _loopBodyCharge(
		const DoubleVec& m1_r_x, const DoubleVec& m1_r_y, const DoubleVec& m1_r_z,
		const DoubleVec& r1_x, const DoubleVec& r1_y, const DoubleVec& r1_z,
		const DoubleVec& qii,
		const DoubleVec& m2_r_x, const DoubleVec& m2_r_y, const DoubleVec& m2_r_z,
		const DoubleVec& r2_x, const DoubleVec& r2_y, const DoubleVec& r2_z,
		const DoubleVec& qjj,
		DoubleVec& f_x, DoubleVec& f_y, DoubleVec& f_z,
		DoubleVec& V_x, DoubleVec& V_y, DoubleVec& V_z,
		DoubleVec& sum_upotXpoles, DoubleVec& sum_virial,
		const MaskVec& forceMask);

	template<bool calculateMacroscopic>
	inline void _loopBodyChargeDipole(
		const DoubleVec& m1_r_x, const DoubleVec& m1_r_y, const DoubleVec& m1_r_z,
		const DoubleVec& r1_x, const DoubleVec& r1_y, const DoubleVec& r1_z,
		const DoubleVec& q,
		const DoubleVec& m2_r_x, const DoubleVec& m2_r_y, const DoubleVec& m2_r_z,
		const DoubleVec& r2_x, const DoubleVec& r2_y, const DoubleVec& r2_z,
		const DoubleVec& e_x, const DoubleVec& e_y, const DoubleVec& e_z,
		const DoubleVec& p,
		DoubleVec& f_x, DoubleVec& f_y, DoubleVec& f_z,
		DoubleVec& V_x, DoubleVec& V_y, DoubleVec& V_z,
		DoubleVec& M_x, DoubleVec& M_y, DoubleVec& M_z,
		DoubleVec& sum_upotXpoles, DoubleVec& sum_virial,
		const MaskVec& forceMask);

	template<bool calculateMacroscopic>
	inline void _loopBodyDipole(
		const DoubleVec& m1_r_x, const DoubleVec& m1_r_y, const DoubleVec& m1_r_z,
		const DoubleVec& r1_x, const DoubleVec& r1_y, const DoubleVec& r1_z,
		const DoubleVec& eii_x, const DoubleVec& eii_y, const DoubleVec& eii_z,
		const DoubleVec& pii,
		const DoubleVec& m2_r_x, const DoubleVec& m2_r_y, const DoubleVec& m2_r_z,
		const DoubleVec& r2_x, const DoubleVec& r2_y, const DoubleVec& r2_z,
		const DoubleVec& ejj_x, const DoubleVec& ejj_y, const DoubleVec& ejj_z,
		const DoubleVec& pjj,
		DoubleVec& f_x, DoubleVec& f_y, DoubleVec& f_z,
		DoubleVec& V_x, DoubleVec& V_y, DoubleVec& V_z,
		DoubleVec& M1_x, DoubleVec& M1_y, DoubleVec& M1_z,
		DoubleVec& M2_x, DoubleVec& M2_y, DoubleVec& M2_z,
		DoubleVec& sum_upotXpoles, DoubleVec& sum_virial, DoubleVec& sum_myRF,
		const MaskVec& forceMask,
		const DoubleVec& epsRFInvrc3);

	template<bool calculateMacroscopic>
	inline void _loopBodyChargeQuadrupole(
		const DoubleVec& m1_r_x, const DoubleVec& m1_r_y, const DoubleVec& m1_r_z,
		const DoubleVec& r1_x, const DoubleVec& r1_y, const DoubleVec& r1_z,
		const DoubleVec& q,
		const DoubleVec& m2_r_x, const DoubleVec& m2_r_y, const DoubleVec& m2_r_z,
		const DoubleVec& r2_x, const DoubleVec& r2_y, const DoubleVec& r2_z,
		const DoubleVec& ejj_x, const DoubleVec& ejj_y, const DoubleVec& ejj_z,
		const DoubleVec& m,
		DoubleVec& f_x, DoubleVec& f_y, DoubleVec& f_z,
		DoubleVec& V_x, DoubleVec& V_y, DoubleVec& V_z,
		DoubleVec& M_x, DoubleVec& M_y, DoubleVec& M_z,
		DoubleVec& sum_upotXpoles, DoubleVec& sum_virial,
		const MaskVec& forceMask);

	template<bool calculateMacroscopic>
	inline void _loopBodyDipoleQuadrupole(
		const DoubleVec& m1_r_x, const DoubleVec& m1_r_y, const DoubleVec& m1_r_z,
		const DoubleVec& r1_x, const DoubleVec& r1_y, const DoubleVec& r1_z,
		const DoubleVec& eii_x, const DoubleVec& eii_y, const DoubleVec& eii_z,
		const DoubleVec& p,
		const DoubleVec& m2_r_x, const DoubleVec& m2_r_y, const DoubleVec& m2_r_z,
		const DoubleVec& r2_x, const DoubleVec& r2_y, const DoubleVec& r2_z,
		const DoubleVec& ejj_x, const DoubleVec& ejj_y, const DoubleVec& ejj_z,
		const DoubleVec& m,
		DoubleVec& f_x, DoubleVec& f_y, DoubleVec& f_z,
		DoubleVec& V_x, DoubleVec& V_y, DoubleVec& V_z,
		DoubleVec& M1_x, DoubleVec& M1_y, DoubleVec& M1_z,
		DoubleVec& M2_x, DoubleVec& M2_y, DoubleVec& M2_z,
		DoubleVec& sum_upotXpoles, DoubleVec& sum_virial,
		const MaskVec& forceMask);

	template<bool calculateMacroscopic>
	inline void _loopBodyQuadrupole(
		const DoubleVec& m1_r_x, const DoubleVec& m1_r_y, const DoubleVec& m1_r_z,
		const DoubleVec& r1_x, const DoubleVec& r1_y, const DoubleVec& r1_z,
		const DoubleVec& eii_x, const DoubleVec& eii_y, const DoubleVec& eii_z,
		const DoubleVec& mii,
		const DoubleVec& m2_r_x, const DoubleVec& m2_r_y, const DoubleVec& m2_r_z,
		const DoubleVec& r2_x, const DoubleVec& r2_y, const DoubleVec& r2_z,
		const DoubleVec& ejj_x, const DoubleVec& ejj_y, const DoubleVec& ejj_z,
		const DoubleVec& mjj,
		DoubleVec& f_x, DoubleVec& f_y, DoubleVec& f_z,
		DoubleVec& V_x, DoubleVec& V_y, DoubleVec& V_z,
		DoubleVec& Mii_x, DoubleVec& Mii_y, DoubleVec& Mii_z,
		DoubleVec& Mjj_x, DoubleVec& Mjj_y, DoubleVec& Mjj_z,
		DoubleVec& sum_upotXpoles, DoubleVec& sum_virial,
		const MaskVec& forceMask);

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
	 * static DoubleVec GetForceMask(DoubleVec m_r2, DoubleVec rc2);<br>
	 * Returns the mask indicating which pairs to calculate in the vectorized code.<br>
	 * <br>
	 * The boolean CalculateMacroscopic should specify, whether macroscopic values are to be calculated or not.
	 * <br>
	 * The class MaskGatherChooser is a class, that specifies the used loading,storing and masking routines.
	 */
	template<class ForcePolicy, bool CalculateMacroscopic, class MaskGatherChooser>
	void _calculatePairs(const CellDataSoA & soa1, const CellDataSoA & soa2);

}; /* end of class VectorizedCellProcessor */

#endif /* VECTORIZEDCELLPROCESSOR_H_ */
