/**
 * \file
 * \brief VCPADR.h
 * \author Johannes Heckl, Wolfgang Eckhardt, Uwe Ehmann, Steffen Seckler, Alex Hocks
 */

#ifndef VCPADR_H_
#define VCPADR_H_

#include <iostream>
#include <vector>
#include <array>
#include <cmath>

#include "particleContainer/adapter/CellProcessor.h"
#include "particleContainer/adapter/vectorization/SIMD_VectorizedCellProcessorHelpers.h"

class Component;
class Domain;
class Comp2Param;
class CellDataSoA;
namespace Resolution {
	class Handler;
}

/**
 * \brief Vectorized calculation of the force.
 * \author Johannes Heckl
 */
class VCPADR : public CellProcessor {
private:
	enum SiteType {
		LJ = 0, CHARGE = 1, DIPOLE = 2, QUADRUPOLE = 3, NUM_TYPES
	};
public:
	typedef std::vector<Component> ComponentList;

	VCPADR& operator=(const VCPADR&) = delete;

	/**
	 * \brief Construct and set up the internal parameter table.
	 * \details Components and parameters should be finalized before this call.
	 */
	VCPADR(Domain & domain, double cutoffRadius, double LJcutoffRadius, const Resolution::Handler& resolutionHandler);

	~VCPADR();

	/**
	 * Initializes parameter tables. Must be called after Comp2Param was initialized, but before first force computation
	 * during simulation initialization.
	 * */
	void init();

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

        void processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll = false);

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
	 * \brief Masks to only allow same resolution interaction within hybrid
	 * */
	std::array<AlignedArray<vcp_mask_single>, SiteType::NUM_TYPES> _compMask;
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
		 * To utilize the gather operations of the KNL architecture, the dist_lookup is able to store the indices of the required particles.
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

		AlignedArray<vcp_real_accum> _upot6ljV, _upotXpolesV, _virialV, _myRFV;
	};

	std::vector<VLJCPThreadData *> _threadData;

	static const size_t _numVectorElements = VCP_VEC_SIZE;
	size_t _numThreads;

	//! @brief AdResS resolution handler that manages the different regions and their corresponding resolutions
	const Resolution::Handler& _resolutionHandler;

	template<bool calculateMacroscopic>
	inline void _loopBodyLJ(
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

	template<bool calculateMacroscopic>
	inline void _loopBodyCharge(
		const RealCalcVec& m1_r_x, const RealCalcVec& m1_r_y, const RealCalcVec& m1_r_z,
		const RealCalcVec& r1_x, const RealCalcVec& r1_y, const RealCalcVec& r1_z,
		const RealCalcVec& qii,
		const RealCalcVec& m2_r_x, const RealCalcVec& m2_r_y, const RealCalcVec& m2_r_z,
		const RealCalcVec& r2_x, const RealCalcVec& r2_y, const RealCalcVec& r2_z,
		const RealCalcVec& qjj,
		RealCalcVec& f_x, RealCalcVec& f_y, RealCalcVec& f_z,
		RealAccumVec& V_x, RealAccumVec& V_y, RealAccumVec& V_z,
		RealAccumVec& sum_upotXpoles, RealAccumVec& sum_virial,
		const MaskCalcVec& forceMask);

	template<bool calculateMacroscopic>
	inline void _loopBodyChargeDipole(
		const RealCalcVec& m1_r_x, const RealCalcVec& m1_r_y, const RealCalcVec& m1_r_z,
		const RealCalcVec& r1_x, const RealCalcVec& r1_y, const RealCalcVec& r1_z,
		const RealCalcVec& q,
		const RealCalcVec& m2_r_x, const RealCalcVec& m2_r_y, const RealCalcVec& m2_r_z,
		const RealCalcVec& r2_x, const RealCalcVec& r2_y, const RealCalcVec& r2_z,
		const RealCalcVec& e_x, const RealCalcVec& e_y, const RealCalcVec& e_z,
		const RealCalcVec& p,
		RealCalcVec& f_x, RealCalcVec& f_y, RealCalcVec& f_z,
		RealAccumVec& V_x, RealAccumVec& V_y, RealAccumVec& V_z,
		RealAccumVec& M_x, RealAccumVec& M_y, RealAccumVec& M_z,
		RealAccumVec& sum_upotXpoles, RealAccumVec& sum_virial,
		const MaskCalcVec& forceMask);

	template<bool calculateMacroscopic>
	inline void _loopBodyDipole(
		const RealCalcVec& m1_r_x, const RealCalcVec& m1_r_y, const RealCalcVec& m1_r_z,
		const RealCalcVec& r1_x, const RealCalcVec& r1_y, const RealCalcVec& r1_z,
		const RealCalcVec& eii_x, const RealCalcVec& eii_y, const RealCalcVec& eii_z,
		const RealCalcVec& pii,
		const RealCalcVec& m2_r_x, const RealCalcVec& m2_r_y, const RealCalcVec& m2_r_z,
		const RealCalcVec& r2_x, const RealCalcVec& r2_y, const RealCalcVec& r2_z,
		const RealCalcVec& ejj_x, const RealCalcVec& ejj_y, const RealCalcVec& ejj_z,
		const RealCalcVec& pjj,
		RealCalcVec& f_x, RealCalcVec& f_y, RealCalcVec& f_z,
		RealAccumVec& V_x, RealAccumVec& V_y, RealAccumVec& V_z,
		RealAccumVec& M1_x, RealAccumVec& M1_y, RealAccumVec& M1_z,
		RealAccumVec& M2_x, RealAccumVec& M2_y, RealAccumVec& M2_z,
		RealAccumVec& sum_upotXpoles, RealAccumVec& sum_virial, RealAccumVec& sum_myRF,
		const MaskCalcVec& forceMask,
		const RealCalcVec& epsRFInvrc3);

	template<bool calculateMacroscopic>
	inline void _loopBodyChargeQuadrupole(
		const RealCalcVec& m1_r_x, const RealCalcVec& m1_r_y, const RealCalcVec& m1_r_z,
		const RealCalcVec& r1_x, const RealCalcVec& r1_y, const RealCalcVec& r1_z,
		const RealCalcVec& q,
		const RealCalcVec& m2_r_x, const RealCalcVec& m2_r_y, const RealCalcVec& m2_r_z,
		const RealCalcVec& r2_x, const RealCalcVec& r2_y, const RealCalcVec& r2_z,
		const RealCalcVec& ejj_x, const RealCalcVec& ejj_y, const RealCalcVec& ejj_z,
		const RealCalcVec& m,
		RealCalcVec& f_x, RealCalcVec& f_y, RealCalcVec& f_z,
		RealAccumVec& V_x, RealAccumVec& V_y, RealAccumVec& V_z,
		RealAccumVec& M_x, RealAccumVec& M_y, RealAccumVec& M_z,
		RealAccumVec& sum_upotXpoles, RealAccumVec& sum_virial,
		const MaskCalcVec& forceMask);

	template<bool calculateMacroscopic>
	inline void _loopBodyDipoleQuadrupole(
		const RealCalcVec& m1_r_x, const RealCalcVec& m1_r_y, const RealCalcVec& m1_r_z,
		const RealCalcVec& r1_x, const RealCalcVec& r1_y, const RealCalcVec& r1_z,
		const RealCalcVec& eii_x, const RealCalcVec& eii_y, const RealCalcVec& eii_z,
		const RealCalcVec& p,
		const RealCalcVec& m2_r_x, const RealCalcVec& m2_r_y, const RealCalcVec& m2_r_z,
		const RealCalcVec& r2_x, const RealCalcVec& r2_y, const RealCalcVec& r2_z,
		const RealCalcVec& ejj_x, const RealCalcVec& ejj_y, const RealCalcVec& ejj_z,
		const RealCalcVec& m,
		RealCalcVec& f_x, RealCalcVec& f_y, RealCalcVec& f_z,
		RealAccumVec& V_x, RealAccumVec& V_y, RealAccumVec& V_z,
		RealAccumVec& M1_x, RealAccumVec& M1_y, RealAccumVec& M1_z,
		RealAccumVec& M2_x, RealAccumVec& M2_y, RealAccumVec& M2_z,
		RealAccumVec& sum_upotXpoles, RealAccumVec& sum_virial,
		const MaskCalcVec& forceMask);

	template<bool calculateMacroscopic>
	inline void _loopBodyQuadrupole(
		const RealCalcVec& m1_r_x, const RealCalcVec& m1_r_y, const RealCalcVec& m1_r_z,
		const RealCalcVec& r1_x, const RealCalcVec& r1_y, const RealCalcVec& r1_z,
		const RealCalcVec& eii_x, const RealCalcVec& eii_y, const RealCalcVec& eii_z,
		const RealCalcVec& mii,
		const RealCalcVec& m2_r_x, const RealCalcVec& m2_r_y, const RealCalcVec& m2_r_z,
		const RealCalcVec& r2_x, const RealCalcVec& r2_y, const RealCalcVec& r2_z,
		const RealCalcVec& ejj_x, const RealCalcVec& ejj_y, const RealCalcVec& ejj_z,
		const RealCalcVec& mjj,
		RealCalcVec& f_x, RealCalcVec& f_y, RealCalcVec& f_z,
		RealAccumVec& V_x, RealAccumVec& V_y, RealAccumVec& V_z,
		RealAccumVec& Mii_x, RealAccumVec& Mii_y, RealAccumVec& Mii_z,
		RealAccumVec& Mjj_x, RealAccumVec& Mjj_y, RealAccumVec& Mjj_z,
		RealAccumVec& sum_upotXpoles, RealAccumVec& sum_virial,
		const MaskCalcVec& forceMask);

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
	void _calculatePairs(CellDataSoA & soa1, CellDataSoA & soa2);

	/**
 	 * unpacks mask i and mask j from the compMap array according to the index array id_i and id_j (for mic+avx2: use gather)
 	 * @param mask_i mask i dst
 	 * @param mask_j mask j dst
 	 * @param compMap initial comp array
 	 * @param id_i array of displacements for i
 	 * @param offset_i offset of the id_i array
 	 * @param id_j array of displacements for j
 	 * @param offset_j offset of the id_j array
 	 */
	static vcp_inline
	void unpackComp(MaskCalcVec& mask_i, MaskCalcVec& mask_j, AlignedArray<vcp_mask_single>& compMap_i, AlignedArray<vcp_mask_single>& compMap_j,
					const vcp_center_id_t* const id_i, const vcp_center_id_t& offset_i,
					const vcp_center_id_t* const id_j, const vcp_center_id_t& offset_j,
					const vcp_lookupOrMask_vec& lookupORforceMask __attribute__((unused)));

}; /* end of class VectorizedCellProcessor */

#endif /* VECTORIZEDCELLPROCESSOR_H_ */
