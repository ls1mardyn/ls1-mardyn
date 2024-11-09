//
// Created by alex on 09.05.23.
//

#ifndef MARDYN_ADRESSFORCEADAPTER_H
#define MARDYN_ADRESSFORCEADAPTER_H

#include "molecules/Molecule.h"
#include "molecules/Comp2Param.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"

#include "plugins/AdResS/util/Region.h"
#include "plugins/AdResS/features/Resolution.h"

#include "plugins/PotentialOfMeanForce/IBI_Math.h"

/**
 * Handles all force calculation done in AdResS method.
 * */
class AdResSForceAdapter : public ParticlePairsHandler {
public:
    /**
     * Constructor
     * @param resolutionHandler Resolution Handler reference from AdResS plugin
     * */
    explicit AdResSForceAdapter(Resolution::Handler& resolutionHandler, const std::string& cgForcePath, const std::string& cgPotPath);

    /**
     * Destructor: clears allocated memory
     * */
    ~AdResSForceAdapter() noexcept;

	/**
	 * @brief initializes the local Comp2Param instance by using the provided domain.
	 * _c2p is used to force calculation
	 * */
    void init() override;

    /**
     * @brief reduction of all thread local buffers, called when all forces are computed.
     * Will set results in _mesoValues.
     * */
    void finish();

    /**
     * public interface method to compute force between two molecules
     * will find in which FPRegion the molecules are and delegate to the corresponding AdResS processPair function
     * @param molecule1       molecule 1
     * @param molecule2       molecule 2
     * @param distanceVector  distance vector from molecule 2 to molecule 1
     * @param pairType        molecule pair type (see PairType)
     * @param dd              square of the distance between the two molecules
     * @param calculateLJ     true if we shall calculate the LJ interaction, otherwise false (default true)
     * */
    double processPair(Molecule& particle1, Molecule& particle2, double distanceVector[3], PairType pairType, double dd, bool calculateLJ) override;

    /** @brief based on ParticlePairs2PotForceAdapter.h -> Adapted to support AdResS
     * calculate force between pairs and collect macroscopic contribution
     *
     * For all pairs, the force between the two Molecules has to be calculated
     * and stored in the molecules. For original pairs(pairType 0), the contributions
     * to the macroscopic values have to be collected
     *
     * @param molecule1       molecule 1
     * @param molecule2       molecule 2
     * @param distanceVector  distance vector from molecule 2 to molecule 1
     * @param pairType        molecule pair type (see PairType)
     * @param dd              square of the distance between the two molecules
     * @param calculateLJ     true if we shall calculate the LJ interaction, otherwise false (default true)
     * @param compResMap      a map from AdResS plugin instance that maps all component ids to the resolution
     * @param noHybrid        flag if there exists a hybrid molecule, computed in interface processPair
     *
     * @return                interaction energy
     */
    double
    processPair(Molecule &molecule1, Molecule &molecule2, double * distanceVector, PairType pairType,
                double dd, bool calculateLJ, const std::vector<Resolution::ResolutionType> &compResMap, bool noHybrid,
				const Resolution::FPRegion &region);

    /** @brief Based on potforce.h::PotForce
     * Calculate potential and force between two molecules including all site-site interactions.
     *
     * Calculates the potential energy and force between two molecules i and j.
     * The interaction paramters are held in precomputed streams, which are initialized in Comp2Param::initialize
     *
     * @param[in]  mi   molecule i
     * @param[in]  mj   molecule j
     * @param[in]  params    reference to the corresponding interaction parameter stream initialized via Comp2Param::initialize
     * @param[in]  drm   distance vector from molecule j to molecule i
     * @param[out] Upot6LJ   potential energy resulting from Lennard Jones interactions
     * @param[out] UpotXpoles   potential energy resulting from Charge, Dipole and Quadrupole interactions
     * @param[out] MyRF
     * @todo Document parameter (is this reaction field?)
     * @param[out] Virial   Virial
     * @param[in]  caculateLJ    enable or disable calculation of Lennard Jones interactions
     * @param[in] noHybrid if true both molecules must be non hybrid
     * @param[in] compResMap      a map from AdResS plugin instance that maps all component ids to the resolution
     */
    void potForce(Molecule &mi, Molecule &mj, ParaStrm &params, ParaStrm &paramInv, double * drm, double &Upot6LJ,
                         double &UpotXpoles,
                         double &MyRF, double Virial[3], bool calculateLJ, bool noHybrid,
                         const std::vector<Resolution::ResolutionType> &compResMap, const Resolution::FPRegion &region);

    /** @brief Based on potforce.h::FluidPot
     * Calculate potential between two molecules including all site-site interactions.
     *
     * Calculates the potential energy and force between two molecules i and j.
     * The interaction paramters are held in precomputed streams, which are initialized in Comp2Param::initialize
     *
     * @param[in]  mi   molecule i
     * @param[in]  mj   molecule j
     * @param[in]  params    reference to the corresponding interaction parameter stream initialized via Comp2Param::initialize
     * @param[in]  drm   distance vector from molecule j to molecule i
     * @param[out] Upot6LJ   potential energy resulting from Lennard Jones interactions
     * @param[out] UpotXpoles   potential energy resulting from Charge, Dipole and Quadrupole interactions
     * @param[out] MyRF
     * @todo Document parameter (is this reaction field?)
     * @param[out] Virial   Virial
     * @param[in]  caculateLJ    enable or disable calculation of Lennard Jones interactions
     * @param[in] noHybrid if true both molecules must be non hybrid
     * @param[in] compResMap      a map from AdResS plugin instance that maps all component ids to the resolution
     */
    void fluidPot(Molecule &mi, Molecule &mj, ParaStrm &params, ParaStrm &paramInv, double * drm, double &Upot6LJ,
                         double &UpotXpoles,
                         double &MyRF, bool calculateLJ, bool noHybrid,
                         const std::vector<Resolution::ResolutionType> &compResMap, const Resolution::FPRegion &region);

private:
    using PP2PFAThreadData = ParticlePairs2PotForceAdapter::PP2PFAThreadData;
    //! @brief meso-buffer and params for each thread
    std::vector<PP2PFAThreadData*> _threadData;
    //! @brief Reference to AdResS submodule responsible for resolution management
    Resolution::Handler& _resolutionHandler;
    //! @brief Force Function from IBI
    FunctionPL _ibiForce;
    //! @brief Potential Function from IBI
    FunctionPL _ibiPot;
    //! @brief flag to enable _ibiForce and _ibiPot
    bool _useIBIFunctions;

	/**
 * Simple container for mesoscopic values.
 * During normal force computation of ls1-mardyn mesoscopic values are computed. In Hybrid regions those are wrong.
 * AdResS needs to recompute those. This struct stores those values.
 * */
	struct MesoValues {
		//! @brief variable used to sum the virial contribution of all pairs
		double _virial;
		//! @brief variable used to sum the Upot6LJ contribution of all pairs
		double _upot6LJ;
		//! @brief variable used to sum the UpotXpoles contribution of all pairs
		double _upotXpoles;
		//! @brief variable used to sum the MyRF contribution of all pairs
		double _myRF;

		/**
		 * @brief Constructs a MesoValues Container
		 * @param v virial
		 * @param lj upot6LJ
		 * @param pole upotXpoles
		 * @param rf myRF
		 * */
		explicit MesoValues(double v = 0.f, double lj = 0.f, double pole = 0.f, double rf = 0.f) : _virial(v), _upot6LJ(lj), _upotXpoles(pole), _myRF(rf) {}

		/**
		 * @brief sets all values to 0
		 * */
		void clear() {
			_virial = 0;
			_upot6LJ = 0;
			_upotXpoles = 0;
			_myRF = 0;
		}

		/**
		 * @brief Sets the mesoscopic values in the passed @param domain.
		 * */
		void setInDomain(Domain* domain) const {
			domain->setLocalUpot(_upot6LJ / 6. + _upotXpoles + _myRF + domain->getLocalUpot());
			domain->setLocalVirial(_virial + 3.0 * _myRF + domain->getLocalVirial());
		}
	} _mesoValues;

    /** @brief Computes potential and force where both molecules are hybrid molecules.
     *
     * Calculate potential and force between two molecules including all site-site interactions.
     *
     * Calculates the potential energy and force between two molecules i and j.
     * The interaction paramters are held in precomputed streams, which are initialized in Comp2Param::initialize
     *
     * @param[in]  mi   molecule i
     * @param[in]  mj   molecule j
     * @param[in]  params    reference to the corresponding interaction parameter stream initialized via Comp2Param::initialize
     * @param[in]  drm   distance vector from molecule j to molecule i
     * @param[out] Upot6LJ   potential energy resulting from Lennard Jones interactions
     * @param[out] UpotXpoles   potential energy resulting from Charge, Dipole and Quadrupole interactions
     * @param[out] MyRF
     * @todo Document parameter (is this reaction field?)
     * @param[out] Virial   Virial
     * @param[in]  caculateLJ    enable or disable calculation of Lennard Jones interactions
     */
    void
    potForceFullHybrid(Molecule &mi, Molecule &mj, ParaStrm &params, double * drm, double &Upot6LJ,
                       double &UpotXpoles, double &MyRF, double Virial[3], bool calculateLJ, const Resolution::FPRegion &region);
    /** @brief Computes potential and force where only one molecule is hybrid.
     * Calculate potential and force between two molecules including all site-site interactions.
     *
     * Calculates the potential energy and force between two molecules i and j.
     * The interaction paramters are held in precomputed streams, which are initialized in Comp2Param::initialize
     *
     * @param[in]  mi   molecule i -> hybrid molecule
     * @param[in]  mj   molecule j
     * @param[in]  params    reference to the corresponding interaction parameter stream initialized via Comp2Param::initialize
     * @param[in]  drm   distance vector from molecule j to molecule i
     * @param[out] Upot6LJ   potential energy resulting from Lennard Jones interactions
     * @param[out] UpotXpoles   potential energy resulting from Charge, Dipole and Quadrupole interactions
     * @param[out] MyRF
     * @todo Document parameter (is this reaction field?)
     * @param[out] Virial   Virial
     * @param[in]  caculateLJ    enable or disable calculation of Lennard Jones interactions
     */
    void
    potForceSingleHybrid(Molecule &mi, Molecule &mj, ParaStrm &params, double * drm, double &Upot6LJ,
                         double &UpotXpoles, double &MyRF, double Virial[3], bool calculateLJ, const Resolution::FPRegion &region,
						 Resolution::ResolutionType resolutionJ);

    /** @brief Computes potential where both molecules are hybrid molecules.
     *
     * Calculate potential and force between two molecules including all site-site interactions.
     *
     * Calculates the potential energy and force between two molecules i and j.
     * The interaction paramters are held in precomputed streams, which are initialized in Comp2Param::initialize
     *
     * @param[in]  mi   molecule i
     * @param[in]  mj   molecule j
     * @param[in]  params    reference to the corresponding interaction parameter stream initialized via Comp2Param::initialize
     * @param[in]  drm   distance vector from molecule j to molecule i
     * @param[out] Upot6LJ   potential energy resulting from Lennard Jones interactions
     * @param[out] UpotXpoles   potential energy resulting from Charge, Dipole and Quadrupole interactions
     * @param[out] MyRF
     * @todo Document parameter (is this reaction field?)
     * @param[out] Virial   Virial
     * @param[in]  caculateLJ    enable or disable calculation of Lennard Jones interactions
     */
    void
    fluidPotFullHybrid(Molecule &mi, Molecule &mj, ParaStrm &params, double * drm, double &Upot6LJ,
                       double &UpotXpoles, double &MyRF, bool calculateLJ, const Resolution::FPRegion &region);
    /** @brief Computes potential where only one molecule is hybrid.
     * Calculate potential and force between two molecules including all site-site interactions.
     *
     * Calculates the potential energy and force between two molecules i and j.
     * The interaction paramters are held in precomputed streams, which are initialized in Comp2Param::initialize
     *
     * @param[in]  mi   molecule i -> hybrid molecule
     * @param[in]  mj   molecule j
     * @param[in]  params    reference to the corresponding interaction parameter stream initialized via Comp2Param::initialize
     * @param[in]  drm   distance vector from molecule j to molecule i
     * @param[out] Upot6LJ   potential energy resulting from Lennard Jones interactions
     * @param[out] UpotXpoles   potential energy resulting from Charge, Dipole and Quadrupole interactions
     * @param[out] MyRF
     * @todo Document parameter (is this reaction field?)
     * @param[out] Virial   Virial
     * @param[in]  caculateLJ    enable or disable calculation of Lennard Jones interactions
     */
    void
    fluidPotSingleHybrid(Molecule &mi, Molecule &mj, ParaStrm &params, double * drm, double &Upot6LJ,
                         double &UpotXpoles, double &MyRF, bool calculateLJ, const Resolution::FPRegion &region,
						 Resolution::ResolutionType resolutionJ);

    void PotForceIBI(Molecule& mi, Molecule& mj, ParaStrm& params, double drm[3], double& Upot6LJ, double& UpotXpoles, double& MyRF, double Virial[3]);

    void FluidPotIBI(Molecule& mi, Molecule& mj, ParaStrm& params, double /*drm*/[3], double& Upot6LJ, double& UpotXpoles, double& MyRF);
};


#endif //MARDYN_ADRESSFORCEADAPTER_H
