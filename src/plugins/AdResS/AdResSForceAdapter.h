//
// Created by alex on 09.05.23.
//

#ifndef MARDYN_ADRESSFORCEADAPTER_H
#define MARDYN_ADRESSFORCEADAPTER_H

#include "molecules/Comp2Param.h"
#include "AdResSData.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"

/**
 * Handles all force calculation done in AdResS method.
 * */
class AdResSForceAdapter {
public:
    /**
     * Constructor
     * @param mesoValues reference to the AdResS::_mesoVals instance
     * */
    explicit AdResSForceAdapter(MesoValues& mesoValues);

    /**
     * Destructor: clears allocated memory
     * */
    ~AdResSForceAdapter() noexcept;

    /**
     * @brief initializes the local Comp2Param instance by using the provided domain.
     * _c2p is used to force calculation
     * */
    void init(Domain* domain);

    /**
     * @brief reduction of all thread local buffers, called when all forces are computed.
     * Will set results in _mesoValues.
     * */
    void finish();

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
     * @param invert          true if inverted forces, potentials should be applied and weights ignored
     * @param compResMap      a map from AdResS plugin instance that maps all component ids to the resolution
     *
     * @return                interaction energy
     */
    double
    processPair(Molecule &molecule1, Molecule &molecule2, std::array<double, 3> distanceVector, PairType pairType,
                double dd, bool calculateLJ, bool invert, std::unordered_map<unsigned long, Resolution> &compResMap,
                FPRegion &region);

    /** @brief Based on potforce.h::PotForce
     *  The main change is the addition of the parameter invert.
     *  When it is true, then the computation is inverted and neglects hybrid weights to undo prior force computation.
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
     * @param[in] invert if true undo the non-AdResS force calculation
     * @param[in] compResMap      a map from AdResS plugin instance that maps all component ids to the resolution
     */
    static void potForce(Molecule &mi, Molecule &mj, ParaStrm &params, std::array<double, 3> drm, double &Upot6LJ,
                         double &UpotXpoles,
                         double &MyRF, double Virial[3], bool calculateLJ, bool invert,
                         std::unordered_map<unsigned long, Resolution> &compResMap, FPRegion &region);

private:
    using PP2PFAThreadData = ParticlePairs2PotForceAdapter::PP2PFAThreadData;
    //! @brief mesoscopic value buffer
    MesoValues& _mesoValues;
    //! @brief meso-buffer and params for each thread
    std::vector<PP2PFAThreadData*> _threadData;

    /** @brief Computes the inverted potential and force to reverse prior forces. Does not consider weight function.
     *
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
    static void inline invertedPotForce(Molecule& mi, Molecule& mj, ParaStrm& params, std::array<double,3> drm,
                                        double& Upot6LJ, double& UpotXpoles, double& MyRF, double Virial[3],
                                        bool calculateLJ);

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
    static void
    potForceFullHybrid(Molecule &mi, Molecule &mj, ParaStrm &params, std::array<double, 3> drm, double &Upot6LJ,
                       double &UpotXpoles, double &MyRF, double Virial[3], bool calculateLJ, FPRegion &region);
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
    static void
    potForceSingleHybrid(Molecule &mi, Molecule &mj, ParaStrm &params, std::array<double, 3> drm, double &Upot6LJ,
                         double &UpotXpoles, double &MyRF, double Virial[3], bool calculateLJ, FPRegion &region,
                         Resolution resolutionJ);
};


#endif //MARDYN_ADRESSFORCEADAPTER_H
