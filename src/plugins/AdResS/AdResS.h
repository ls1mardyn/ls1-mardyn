//
// Created by Alex Hocks on 27.04.23.
//

#ifndef MARDYN_ADRESS_H
#define MARDYN_ADRESS_H

#include <memory>
#include "plugins/PluginBase.h"
#include "utils/Logger.h"
#include "particleContainer/RegionParticleIterator.h"
#include "molecules/potforce.h"
#include "AdResSData.h"

class ParticleContainer;
class DomainDecompBase;
class Domain;

/**
 * With the introduction of AdResS to ls1-MarDyn the simulation class owns two particle containers,
 * each of which handle their own respective domain.
 * Simulation::_moleculeContainerAlt contains particles in a CG resolution.
 * This plugin moves and transforms molecules upon entering a domain of another resolution.
 * Currently supports only split domain along x-Axis
 *
 * BEWARE:
 * When creating components for AdResS,
 * the full particle component must have id n,
 * the hybrid component id n+1 and the CG component id n+2!
 * The names of the components must be prefixed with "FP_", "H_" or "CG_"!
 * */
class AdResS : public PluginBase {
public:
    /**
     * Constructor, no params needed
     * */
    AdResS();

    /**
     * Destructor for inheritance
     * */
    ~AdResS() override;

    /**
     * Initializes local fields and checks if all components have the correct format
     * */
    void init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) override;

    /**
     * @brief Read in XML configuration for AdResS
     *
     * The following XML object structure is handled by this method:
     * \code{.xml}
     *  <plugin name="AdResS">
     *      <fpregions>
     *          <region id="1">
     *              <lowX>0.0</lowX><lowY>0.0</lowY><lowZ>0.0</lowZ>
     *              <highX>0.0</highX><highY>0.0</highY><highZ>0.0</highZ>
     *              <hybridDim>10.0</hybridDim>
     *          </region>
     *      </fpregions>
     *  </plugin>
     * \endcode
     *
     * fpregions can contain N region objects. The first region object must have id 1.
     * The id must grow in increasing order. If N regions are used, then the largest id must be N and the smallest 1.
     * */
    void readXML(XMLfileUnits &xmlconfig) override;

    /**
     * Does nothing
     * */
    void endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                 unsigned long simstep) override;

    /**
     * Does nothing
     * */
    void finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) override;

    /**
     * @brief Performs the shifting from one LOD to another. Currently exchanges the components of the molecules.
     * @param container particle container
     * @param base todo Alex ... idk?
     * @param simstep current simulation step
     * */
    void beforeForces(ParticleContainer *container, DomainDecompBase *base, unsigned long simstep) override;

    /**
     * @brief computes the forces in the direct connection areas between the different particle containers
     * */
    void siteWiseForces(ParticleContainer *container, DomainDecompBase *base, unsigned long i) override;

    std::string getPluginName() override;

    /**
     * @brief Creates a new instance of this plugin. Used in the PluginFactory.
     * */
    static PluginBase* createInstance() {
        return new AdResS();
    }

private:
    /**
     * @brief Container for mesoscopic values during force recomputation.
     * */
    MesoValues _mesoVals;

    /**
     * @brief Container for all areas of interest for AdResS. The ares are currently boxes.
     * */
    std::vector<FPRegion> _fpRegions;

    /**
     * pointer to particle container, is set up during AdResS::init
     * */
    ParticleContainer* _particleContainer;

    /**
     * Ptr to all saved components in the current simulation
     * */
    std::vector<Component>* _components;

    /**
     * Maps each component id to its resolution for faster component switching.
     * */
    std::unordered_map<unsigned long, Resolution> _comp_to_res;

    //! @brief reference to the domain is needed to store the calculated macroscopic values
    Domain* _domain;

    /**
     * @brief checks the component of @param molecule and sets it to the correct LOD depending the @param targetRes.
     * */
    void checkMoleculeLOD(Molecule& molecule, Resolution targetRes);

    /**
     * @brief Weighting function for AdResS force computation.
     * For smooth transition from full particle to coarse grain area.
     * Hybrid area is a wall surrounding full particle area. The weight is then based on the position of the
     * viewed molecules site along the axis from the molecule with the shortest distance to the other area.
     * @param r position of the site
     * @param region region of FP
     * */
    static double weight(std::array<double, 3> r, FPRegion& region);
};

#endif //MARDYN_ADRESS_H

