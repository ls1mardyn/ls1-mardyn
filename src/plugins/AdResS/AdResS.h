//
// Created by Alex Hocks on 27.04.23.
//

#ifndef MARDYN_ADRESS_H
#define MARDYN_ADRESS_H

#include "plugins/PluginBase.h"

#include "plugins/AdResS/adapters/AdResSForceAdapter.h"
#include "features/Resolution.h"
#include "features/FTH.h"
#include "util/WeightFunction.h"

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
 * The hybrid component must first contain all CG parts with mass = 0 and then the FP part.
 * Regions cannot overlap, otherwise forces are computed multiple times.
 * */
class AdResS : public PluginBase {
    friend class AdResSForceAdapterTest;
    friend class AdResSTest;
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
     *      <weightImpl>euclid</weightImpl>
     *      <enableFTH>
     *          <!--createFTH>
     *              <sampleGap>100</sampleGap>
     *              <threshold>0.02</threshold>
     *              <convFactor>0.2</convFactor>
     *              <sampleBinSize>0.2</sampleBinSize>
     *              <logFTH>true</logFTH>
     *              <logDensity>true</logDensity>
     *          </createFTH-->
     *          <forceFunction>
     *              <sampleBinSize>100</sampleBinSize>
     *              <startX>20.0</startX>
     *              <logFTH>true</logFTH>
     *              <logDensity>true</logDensity>
     *              <samplePoint id="1">
     *                  <grad>2.0</grad>
     *                  <func>1.0</func>
     *              </samplePoint>
     *          </forceFunction>
     *      </enableFTH>
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
     * WeightImpl options: "euclid", "manhattan", "component", "near"; Default: "euclid"
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
        return dynamic_cast<PluginBase*>(new AdResS());
    }

	/**
 	* @brief Weighting function for AdResS force computation.
 	* For smooth transition from full particle to coarse grain area.
 	* Hybrid area is a wall surrounding full particle area. The weight is then based on the position of the
 	* viewed molecule along the axis from the molecules position to the center of the FPRegion.
 	* Depending on the initialization of the AdResS Plugin this will point to different weight function implementations.
 	* @param r position of the site
	* @param region region of FP
 	* */
	static Weight::function_t weight;

private:
	/**
	 * Handles all resolution aspects of AdResS
	 * */
	 Resolution::Handler _resolutionHandler;

	/**
	* Handles all FTH aspects of AdResS
	* */
	 FTH::Handler _fthHandler;

    /**
     * Handles force computation
     * */
    AdResSForceAdapter* _forceAdapter;
};

#endif //MARDYN_ADRESS_H

