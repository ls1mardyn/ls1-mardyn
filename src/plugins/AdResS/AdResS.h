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
#include "AdResSForceAdapter.h"
#include "DensityProfile3D.h"
#include "Interpolation.h"

class ParticleContainer;
class DomainDecompBase;
class Domain;
class AdResSForceAdapter;
class AdResSForceAdapterTest;
class AdResSTest;
class AdResSKDDecomposition;

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
    friend class AdResSForceAdapter;
    friend class AdResSKDDecomposition;
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

    void afterForces(ParticleContainer *container, DomainDecompBase *base, unsigned long i) override;

    std::string getPluginName() override;

    /**
     * @brief Creates a new instance of this plugin. Used in the PluginFactory.
     * */
    static PluginBase* createInstance() {
        return new AdResS();
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
static double (*weight)(std::array<double, 3> r, FPRegion& region);

private:
    /**
     * @brief Container for mesoscopic values during force recomputation.
     * */
    MesoValues _mesoVals;

    /**
     * Handles force computation
     * */
    AdResSForceAdapter* _forceAdapter;

    /**
     * @brief Container for all areas of interest for AdResS. The ares are currently boxes.
     * */
    std::vector<FPRegion> _fpRegions;

    /**
     * pointer to particle container, is set up during AdResS::init
     * */
    ParticleContainer* _particleContainer;

    /**
     * Ptr to all saved components in the current simulation.
     * This resource is owned by the active domain instance.
     * */
    std::vector<Component>* _components;

    /**
     * Maps each component id to its resolution for faster component switching.
     * */
    std::vector<Resolution> _comp_to_res;

    //! @brief reference to the domain is needed to store the calculated macroscopic values
    Domain* _domain;

    //! @brief True iff F_th functionality should be used
    bool _enableThermodynamicForce;

    //! @brief True iff F_th must be sampled first
    bool _createThermodynamicForce;

    //! @brief Simulation iterations between each F_th recalculation
    int _thermodynamicForceSampleGap;

    //! @brief Iterations since last sampling
    int _thermodynamicForceSampleCounter;

    //! @brief Class to create 3D density profiles
    DensityProfile3D _densityProfiler;

    //! @brief Thermodynamic force used to correct the density difference created by plain AdResS
    Interpolation::Function _thermodynamicForce;
    Interpolation::Function _thermodynamicForceHist;

    //! @brief Gradient of density distribution, used for convergence checking
    Interpolation::Function _lastGradient;

    //! @brief Density function of a FP simulation, used to find convergence of F_th
    std::vector<double> _targetDensity;

    //! @brief Describes a accuracy measure in terms of: max|rho(x)-rho_target(x)|/rho_target(x) <= threshold
    //! used to find convergence of F_th
    double _convergenceThreshold;

    //! @brief Controls speed of convergence for F_th
    double _convergenceFactor;

    //! @brief step size of density sampling
    double _samplingStepSize;

    //! @brief enables logging of FTH to file
    bool _logFTH;

    //! @brief enables logging of the current simulation densities to file
    bool _logDensities;

    //! @brief maximum allowed force
    double _forceMax;

    //! @brief initial density
    double _rho0;

    //! @brief interpolation smoothness
    double _smoothingFactor;

    /**
     * Writes the function object into the required XML format for input files.
     * @param filename Output filename
     * @param fun function object
     * */
    void writeFunctionToXML(const std::string& filename, Interpolation::Function& fun);

    /**
     * Writes the provided densities as a sequence split by separator into a file.
     * */
    void writeDensities(const std::string& filename, std::vector<double>& densities, const std::string& separator = " ");

    /**
     * Recomputes F_th in _thermodynamicForce by sampling densities and interpolating the gradient.
     * According to: F_k+1(x) = F_k(x) - c * d'(x)
     * */
    void computeF_TH();

    /**
     * Checks if the current simulation density is at most _convergenceThreshold apart from _targetDensity.
     * @returns true iff converged
     * */
    bool checkF_TH_Convergence();

    /**
     * Applies the thermodynamic force to ?all? molecules.
     * */
    void applyF_TH();

    /**
     * @brief checks the component of @param molecule and sets it to the correct LOD depending the @param targetRes.
     * */
    void checkMoleculeLOD(Molecule& molecule, Resolution targetRes);

    /**
     * Computes all forces. Depending on "invert" all forces are applied with AdResS concepts or inverted with no weight.
     *
     * Invert = true:
     * During the traversal of all cells, the default cell processor handles the cells without knowledge of AdResS.
     * Therefore, all forces and mesoscopic values computed for FPRegions are incorrect.
     * This method removes all force and meso contributions all FPRegions made.
     * */
    void computeForce(bool invert);

    /**
     * @brief Weighting function for AdResS force computation.
     * Implementation computes axis intersection points and uses euclidean distance to determine the period of
     * the underlying cosine function.
     * @param r position of the site
     * @param region region of FP
     * */
    static double weightEuclid(std::array<double, 3> r, FPRegion& region);

    /**
     * @brief Weighting function for AdResS force computation.
     * Implementation computes axis intersection points and uses manhattan distance to determine the period of
     * the underlying cosine function.
     * @param r position of the site
     * @param region region of FP
     * */
    static double weightManhattan(std::array<double, 3> r, FPRegion& region);

    /**
     * @brief Weighting function for AdResS force computation.
     * Implementation computes weight percentage for each component by checking where each component is in the hybrid region.
     * All component weights are multiplied together.
     * This results in a weight function, where each component does its own contribution.
     * @param r position of the site
     * @param region region of FP
     * */
    static double weightComponent(std::array<double, 3> r, FPRegion& region);

    /**
     * @brief Weighting function for AdResS force computation.
     * Implementation conceptually find the nearest point on the surface of the inner region in respect to r and
     * computes the distance between r and this point.
     * By doing so, this weight function treats the FPRegion as a box with rounded corners and edges.
     * @param r position of the site
     * @param region region of FP
     * */
    static double weightNearest(std::array<double, 3> r, FPRegion& region);

    /**
     * @brief Weighting function for AdResS force computation.
     * Implementation disables Hybrid region. Only for testing purposes, do not use in production.
     * @param r position of the site
     * @param region region of FP
     * */
    static double weightFlat(std::array<double, 3> r, FPRegion& region);
};

#endif //MARDYN_ADRESS_H

