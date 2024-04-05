//
// Created by alex on 27.04.23.
//

#include "AdResS.h"

#include "Simulation.h"
#include "Domain.h"
#include "utils/mardyn_assert.h"
#include "ensemble/EnsembleBase.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"

#include "plugins/AdResS/util/AdResSRegionTraversal.h"
#include "plugins/AdResS/parallel/AdResSKDDecomposition.h"
#include "Interpolation.h"

#include <cmath>

using namespace std;
using Log::global_log;

Weight::function_t AdResS::weight = nullptr;

AdResS::AdResS() : _resolutionHandler(), _forceAdapter(nullptr) {};

AdResS::~AdResS() = default;

void AdResS::init(ParticleContainer*, DomainDecompBase*, Domain*) {
    global_log->debug() << "[AdResS] Enabled " << std::endl;
}

void AdResS::readXML(XMLfileUnits &xmlconfig) {
    static const std::array<std::string, 5> weight_impls {"euclid", "manhattan", "component", "near", "flat"};
    std::string impl = weight_impls[0];
    xmlconfig.getNodeValue("weightImpl", impl);
    unsigned long index = std::distance(weight_impls.begin(), std::find(weight_impls.begin(), weight_impls.end(), impl));
    AdResS::weight = Weight::functions[index >= 5 ? 0 : index];
    global_log->info() << "[AdResS] Using weight implementation " << weight_impls[index >= 5 ? 0 : index] << std::endl;

	FTH::Config fthConf { };
    XMLfile::Query query = xmlconfig.query("enableFTH");
    fthConf._enableThermodynamicForce = false;
    unsigned long count = query.card();
    if (count > 1) {
        global_log->fatal() << "[AdResS] Only specify one enableFTH block in config file!" << std::endl;
        Simulation::exit(668);
    }
    //F_TH enabled
    if(count == 1) {
		fthConf._enableThermodynamicForce = true;
		fthConf._logFTH = xmlconfig.getNodeValue_bool("enableFTH/logFTH", false);
		fthConf._logDensities = xmlconfig.getNodeValue_bool("enableFTH/logDensity", false);
		fthConf._rho0 = xmlconfig.getNodeValue_double("enableFTH/rho", 0.01);
		fthConf._smoothingFactor = xmlconfig.getNodeValue_double("enableFTH/smoothingFactor", 10);

        query = xmlconfig.query("enableFTH/createFTH");
        count = query.card();
        if (count > 1) {
            global_log->fatal() << "[AdResS] Only specify one createFTH block in config file!" << std::endl;
            Simulation::exit(668);
        }

        if(count == 1) { // sample FTH function
            fthConf._createThermodynamicForce = true;
            fthConf._thermodynamicForceSampleGap = xmlconfig.getNodeValue_int("enableFTH/createFTH/sampleGap", 200);
            fthConf._convergenceThreshold = xmlconfig.getNodeValue_double("enableFTH/createFTH/threshold", 0.02);
            fthConf._convergenceFactor = xmlconfig.getNodeValue_double("enableFTH/createFTH/convFactor", 0.2);
            fthConf._samplingStepSize = xmlconfig.getNodeValue_double("enableFTH/createFTH/sampleBinSize", 0.2);

        }
        else { // use existing FTH function
            query = xmlconfig.query("enableFTH/forceFunction");
            count = query.card();
            if (count != 1) {
                global_log->fatal() << "[AdResS] Must specify one forceFunction block in config file!" << std::endl;
                Simulation::exit(668);
            }
			fthConf._logFTH = false;
			fthConf._thermodynamicForce.begin = xmlconfig.getNodeValue_double("enableFTH/forceFunction/startX");
            query = xmlconfig.query("enableFTH/forceFunction/samplePoint");
            count = query.card();
            if(count < 5) {
                global_log->fatal() << "[AdResS] Force function must have at least 5 sample points!" << std::endl;
                Simulation::exit(668);
            }

            fthConf._thermodynamicForce.n = count;
            fthConf._thermodynamicForce.gradients.resize(count, 0.0);
            fthConf._thermodynamicForce.function_values.resize(count, 0.0);
            fthConf._thermodynamicForce.step_width.resize(count-1, 0.0);
            XMLfile::Query::const_iterator sampleIter;
            std::string oldpath = xmlconfig.getcurrentnodepath();
            for(sampleIter = query.begin(); sampleIter; sampleIter++) {
                xmlconfig.changecurrentnode(sampleIter);
                unsigned long id = 0;
                xmlconfig.getNodeValue("@id", id);
                xmlconfig.getNodeValue("grad", fthConf._thermodynamicForce.gradients[id - 1]);
                xmlconfig.getNodeValue("func", fthConf._thermodynamicForce.function_values[id - 1]);
                if(id < count) xmlconfig.getNodeValue("step", fthConf._thermodynamicForce.step_width[id - 1]);
            }
            xmlconfig.changecurrentnode(oldpath);
			fthConf._samplingStepSize = fthConf._thermodynamicForce.step_width[0];
			fthConf._thermodynamicForceSampleGap = 200;
			fthConf._createThermodynamicForce = false;
        }
    }
	_fthHandler = FTH::Handler(fthConf);

	Resolution::Config resConf { };
	resConf.components = _simulation.getEnsemble()->getComponents();
	resConf.domain = _simulation.getDomain();
	resConf.comp_to_res.resize(resConf.components->size(), Resolution::FullParticle);

    long numRegions = 0;
    query = xmlconfig.query("fpregions/region");
    numRegions = query.card();
    if (numRegions == 0) {
        global_log->fatal() << "No AdResS regions specified: Falling back to dynamic region selection. ERROR: not implemented yet!" << std::endl;
        Simulation::exit(668);
    }
    if (numRegions != 1) {
        global_log->fatal() << "AdResS currently only supports a single FP Region. For more general support select a previous version." << std::endl;
        Simulation::exit(-1);
    }
    resConf.fpRegions.resize(numRegions);

    XMLfile::Query::const_iterator regionIter;
    std::string oldpath = xmlconfig.getcurrentnodepath();
    for(regionIter = query.begin(); regionIter; regionIter++) {
        xmlconfig.changecurrentnode(regionIter);
        unsigned int id = 0;
        xmlconfig.getNodeValue("@id", id);
        resConf.fpRegions[id - 1].readXML(xmlconfig);
    }
    xmlconfig.changecurrentnode(oldpath);
	_resolutionHandler = Resolution::Handler(resConf);
    // todo add check that no region overlap even in hybrid considering periodic bounds

    _forceAdapter = new AdResSForceAdapter(_resolutionHandler);
    _simulation.setParticlePairsHandler(_forceAdapter);
    _simulation.setCellProcessor(new LegacyCellProcessor(_simulation.getcutoffRadius(), _simulation.getLJCutoff(), _forceAdapter));
    if(auto decomp = dynamic_cast<AdResSKDDecomposition*>(&_simulation.domainDecomposition())) {
        decomp->setResolutionHandler(_resolutionHandler);
    }
}

void AdResS::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                     unsigned long simstep) {
	_fthHandler.writeLogs(*particleContainer, *domainDecomp, *domain, simstep);
}

void AdResS::finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {
	_fthHandler.writeLogs(*particleContainer, *domainDecomp, *domain, _simulation.getNumTimesteps());
}

std::string AdResS::getPluginName() {
    return {"AdResS"};
}

void AdResS::beforeForces(ParticleContainer *container, DomainDecompBase *, unsigned long) {
    _resolutionHandler.checkResolution(*container);
	_fthHandler.step(*container, _resolutionHandler.getRegions());
}

void AdResS::siteWiseForces(ParticleContainer *container, DomainDecompBase *base, unsigned long i) {
    _fthHandler.apply(*container, _resolutionHandler.getRegions(), _resolutionHandler.getCompResMap());
}
