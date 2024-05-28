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
#include "plugins/AdResS/adapters/VCPADR.h"
#include "Interpolation.h"

#include <cmath>

using namespace std;
using Log::global_log;

Weight::function_t AdResS::weight = nullptr;

AdResS::AdResS() : _resolutionHandler(), _forceAdapter(nullptr), _density_sampler(nullptr), _grid(nullptr) {};

AdResS::~AdResS() {
	delete _density_sampler;
	delete _grid;
}

void AdResS::init(ParticleContainer* pc, DomainDecompBase* dd, Domain* dom) {
    
    //sampler->init();
    //output mesh once at the beginning to a different file
    //ofstream mesh("mesh.txt");
    //mesh << grid;
    //mesh.close();
    global_log->debug() << "[AdResS] Enabled " << std::endl;
}

void AdResS::readXML(XMLfileUnits &xmlconfig) {
    static const std::array<std::string, 5> weight_impls {"euclid", "manhattan", "component", "near", "flat"};
    std::string impl = weight_impls[0];
    xmlconfig.getNodeValue("weightImpl", impl);
    unsigned long index = std::distance(weight_impls.begin(), std::find(weight_impls.begin(), weight_impls.end(), impl));
    AdResS::weight = Weight::functions[index >= 5 ? 0 : index];
    global_log->info() << "[AdResS] Using weight implementation " << weight_impls[index >= 5 ? 0 : index] << std::endl;

    //Check which type of sampler is selected
    std::string sample_type;
    xmlconfig.changecurrentnode("sampler");
    xmlconfig.getNodeValue("@type",sample_type);
    if (sample_type == "Grid") {
        global_log->info() << "[AdResS] Using sampler implementation " << sample_type << std::endl;
        int x = xmlconfig.getNodeValue_int("gridSize/xElements",1);
        int y = xmlconfig.getNodeValue_int("gridSize/yElements",1);
        int z = xmlconfig.getNodeValue_int("gridSize/zElements",1);
        global_log->info() << "[AdResS] Grid size is" << x<< std::endl;

		_grid = new FTH::grid_t();
		_grid->init({_simulation.domainDecomposition().getBoundingBoxMin(0, _simulation.getDomain()),
					 _simulation.domainDecomposition().getBoundingBoxMin(1, _simulation.getDomain()),
					 _simulation.domainDecomposition().getBoundingBoxMin(2, _simulation.getDomain())},
					{_simulation.domainDecomposition().getBoundingBoxMax(0, _simulation.getDomain()),
					 _simulation.domainDecomposition().getBoundingBoxMax(1, _simulation.getDomain()),
					 _simulation.domainDecomposition().getBoundingBoxMax(2, _simulation.getDomain())},
					 x, y, z);

        double rad = xmlconfig.getNodeValue_double("sampleRadius", 1.0);
        _density_sampler = new AveragedGridSampler{_grid, rad};
    } else if (sample_type == "Direct") {
		global_log->info() << "[AdResS] Using sampler implementation " << sample_type << std::endl;
		int dim = xmlconfig.getNodeValue_int("dim", 0);
		double bin_width = xmlconfig.getNodeValue_double("binWidth", 1);

		_density_sampler = new DirectProjectedSampler(dim, bin_width);
	} else if (sample_type == "Smooth") {
		global_log->info() << "[AdResS] Using sampler implementation " << sample_type << std::endl;
		int dim = xmlconfig.getNodeValue_int("dim", 0);
		double bin_width = xmlconfig.getNodeValue_double("binWidth", 1);
		double smoothingStrength = xmlconfig.getNodeValue_double("smoothingStrength", 1);

		_density_sampler = new SmoothedProjectedSampler(dim, bin_width, smoothingStrength);
	} else if (sample_type == "FT") {
		global_log->info() << "[AdResS] Using sampler implementation " << sample_type << std::endl;
		int dim = xmlconfig.getNodeValue_int("dim", 0);
		int frequencies = xmlconfig.getNodeValue_int("frequencies", 500);
		int samples = xmlconfig.getNodeValue_int("samples", 100);

		_density_sampler = new FTProjectedSampler(dim, frequencies, samples);
	} else if (sample_type == "GMM") {
		global_log->info() << "[AdResS] Using sampler implementation " << sample_type << std::endl;
		int dim = xmlconfig.getNodeValue_int("dim", 0);
		int samples = xmlconfig.getNodeValue_int("samples", 100);
		double smoothingStrength = xmlconfig.getNodeValue_double("smoothingStrength", 1);

		_density_sampler = new GMMProjectedSampler(dim, samples, smoothingStrength);
	} else {
		global_log->info() << "[AdResS] Falling back to Direct sampler implementation " << std::endl;
		int dim = 0;
		double bin_width = 1;

		_density_sampler = new DirectProjectedSampler(dim, bin_width);
	}
	_density_sampler->init(_simulation.getDomain());
	xmlconfig.changecurrentnode("..");

    //TODO:this mixes the config of 2 separate things, misleading
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
			fthConf._thermodynamicForceSampleGap = 200;
			fthConf._createThermodynamicForce = false;
        }
    }
	fthConf._grid = _grid;
	fthConf._density_sampler = _density_sampler;
	_fthHandler.init(fthConf);

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
	_resolutionHandler.init(resConf);
    // todo add check that no region overlap even in hybrid considering periodic bounds

	// we set the pair handler regardless, in case some other plugin requires pair wise interactions using AdResS
	_forceAdapter = new AdResSForceAdapter(_resolutionHandler);
	_simulation.setParticlePairsHandler(_forceAdapter);
	if(_simulation.usingLegacyCellProcessor()) {
		_simulation.setCellProcessor(new LegacyCellProcessor(_simulation.getcutoffRadius(), _simulation.getLJCutoff(), _forceAdapter));
	}
	else {
		_simulation.setCellProcessor(new VCPADR(*_simulation.getDomain(), _simulation.getcutoffRadius(), _simulation.getLJCutoff(), _resolutionHandler));
	}

    if(auto decomp = dynamic_cast<AdResSKDDecomposition*>(&_simulation.domainDecomposition())) {
        decomp->setResolutionHandler(_resolutionHandler);
    }
}

void AdResS::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                     unsigned long simstep) {
    
    //Sampler output 
    //std::string sampler_file{"density"+std::to_string(simstep)+".txt"};
    //sampler->writeSample(sampler_file);
	_fthHandler.writeLogs(*particleContainer, *domainDecomp, *domain, simstep);
}

void AdResS::finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {
	_fthHandler.writeFinalFTH();
}

std::string AdResS::getPluginName() {
    return {"AdResS"};
}

void AdResS::beforeForces(ParticleContainer *container, DomainDecompBase *, unsigned long) {
    //sampler->SampleData(container);
    _resolutionHandler.checkResolution(*container);
	_fthHandler.computeSingleIteration(*container, _resolutionHandler.getRegions());
}

void AdResS::siteWiseForces(ParticleContainer *container, DomainDecompBase *base, unsigned long i) {
	_fthHandler.applyForce(*container, _resolutionHandler.getRegions(), _resolutionHandler.getCompResMap());
}
