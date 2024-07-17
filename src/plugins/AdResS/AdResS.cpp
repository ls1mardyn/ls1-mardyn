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

AdResS::AdResS() : _resolutionHandler(), _fthHandler(nullptr), _forceAdapter(nullptr), _density_sampler(nullptr), _grid(nullptr), _samplingGap(100), _samlingCounter(0) {};

AdResS::~AdResS() {
	delete _fthHandler;
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
	_samplingGap = xmlconfig.getNodeValue_int("sampleGap", 100);

    if (sample_type == "Grid") {
        global_log->info() << "[AdResS] Using sampler implementation " << sample_type << std::endl;
        int x = xmlconfig.getNodeValue_int("gridSize/nx",1);
        int y = xmlconfig.getNodeValue_int("gridSize/ny",1);
        int z = xmlconfig.getNodeValue_int("gridSize/nz",1);
        global_log->info() << "[AdResS] Grid size is" << x << " " << y << " " << z << std::endl;

		_grid = new FTH::grid_t();
        // we want the grid to also extend to the bounds, so we do not have to manage grid information exchange
        const double r_cutoff = _simulation.getcutoffRadius();
		_grid->init({_simulation.domainDecomposition().getBoundingBoxMin(0, _simulation.getDomain()) - r_cutoff,
					 _simulation.domainDecomposition().getBoundingBoxMin(1, _simulation.getDomain()) - r_cutoff,
					 _simulation.domainDecomposition().getBoundingBoxMin(2, _simulation.getDomain()) - r_cutoff},
					{_simulation.domainDecomposition().getBoundingBoxMax(0, _simulation.getDomain()) + r_cutoff,
					 _simulation.domainDecomposition().getBoundingBoxMax(1, _simulation.getDomain()) + r_cutoff,
					 _simulation.domainDecomposition().getBoundingBoxMax(2, _simulation.getDomain()) + r_cutoff},
					 x, y, z);

        double rad = xmlconfig.getNodeValue_double("sampleRadius", 1.0);
        _density_sampler = new AveragedGridSampler{_grid, rad};
		_fthHandler = new FTH::Grid3DHandler();
    } else if (sample_type == "Direct") {
		global_log->info() << "[AdResS] Using sampler implementation " << sample_type << std::endl;
		int dim = xmlconfig.getNodeValue_int("dim", 0);
		double bin_width = xmlconfig.getNodeValue_double("binWidth", 1);

		_density_sampler = new DirectProjectedSampler(dim, bin_width);
		_fthHandler = new FTH::Grid1DHandler();
	} else if (sample_type == "Smooth") {
		global_log->info() << "[AdResS] Using sampler implementation " << sample_type << std::endl;
		int dim = xmlconfig.getNodeValue_int("dim", 0);
		double bin_width = xmlconfig.getNodeValue_double("binWidth", 1);
		double smoothingStrength = xmlconfig.getNodeValue_double("smoothingStrength", 1);

		_density_sampler = new SmoothedProjectedSampler(dim, bin_width, smoothingStrength);
		_fthHandler = new FTH::Grid1DHandler();
	} else if (sample_type == "FT") {
		global_log->info() << "[AdResS] Using sampler implementation " << sample_type << std::endl;
		int dim = xmlconfig.getNodeValue_int("dim", 0);
		int frequencies = xmlconfig.getNodeValue_int("frequencies", 500);
		int samples = xmlconfig.getNodeValue_int("samples", 100);

		_density_sampler = new FTProjectedSampler(dim, frequencies, samples);
		_fthHandler = new FTH::Function1DHandler();
	} else if (sample_type == "GMM") {
		global_log->info() << "[AdResS] Using sampler implementation " << sample_type << std::endl;
		int dim = xmlconfig.getNodeValue_int("dim", 0);
		int samples = xmlconfig.getNodeValue_int("samples", 100);
		double smoothingStrength = xmlconfig.getNodeValue_double("smoothingStrength", 1);

		_density_sampler = new GMMProjectedSampler(dim, samples, smoothingStrength);
		_fthHandler = new FTH::Function1DHandler();
	} else {
		global_log->info() << "[AdResS] Falling back to Direct sampler implementation " << std::endl;
		int dim = 0;
		double bin_width = 1;

		_density_sampler = new DirectProjectedSampler(dim, bin_width);
		_fthHandler = new FTH::Grid1DHandler();
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
        	fthConf._fth_file_path = xmlconfig.getNodeValue_string("enableFTH/pathFTH", ""); // reload from file and keep refining
        }
        else { // use existing FTH function
			fthConf._fth_file_path = xmlconfig.getNodeValue_string("enableFTH/pathFTH", "");
			fthConf._logFTH = false;
            fthConf._thermodynamicForceSampleGap = 200;
			fthConf._createThermodynamicForce = false;
        }
    }
	fthConf._grid = _grid;
	fthConf._density_sampler = _density_sampler;
	_fthHandler->init(fthConf);

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
	_fthHandler->writeLogs(*particleContainer, *domainDecomp, *domain, simstep);
}

void AdResS::finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {
	_fthHandler->writeFinalFTH();
}

std::string AdResS::getPluginName() {
    return {"AdResS"};
}

void AdResS::beforeForces(ParticleContainer *container, DomainDecompBase *dd, unsigned long) {
	_samlingCounter++;
	if ((_samlingCounter % _samplingGap) == 0) {
		_density_sampler->sampleData(container, dd, _simulation.getDomain());
		_samlingCounter = 0;
	}
    _resolutionHandler.checkResolution(*container);
	_fthHandler->computeSingleIteration(*container, _resolutionHandler.getRegions());
}

void AdResS::siteWiseForces(ParticleContainer *container, DomainDecompBase *base, unsigned long i) {
	_fthHandler->applyForce(*container, _resolutionHandler.getRegions(), _resolutionHandler.getCompResMap());
}
