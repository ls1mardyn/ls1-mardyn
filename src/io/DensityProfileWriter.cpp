#include "io/DensityProfileWriter.h"

#include <iomanip>

#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"

#include "utils/FileUtils.h"

void DensityProfileWriter::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "[DensityProfileWriter] Write frequency: " << _writeFrequency << endl;
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "[DensityProfileWriter] Output prefix: " << _outputPrefix << endl;

	int doRecordVirialProfile = 0;
	xmlconfig.getNodeValue("options/option@[keyword='profileVirial']", doRecordVirialProfile);
	_doRecordVirialProfile = (doRecordVirialProfile > 0) ? true : false;
	global_log->info() << "[DensityProfileWriter] Record Virial: " << _doRecordVirialProfile << endl;

	xmlconfig.getNodeValue("timesteps/init", _initStatistics);
	global_log->info() << "[DensityProfileWriter] init statistics: " << _initStatistics << endl;
	xmlconfig.getNodeValue("timesteps/recording", _profileRecordingTimesteps);
	global_log->info() << "[DensityProfileWriter] profile recording timesteps: " << _profileRecordingTimesteps << endl;
}

void DensityProfileWriter::init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {
	
}

void DensityProfileWriter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                                   unsigned long simstep, std::list<ChemicalPotential> *lmu,
                                   std::map<unsigned int, CavityEnsemble> *mcav) {
	int mpi_rank = domainDecomp->getRank();
	if ((simstep >= _initStatistics) && (simstep % _profileRecordingTimesteps == 0)) {
		domain->recordProfile(particleContainer, _doRecordVirialProfile);
	}
	if ((simstep >= _initStatistics) && (simstep % _writeFrequency == 0)) {
		domain->collectProfile(domainDecomp, _doRecordVirialProfile);
		if (mpi_rank == 0) {
			ostringstream osstrm;
			osstrm << _outputPrefix << "." << fill_width('0', 9) << simstep;
			//edited by Michaela Heier
			if(domain->isCylindrical()){
				domain->outputCylProfile(osstrm.str().c_str(),_doRecordVirialProfile);
			}
			else{
				domain->outputProfile(osstrm.str().c_str(), _doRecordVirialProfile);
			}
			osstrm.str("");
			osstrm.clear();
		}
		domain->resetProfile(_doRecordVirialProfile);
	}
}

void DensityProfileWriter::finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {
	
}
