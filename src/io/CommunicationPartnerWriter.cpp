#include "io/CommunicationPartnerWriter.h"

#include <sstream>
#include <string>

#include "Common.h"
#include "Domain.h"
#include "Simulation.h"
#include "utils/mardyn_assert.h"
#include "utils/Logger.h"
#include "parallel/DomainDecompBase.h"


void CommunicationPartnerWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	Log::global_log->info() << "Write frequency: " << _writeFrequency << std::endl;

	if(_writeFrequency == 0) {
		Log::global_log->error() << "Write frequency must be a positive nonzero integer, but is " << _writeFrequency << std::endl;
		MARDYN_EXIT(-1);
	}

	std::string HaloParticleType = "unknown";
	xmlconfig.getNodeValue("type", HaloParticleType);

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	Log::global_log->info() << "Output prefix: " << _outputPrefix << std::endl;

	int incremental = 1;
	xmlconfig.getNodeValue("incremental", incremental);
	_incremental = (incremental != 0);
	Log::global_log->info() << "Incremental numbers: " << _incremental << std::endl;

	int appendTimestamp = 0;
	xmlconfig.getNodeValue("appendTimestamp", appendTimestamp);
	if(appendTimestamp > 0) {
		_appendTimestamp = true;
	}else{
		_appendTimestamp = false;
	}
	Log::global_log->info() << "Append timestamp: " << _appendTimestamp << std::endl;
}

void CommunicationPartnerWriter::init(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
                            Domain * /*domain*/) {}

void CommunicationPartnerWriter::afterForces(ParticleContainer *particleContainer, DomainDecompBase* domainDecomp,
                               unsigned long simstep) {
	if( simstep % _writeFrequency == 0 ) {
		std::stringstream filenamestream;
		filenamestream << _outputPrefix << "-rank" << domainDecomp->getRank();

		if(_incremental) {
			/* align file numbers with preceding '0's in the required range from 0 to _numberOfTimesteps. */
			unsigned long numTimesteps = _simulation.getNumTimesteps();
			int num_digits = (int) ceil( log( double( numTimesteps / _writeFrequency ) ) / log(10.) );
			filenamestream << "-" << aligned_number( simstep / _writeFrequency, num_digits, '0' );
		}
		if(_appendTimestamp) {
			filenamestream << "-" << gettimestring();
		}

		filenamestream << ".commPartners.dat";

		std::string filename = filenamestream.str();

		domainDecomp->printCommunicationPartners(filename);

	}
}

void CommunicationPartnerWriter::finish(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
							  Domain * /*domain*/) {
}
