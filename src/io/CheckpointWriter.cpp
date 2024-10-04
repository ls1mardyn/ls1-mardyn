#include "io/CheckpointWriter.h"


#include <sstream>
#include <string>
#include <cstring>

#include "Common.h"
#include "Domain.h"
#include "Simulation.h"
#include "parallel/DomainDecompBase.h"
#include "utils/Logger.h"
#include "utils/mardyn_assert.h"


void CheckpointWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	Log::global_log->info() << "Write frequency: " << _writeFrequency << std::endl;

	if(_writeFrequency == 0) {
		std::ostringstream error_message;
		error_message << "Write frequency must be a positive nonzero integer, but is " << _writeFrequency << std::endl;
		MARDYN_EXIT(error_message);
	}

	std::string checkpointType = "unknown";
	xmlconfig.getNodeValue("type", checkpointType);
	if("ASCII" == checkpointType) {
		_useBinaryFormat = false;
	}
	else if("binary" == checkpointType) {
		_useBinaryFormat = true;
	}
	else {
		std::ostringstream error_message;
		error_message << "Unknown CheckpointWriter type '" << checkpointType << "', expected: ASCII|binary." << std::endl;
		MARDYN_EXIT(error_message);
	}

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

void CheckpointWriter::init(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
                            Domain * /*domain*/) {}

void CheckpointWriter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                               unsigned long simstep) {
	if( simstep % _writeFrequency == 0 ) {
		std::stringstream filenamestream;
		filenamestream << _outputPrefix;

		if(_incremental) {
			/* align file numbers with preceding '0's in the required range from 0 to _numberOfTimesteps. */
			unsigned long numTimesteps = _simulation.getNumTimesteps();
			int num_digits = (int) ceil( log( double( numTimesteps / _writeFrequency ) ) / log(10.) );
			filenamestream << "-" << aligned_number( simstep / _writeFrequency, num_digits, '0' );
		}
		if(_appendTimestamp) {
			filenamestream << "-" << gettimestring();
		}

		if (_useBinaryFormat) {
			filenamestream << ".restart";
		} else { /* ASCII mode */
			filenamestream << ".restart.dat";
		}

		std::string filename = filenamestream.str();
		domain->writeCheckpoint(filename, particleContainer, domainDecomp, _simulation.getSimulationTime(), _useBinaryFormat);
	}
}

void CheckpointWriter::finish(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
							  Domain * /*domain*/) {
}
