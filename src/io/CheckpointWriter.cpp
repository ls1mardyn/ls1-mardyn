#include "io/CheckpointWriter.h"


#include <sstream>
#include <string>
#include <cstring>

#include "Common.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "utils/Logger.h"


using Log::global_log;
using namespace std;

void CheckpointWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << endl;

	if(_writeFrequency == 0) {
		global_log->error() << "Write frequency must be a positive nonzero integer, but is " << _writeFrequency << endl;
		Simulation::exit(-1);
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
		global_log->error() << "Unknown CheckpointWriter type '" << checkpointType << "', expected: ASCII|binary." << endl;
		Simulation::exit(-1);
	}

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << endl;

	int incremental = 1;
	xmlconfig.getNodeValue("incremental", incremental);
	_incremental = (incremental != 0);
	global_log->info() << "Incremental numbers: " << _incremental << endl;

	int appendTimestamp = 0;
	xmlconfig.getNodeValue("appendTimestamp", appendTimestamp);
	if(appendTimestamp > 0) {
		_appendTimestamp = true;
	}else{
		_appendTimestamp = false;
	}
	global_log->info() << "Append timestamp: " << _appendTimestamp << endl;
}

void CheckpointWriter::init(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
                            Domain * /*domain*/) {}

void CheckpointWriter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                               unsigned long simstep) {
#ifdef ENABLE_ADRESS
    if(particleContainer != _simulation.getMoleculeContainer()) return;
#endif
	if( simstep % _writeFrequency == 0 ) {
		stringstream filenamestream;
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

		string filename = filenamestream.str();
		domain->writeCheckpoint(filename, _simulation.getMoleculeContainers(), domainDecomp, _simulation.getSimulationTime(), _useBinaryFormat);
	}
}

void CheckpointWriter::finish(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
							  Domain * /*domain*/) {
}
