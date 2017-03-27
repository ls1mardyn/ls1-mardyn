#include "io/BinaryCheckpointWriter.h"

#include <sstream>
#include <string>

#include "Common.h"
#include "Domain.h"
#include "utils/Logger.h"

using Log::global_log;
using namespace std;

BinaryCheckpointWriter::BinaryCheckpointWriter(unsigned long writeFrequency,
		string outputPrefix, bool incremental) {
	_outputPrefix = outputPrefix;
	_writeFrequency = writeFrequency;
	_incremental = incremental;

	if (outputPrefix == "default") {
		_appendTimestamp = true;
	} else {
		_appendTimestamp = false;
	}
}

BinaryCheckpointWriter::~BinaryCheckpointWriter() {
}

void BinaryCheckpointWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << endl;

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << endl;

	int incremental = 1;
	xmlconfig.getNodeValue("incremental", incremental);
	_incremental = (incremental != 0);
	global_log->info() << "Incremental numbers: " << _incremental << endl;

	int appendTimestamp = 0;
	xmlconfig.getNodeValue("appendTimestamp", appendTimestamp);
	if (appendTimestamp > 0) {
		_appendTimestamp = true;
	}
	global_log->info() << "Append timestamp: " << _appendTimestamp << endl;
}

void BinaryCheckpointWriter::initOutput(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain) {
}

void BinaryCheckpointWriter::doOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
		unsigned long simstep, list<ChemicalPotential>* /*lmu*/, map<unsigned, CavityEnsemble>* /*mcav*/) {
	if (simstep % _writeFrequency == 0) {
		stringstream filenamestream;
		filenamestream << _outputPrefix;

		if (_incremental) {
			/* align file numbers with preceding '0's in the required range from 0 to _numberOfTimesteps. */

			unsigned long numTimesteps = _simulation.getNumTimesteps();
			int num_digits = (int) ceil(
					log(double(numTimesteps / _writeFrequency)) / log(10.));
			filenamestream << "-" << aligned_number(simstep / _writeFrequency,
					num_digits, '0');
		}
		if (_appendTimestamp) {
			filenamestream << "-" << gettimestring();
		}
		filenamestream << ".restart";

		string filename = filenamestream.str();
		domain->writeCheckpoint(filename, particleContainer, domainDecomp, _simulation.getSimulationTime(), true);
	}
}

void BinaryCheckpointWriter::finishOutput(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain) {
}
