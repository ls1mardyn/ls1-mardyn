#include "KDTreePrinter.h"

#include "Common.h"
#include "Simulation.h"
#include "utils/mardyn_assert.h"
#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"
#include "parallel/DomainDecompBase.h"
#ifdef ENABLE_MPI
#include "parallel/KDDecomposition.h"
#endif

void KDTreePrinter::init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {}

void KDTreePrinter::readXML(XMLfileUnits &xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	Log::global_log->info() << "Write frequency: " << _writeFrequency << std::endl;

	if (_writeFrequency == 0) {
		std::ostringstream error_message;error_message << "Write frequency must be a positive nonzero integer, but is " << _writeFrequency << std::endl;
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
	if (appendTimestamp > 0) {
		_appendTimestamp = true;
	} else {
		_appendTimestamp = false;
	}
	Log::global_log->info() << "Append timestamp: " << _appendTimestamp << std::endl;
}

void KDTreePrinter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
							unsigned long simstep) {
	if (simstep % _writeFrequency == 0 and domainDecomp->getRank() == 0) {
#ifdef ENABLE_MPI
		auto kdd = dynamic_cast<KDDecomposition*>(domainDecomp);
		if(kdd == nullptr){
			Log::global_log->warning() << "KDTreePrinter cannot print KDD, as the decomposition is not kdd." << std::endl;
			return;
		}
		std::stringstream filenamestream;
		filenamestream << _outputPrefix;

		if (_incremental) {
			/* align file numbers with preceding '0's in the required range from 0 to _numberOfTimesteps. */
			unsigned long numTimesteps = _simulation.getNumTimesteps();
			int num_digits = (int)ceil(log(double(numTimesteps / _writeFrequency)) / log(10.));
			filenamestream << "-" << aligned_number(simstep / _writeFrequency, num_digits, '0');
		}
		if (_appendTimestamp) {
			filenamestream << "-" << gettimestring();
		}

		filenamestream << ".dat";

		std::string filename = filenamestream.str();
		std::ofstream filestream(filename.c_str(), std::ios_base::app);
		filestream << "step " << simstep << ":" << std::endl;
		kdd->printTree(filestream);
#else
		Log::global_log->warning() << "KDTreePrinter cannot print KDD, as MPI is disabled." << std::endl;
#endif
	}
}

void KDTreePrinter::finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {}

std::string KDTreePrinter::getPluginName() { return "KDTreePrinter"; }
