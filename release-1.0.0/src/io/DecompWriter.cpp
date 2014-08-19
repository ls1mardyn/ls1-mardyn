#include "io/DecompWriter.h"

#include <fstream>
#include <sstream>

#include "Common.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "utils/Logger.h"

using Log::global_log;
using namespace std;

DecompWriter::DecompWriter(unsigned long writeFrequency, string mode, string outputPrefix, bool incremental) {
	_outputPrefix = outputPrefix;
	_mode = mode;
	_writeFrequency = writeFrequency;
	_incremental = incremental;

	if (outputPrefix== "default") {
		_appendTimestamp = true;
	}
	else {
		_appendTimestamp = false;
	}
}

DecompWriter::~DecompWriter(){}

void DecompWriter::readXML(XMLfileUnits& xmlconfig) {
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
	if(appendTimestamp > 0) {
		_appendTimestamp = true;
	}
	global_log->info() << "Append timestamp: " << _appendTimestamp << endl;

	xmlconfig.getNodeValue("mode", _mode);
	global_log->info() << "Mode: " << _mode << endl;
}


void DecompWriter::initOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {}

void DecompWriter::doOutput( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep, list<ChemicalPotential>* lmu ) {
	if(simstep % _writeFrequency == 0) {
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
		filenamestream << ".decomp";

		domainDecomp->printDecomp(filenamestream.str(), domain);
		
		if(_mode=="withParticles"){
			int ownRank = domainDecomp->getRank();
			for(int process = 0; process < domainDecomp->getNumProcs(); process++){
				if(ownRank==process){
					ofstream decompstrm(filenamestream.str().c_str(), ios::app);
					if(ownRank==0) {
						decompstrm << "particleData xyz" << endl;
					}
					Molecule* moleculePtr;
					for(moleculePtr = particleContainer->begin(); moleculePtr != particleContainer->end(); moleculePtr = particleContainer->next()) {
						decompstrm << moleculePtr->r(0) << "\t" << moleculePtr->r(1) << "\t" << moleculePtr->r(2) << endl;
					}
					decompstrm.close();
				}
				domainDecomp->barrier();
			}
		}
		else if(domainDecomp->getRank()==0){
			ofstream decompstrm(filenamestream.str().c_str(), ios::app);
			decompstrm << "particleData none" << endl;
			decompstrm.close();
		}
	}  
}

void DecompWriter::finishOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {}
