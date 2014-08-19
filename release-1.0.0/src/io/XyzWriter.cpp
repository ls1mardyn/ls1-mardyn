#include "io/XyzWriter.h"

#include <fstream>
#include <sstream>

#include "Common.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "Simulation.h"

using Log::global_log;
using namespace std;

XyzWriter::XyzWriter(unsigned long writeFrequency, string outputPrefix, bool incremental) {
	_outputPrefix= outputPrefix;
	_writeFrequency = writeFrequency;
	_incremental = incremental;

	if (outputPrefix == "default") {
		_appendTimestamp = true;
	}
	else {
		_appendTimestamp = false;
	}
}

XyzWriter::~XyzWriter(){}

void XyzWriter::readXML(XMLfileUnits& xmlconfig) {
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
}

void XyzWriter::initOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {}

void XyzWriter::doOutput( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep, list<ChemicalPotential>* lmu ) {
	if( simstep % _writeFrequency == 0) {
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
		filenamestream << ".xyz";
		
		int ownRank = domainDecomp->getRank();
		if( ownRank == 0 ) {
			ofstream xyzfilestream( filenamestream.str(). c_str() );
			xyzfilestream << domain->getglobalNumMolecules() << endl;
			xyzfilestream << "comment line" << endl;
			xyzfilestream.close();
		}
		for( int process = 0; process < domainDecomp->getNumProcs(); process++ ){
			domainDecomp->barrier();
			if( ownRank == process ){
				ofstream xyzfilestream( filenamestream.str().c_str(), ios::app );
				Molecule* tempMol;
				for( tempMol = particleContainer->begin(); tempMol != particleContainer->end(); tempMol = particleContainer->next()){
					if( tempMol->componentid() == 0) { xyzfilestream << "Ar ";}
					else if( tempMol->componentid() == 1 ) { xyzfilestream << "Xe ";}
					else if( tempMol->componentid() == 2 ) { xyzfilestream << "C ";}
					else if( tempMol->componentid() == 3 ) { xyzfilestream << "O ";}
					else { xyzfilestream << "H ";}
					xyzfilestream << tempMol->r(0) << "\t" << tempMol->r(1) << "\t" << tempMol->r(2) << endl;
				}
				xyzfilestream.close();
			}
		}
	}
}

void XyzWriter::finishOutput( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain ) {}
