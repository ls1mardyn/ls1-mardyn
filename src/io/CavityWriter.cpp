#include "io/CavityWriter.h"

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

CavityWriter::CavityWriter(unsigned long writeFrequency, string outputPrefix, bool incremental) {
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

CavityWriter::~CavityWriter(){}

void CavityWriter::readXML(XMLfileUnits& xmlconfig) {
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

void CavityWriter::init(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
						Domain * /*domain*/) {}

void CavityWriter::endStep(ParticleContainer * /*particleContainer*/, DomainDecompBase *domainDecomp,
                           Domain * /*domain*/, unsigned long simstep) {

	map<unsigned, CavityEnsemble> * mcav = global_simulation->getMcav();

	if( simstep % _writeFrequency == 0) {
                map<unsigned, CavityEnsemble>::iterator ceit;
           
                map<unsigned, stringstream*> cav_filenamestream;
                for(ceit = mcav->begin(); ceit != mcav->end(); ceit++)
                {
                   cav_filenamestream[ceit->first] = new stringstream;
                   *cav_filenamestream[ceit->first] << _outputPrefix << "-c" << ceit->first;
                }

		if(_incremental) {
			unsigned long numTimesteps = _simulation.getNumTimesteps();
			int num_digits = (int) ceil( log( double( numTimesteps / _writeFrequency ) ) / log(10.) );
                        for(ceit = mcav->begin(); ceit != mcav->end(); ceit++)
                           *cav_filenamestream[ceit->first] << "-" << aligned_number( simstep / _writeFrequency, num_digits, '0' );
		}
                for(ceit = mcav->begin(); ceit != mcav->end(); ceit++)
                   *cav_filenamestream[ceit->first] << ".cav.xyz";
		
		int ownRank = domainDecomp->getRank();
		if( ownRank == 0 ) {
                        for(ceit = mcav->begin(); ceit != mcav->end(); ceit++)
                        {
                           ofstream cavfilestream( cav_filenamestream[ceit->first]->str().c_str() );
                           cavfilestream << ceit->second.numCavities() << endl;
                           cavfilestream << "comment line" << endl;
                           cavfilestream.close();
                        }
		}
		for( int process = 0; process < domainDecomp->getNumProcs(); process++ ){
			domainDecomp->barrier();
			if( ownRank == process ){
                                for(ceit = mcav->begin(); ceit != mcav->end(); ceit++)
                                {
                                   ofstream cavfilestream( cav_filenamestream[ceit->first]->str().c_str(), ios::app );
                                   
                                   map<unsigned long, Molecule*> tcav = ceit->second.activeParticleContainer();
                                   map<unsigned long, Molecule*>::iterator tcit;
                                   for(tcit = tcav.begin(); tcit != tcav.end(); tcit++)
                                   {
                                      if( ceit->first == 0 ) { cavfilestream << "C ";}
                                      else if( ceit->first == 1 ) { cavfilestream << "N ";}
                                      else if( ceit->first == 2 ) { cavfilestream << "O ";}
                                      else if( ceit->first == 3 ) { cavfilestream << "F ";}
                                      else { cavfilestream << "Ne "; }
                                      cavfilestream << tcit->second->r(0) << "\t" << tcit->second->r(1) << "\t" << tcit->second->r(2) << "\n";
                                   }
                                   
                                   cavfilestream.close();
                                }
			}
		}
	}
}

void CavityWriter::finish(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
						  Domain * /*domain*/ ) {}
