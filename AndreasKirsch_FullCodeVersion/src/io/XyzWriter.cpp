// XyzWriter.cpp

#include "io/XyzWriter.h"

#include "Common.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"

#include <fstream>
#include <sstream>

using namespace std;

XyzWriter::XyzWriter(unsigned long writeFrequency, string filename, unsigned long numberOfTimesteps, bool incremental) {
	_filename = filename;
	_writeFrequency = writeFrequency;
	_incremental = incremental;
	_numberOfTimesteps = numberOfTimesteps;

	if (filename == "default")
		_filenameisdate = true;
	else
		_filenameisdate = false;
}

XyzWriter::~XyzWriter(){}

void XyzWriter::initOutput(ParticleContainer* particleContainer,
			   DomainDecompBase* domainDecomp, Domain* domain){
}

void XyzWriter::doOutput( ParticleContainer* particleContainer,
													DomainDecompBase* domainDecomp, Domain* domain,
			  unsigned long simstep, list<ChemicalPotential>* lmu ) 
{
	if( simstep % _writeFrequency == 0) {
		stringstream filenamestream;
		if(_filenameisdate) {
			filenamestream << "mardyn" << gettimestring() << ".out";
		} else {
			filenamestream << _filename;
		}

		if( _incremental ) {
			filenamestream << "-";
			/* align file numbers with preceding '0's in the required range from 0 to _numberOfTimesteps. */
			int num_digits = ceil( log( double( _numberOfTimesteps / _writeFrequency ) ) / log(10.) );
			filenamestream << aligned_number( simstep / _writeFrequency, num_digits, '0' ) << ".xyz";
		} else {
			filenamestream << ".xyz";
		}
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

void XyzWriter::finishOutput( ParticleContainer* particleContainer,
			     DomainDecompBase* domainDecomp, Domain* domain ){
}
