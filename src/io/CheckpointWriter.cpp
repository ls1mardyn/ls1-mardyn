// CheckpointWriter.cpp

#include "io/CheckpointWriter.h"

#include "Common.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"

#include <sstream>

using namespace std;

CheckpointWriter::CheckpointWriter(unsigned long writeFrequency, string filename, unsigned long numberOfTimesteps, bool incremental) {
	_filename = filename;
	_writeFrequency = writeFrequency;
	_incremental = incremental;
	_numberOfTimesteps = numberOfTimesteps;

	if (filename == "default")
		_filenameisdate = true;
	else
		_filenameisdate = false;
}

CheckpointWriter::~CheckpointWriter(){}

void CheckpointWriter::initOutput(ParticleContainer* particleContainer,
					DomainDecompBase* domainDecomp, Domain* domain) {
}

void CheckpointWriter::doOutput( ParticleContainer* particleContainer,
        DomainDecompBase* domainDecomp, 
        Domain* domain,
        unsigned long simstep, 
        list<ChemicalPotential>* lmu )
{
	if( simstep % _writeFrequency != 0 ) return;
	
		stringstream filenamestream;
		if(_filenameisdate) {
			filenamestream << "mardyn" << gettimestring();
		} else {
			filenamestream << _filename;
		}

		if(_incremental) {
			/* align file numbers with preceding '0's in the required range from 0 to _numberOfTimesteps. */
			int num_digits = ceil( log( double( _numberOfTimesteps / _writeFrequency ) ) / log(10.) );
			filenamestream << "-" << aligned_number( simstep / _writeFrequency, num_digits, '0' );
		}
		filenamestream << ".restart.xdr";

	string filename = filenamestream.str();
	domain->writeCheckpoint(filename, particleContainer, domainDecomp);
}

void CheckpointWriter::finishOutput(
	 ParticleContainer* particleContainer,
	 DomainDecompBase* domainDecomp, Domain* domain
) {
}
