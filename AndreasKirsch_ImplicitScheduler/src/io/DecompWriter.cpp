// DecompWriter.cpp

#include "io/DecompWriter.h"

#include "Common.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"

using namespace std;

DecompWriter::DecompWriter(unsigned long writeFrequency, string mode, string filename, unsigned long numberOfTimesteps, bool incremental) {
	_filename = filename;
	_mode = mode;
	_writeFrequency = writeFrequency;
	_incremental = incremental;
	_numberOfTimesteps = numberOfTimesteps;

	if (filename == "default")
		_filenameisdate = true;
	else
		_filenameisdate = false;
}

DecompWriter::~DecompWriter(){}

void DecompWriter::initOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain){
}

void DecompWriter::doOutput( ParticleContainer* particleContainer,
														 DomainDecompBase* domainDecomp, Domain* domain,
			     unsigned long simstep, list<ChemicalPotential>* lmu ) 
{
	if(simstep%_writeFrequency == 0) {
		stringstream filenamestream;
		if(_filenameisdate) {
			filenamestream << gettimestring() << ".out";
		} else {
			filenamestream << _filename;
		}

		if(_incremental) {
			filenamestream << "-";
			/* align file numbers with preceding '0's in the required range from 0 to _numberOfTimesteps. */
			int num_digits = ceil( log( double( _numberOfTimesteps / _writeFrequency ) ) / log(10.) );
			filenamestream << aligned_number( simstep / _writeFrequency, num_digits, '0' );
		}
		filenamestream << ".decomp";

		domainDecomp->printDecomp(filenamestream.str(), domain);
		
		if(_mode=="withParticles"){
			int ownRank = domainDecomp->getRank();
			for(int process = 0; process < domainDecomp->getNumProcs(); process++){
				if(ownRank==process){
					ofstream decompstrm(filenamestream.str().c_str(), ios::app);
					if(ownRank==0) decompstrm << "particleData xyz" << endl;
					Molecule* tempMol;
					for(tempMol = particleContainer->begin(); tempMol != particleContainer->end(); tempMol = particleContainer->next()){
						decompstrm << tempMol->r(0) << "\t" << tempMol->r(1) << "\t" << tempMol->r(2) << endl;
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

void DecompWriter::finishOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain){
}
