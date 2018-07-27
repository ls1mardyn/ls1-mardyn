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

void XyzWriter::initOutput(ParticleContainer* /*particleContainer*/, DomainDecompBase* /*domainDecomp*/, Domain* /*domain*/) {}

void XyzWriter::doOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* /*domain*/,
		unsigned long simstep, list<ChemicalPotential>* /*lmu*/, map<unsigned, CavityEnsemble>* /*mcav*/) {
	if( simstep % _writeFrequency == 0) {
		vector<Component>*  components = _simulation.getEnsemble()->getComponents();
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
			unsigned number = 0;
			for (unsigned i=0; i< components->size(); i++){
				number += (*components)[i].getNumMolecules()*((*components)[i].numLJcenters() + (*components)[i].numDipoles() + (*components)[i].numCharges() + (*components)[i].numQuadrupoles());
			}
			xyzfilestream << number << endl;
			xyzfilestream << "comment line" << endl;
			xyzfilestream.close();
		}
		for( int process = 0; process < domainDecomp->getNumProcs(); process++ ){
			domainDecomp->barrier();
			if( ownRank == process ){
				ofstream xyzfilestream( filenamestream.str().c_str(), ios::app );
				for(ParticleIterator tempMol = particleContainer->iteratorBegin(); tempMol != particleContainer->iteratorEnd(); ++tempMol){
					for (unsigned i=0; i< tempMol->numLJcenters(); i++){
						if( tempMol->componentid() == 0) { xyzfilestream << "Ar ";}
						else if( tempMol->componentid() == 1 ) { xyzfilestream << "Xe ";}
						else if( tempMol->componentid() == 2 ) { xyzfilestream << "C ";}
						else if( tempMol->componentid() == 3 ) { xyzfilestream << "O ";}
						else { xyzfilestream << "H ";}
						xyzfilestream << tempMol->r(0) + tempMol->ljcenter_d(i)[0] << "\t" << tempMol->r(1) + tempMol->ljcenter_d(i)[1] << "\t" << tempMol->r(2) + tempMol->ljcenter_d(i)[2] << endl;
					}
					for (unsigned i=0; i< tempMol->numDipoles(); i++){
						if( tempMol->componentid() == 0) { xyzfilestream << "O ";}
						else if( tempMol->componentid() == 1 ) { xyzfilestream << "H ";}
						else if( tempMol->componentid() == 2 ) { xyzfilestream << "Xe ";}
						else if( tempMol->componentid() == 3 ) { xyzfilestream << "Ar ";}
						else { xyzfilestream << "C ";}
						xyzfilestream << tempMol->r(0) + tempMol->dipole_d(i)[0] << "\t" << tempMol->r(1) + tempMol->dipole_d(i)[1] << "\t" << tempMol->r(2) + tempMol->dipole_d(i)[2] << endl;
					}
					for (unsigned i=0; i< tempMol->numQuadrupoles(); i++){
						if( tempMol->componentid() == 0) { xyzfilestream << "C ";}
						else if( tempMol->componentid() == 1 ) { xyzfilestream << "Ar ";}
						else if( tempMol->componentid() == 2 ) { xyzfilestream << "H ";}
						else if( tempMol->componentid() == 3 ) { xyzfilestream << "Xe ";}
						else { xyzfilestream << "O ";}
						xyzfilestream << tempMol->r(0) + tempMol->quadrupole_d(i)[0] << "\t" << tempMol->r(1) + tempMol->quadrupole_d(i)[1] << "\t" << tempMol->r(2) + tempMol->quadrupole_d(i)[2] << endl;
					}
					for (unsigned i=0; i< tempMol->numCharges(); i++){
						if( tempMol->componentid() == 0) { xyzfilestream << "H ";}
						else if( tempMol->componentid() == 1 ) { xyzfilestream << "O ";}
						else if( tempMol->componentid() == 2 ) { xyzfilestream << "Xe ";}
						else if( tempMol->componentid() == 3 ) { xyzfilestream << "Ar ";}
						else { xyzfilestream << "C ";}
						xyzfilestream << tempMol->r(0) + tempMol->charge_d(i)[0] << "\t" << tempMol->r(1) + tempMol->charge_d(i)[1] << "\t" << tempMol->r(2) + tempMol->charge_d(i)[2] << endl;
					}
				}
				xyzfilestream.close();
			}
		}
	}
}

void XyzWriter::finishOutput( ParticleContainer* /*particleContainer*/, DomainDecompBase* /*domainDecomp*/, Domain* /*domain*/ ) {}
