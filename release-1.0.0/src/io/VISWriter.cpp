#include "io/VISWriter.h"

#include <iomanip>
#include <fstream>
#include <sstream>

#include "Common.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "utils/Logger.h"

using Log::global_log;
using namespace std;

VISWriter::VISWriter(unsigned long writeFrequency, string outputPrefix, bool incremental) {
	_outputPrefix = outputPrefix;
	_writeFrequency = writeFrequency;
	_incremental = incremental;

	if (outputPrefix == "default") {
		_appendTimestamp = true;
	}
	else {
		_appendTimestamp = false;
	}
}

VISWriter::~VISWriter(){}

void VISWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << endl;

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << endl;
	
	int incremental = 1;
	xmlconfig.getNodeValue("incremental", incremental);
	global_log->info() << "Incremental numbers: " << _incremental << endl;
	
	int appendTimestamp = 0;
	xmlconfig.getNodeValue("appendTimestamp", appendTimestamp);
	if(appendTimestamp > 0) {
		_appendTimestamp = true;
	}
	global_log->info() << "Append timestamp: " << _appendTimestamp << endl;
}

void VISWriter::initOutput(ParticleContainer* particleContainer,
                           DomainDecompBase* domainDecomp, Domain* domain)
{
	this->_wroteVIS = false;
}

void VISWriter::doOutput(ParticleContainer* particleContainer,
                         DomainDecompBase* domainDecomp, Domain* domain,
                         unsigned long simstep, list<ChemicalPotential>* lmu) {
	if (simstep % _writeFrequency == 0) {
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
		filenamestream << ".vis";
		
		ofstream visittfstrm(filenamestream.str().c_str());

		int rank = domainDecomp->getRank();
		if ((rank== 0) && (!_wroteVIS)) {
			visittfstrm << "      id t          x          y          z     q0     q1     q2     q3        c\n";
			_wroteVIS = true;
		}
		else if (rank == 0) {
            visittfstrm << "#\n";
        }

		// originally VIS files had a fixed width of 8 (and no t), here I use 12 (with 2 for t)
		//ostrm << "t           x           y           z          q0          q1          q2          q3" << endl;
		for (Molecule* pos = particleContainer->begin(); pos != particleContainer->end(); pos = particleContainer->next()) {
			bool halo = false;
			for (unsigned short d = 0; d < 3; d++) {
				if ((pos->r(d) < particleContainer->getBoundingBoxMin(d)) || (pos->r(d) > particleContainer->getBoundingBoxMax(d))) {
					halo = true;
					break;
				}
			}
			if (!halo) {
				visittfstrm << setiosflags(ios::fixed) << setw(8) << pos->id() << setw(2)
				            << pos->componentid() << setprecision(3);
				for (unsigned short d = 0; d < 3; d++) visittfstrm << setw(11) << pos->r(d);
				visittfstrm << setprecision(3) << setw(7) << pos->q().qw() << setw(7) << pos->q().qx()
				            << setw(7) << pos->q().qy()<< setw(7) << pos->q().qz()
				            << setw(9) << right << 0 << "\n";
			}
		}
		visittfstrm.close();
	}
}

void VISWriter::finishOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {}
