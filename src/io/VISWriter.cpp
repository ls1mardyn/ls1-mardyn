// VISWriter.cpp

#include "io/VISWriter.h"
#include "Common.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"

#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

VISWriter::VISWriter(unsigned long writeFrequency, string filename, unsigned long numberOfTimesteps, bool incremental) {
	_filename = filename;
	_writeFrequency = writeFrequency;
	_incremental = incremental;
	_numberOfTimesteps = numberOfTimesteps;

	if (filename == "default")
		_filenameisdate = true;
	else
		_filenameisdate = false;
}

VISWriter::~VISWriter(){}

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
		if (_filenameisdate) {
			filenamestream << "mardyn" << gettimestring();
		}
		else {
			filenamestream << _filename;
		}

		if (_incremental) {
			/* align file numbers with preceding '0's in the required range from 0 to _numberOfTimesteps. */
			int num_digits = ceil(log(double(_numberOfTimesteps / _writeFrequency)) / log(10.));
			filenamestream << aligned_number(simstep / _writeFrequency, num_digits, '0');
		}
		filenamestream << ".vis";

		ofstream visittfstrm(filenamestream.str().c_str());

		if ((domain->ownrank() == 0) && (!_wroteVIS)) {
			visittfstrm << "      id t          x          y          z     q0     q1     q2     q3        c\n";
			this->_wroteVIS = true;
		}
		else if (domain->ownrank() == 0) {
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

void VISWriter::finishOutput(ParticleContainer* particleContainer,
                             DomainDecompBase* domainDecomp, Domain* domain) {
}
