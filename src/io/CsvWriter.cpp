#include "io/CsvWriter.h"

#include <fstream>
#include <sstream>

#include "Common.h"
#include "Domain.h"
#include "Simulation.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"

using Log::global_log;
using namespace std;

void CsvWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << std::endl;

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << std::endl;

	xmlconfig.getNodeValue("incremental", _incremental);
	global_log->info() << "Incremental numbers: " << _incremental << std::endl;

	xmlconfig.getNodeValue("appendTimestamp", _appendTimestamp);
	global_log->info() << "Append timestamp: " << _appendTimestamp << std::endl;
}

void CsvWriter::init(ParticleContainer* /*particleContainer*/, DomainDecompBase* /*domainDecomp*/, Domain* /*domain*/) {
}

void CsvWriter::endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* /*domain*/,
						unsigned long simstep) {
	if (simstep % _writeFrequency == 0) {
		vector<Component>* components = _simulation.getEnsemble()->getComponents();
		stringstream filenamestream;
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
		filenamestream << ".csv";

		int ownRank = domainDecomp->getRank();
		if (ownRank == 0) {
			ofstream csvfilestream(filenamestream.str().c_str());
			csvfilestream << "atom,x,y,z"
						  << "\n";  
			csvfilestream.close();
		}
		for (int process = 0; process < domainDecomp->getNumProcs(); process++) {
			domainDecomp->barrier();
			if (ownRank == process) {
				ofstream csvfilestream(filenamestream.str().c_str(), ios::app);
                // FIXME: This is just copied from the XyzWriter and it is wrong! (Look at CO2 example. Sites don't (not always at least) translate into atoms)
                /*
                * TODO:
                *   - Do NOT infer atom type, as this code from XyzWriter seems highly dubious
                *   - Output only atom positions? (Use extra configs for representation like MmpldWriter?)
                *   - Output additional information such as forces in seperate columns?
                *   - Maybe make code a bit nicer and more modular
                *   - Provide documentation for visualization in ParaView (open CSVs as Group -> TableToPoints -> Glyph)
                */
				for (ParticleIterator tempMol = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
					 tempMol.isValid(); ++tempMol) {
					for (unsigned i = 0; i < tempMol->numLJcenters(); i++) {
						if (tempMol->componentid() == 0) {
							csvfilestream << "Ar";
						} else if (tempMol->componentid() == 1) {
							csvfilestream << "Xe";
						} else if (tempMol->componentid() == 2) {
							csvfilestream << "C";
						} else if (tempMol->componentid() == 3) {
							csvfilestream << "O";
						} else {
							csvfilestream << "H";
						}
						csvfilestream << _separator << tempMol->r(0) + tempMol->ljcenter_d(i)[0] << _separator
									  << tempMol->r(1) + tempMol->ljcenter_d(i)[1] << _separator
									  << tempMol->r(2) + tempMol->ljcenter_d(i)[2] << "\n";
					}
					for (unsigned i = 0; i < tempMol->numDipoles(); i++) {
						if (tempMol->componentid() == 0) {
							csvfilestream << "O";
						} else if (tempMol->componentid() == 1) {
							csvfilestream << "H";
						} else if (tempMol->componentid() == 2) {
							csvfilestream << "Xe";
						} else if (tempMol->componentid() == 3) {
							csvfilestream << "Ar";
						} else {
							csvfilestream << "C";
						}
						csvfilestream << _separator << tempMol->r(0) + tempMol->dipole_d(i)[0] << _separator
									  << tempMol->r(1) + tempMol->dipole_d(i)[1] << _separator
									  << tempMol->r(2) + tempMol->dipole_d(i)[2] << "\n";
					}
					for (unsigned i = 0; i < tempMol->numQuadrupoles(); i++) {
						if (tempMol->componentid() == 0) {
							csvfilestream << "C";
						} else if (tempMol->componentid() == 1) {
							csvfilestream << "Ar";
						} else if (tempMol->componentid() == 2) {
							csvfilestream << "H";
						} else if (tempMol->componentid() == 3) {
							csvfilestream << "Xe";
						} else {
							csvfilestream << "O";
						}
						csvfilestream << _separator << tempMol->r(0) + tempMol->quadrupole_d(i)[0] << _separator
									  << tempMol->r(1) + tempMol->quadrupole_d(i)[1] << _separator
									  << tempMol->r(2) + tempMol->quadrupole_d(i)[2] << "\n";
					}
					for (unsigned i = 0; i < tempMol->numCharges(); i++) {
						if (tempMol->componentid() == 0) {
							csvfilestream << "H";
						} else if (tempMol->componentid() == 1) {
							csvfilestream << "O";
						} else if (tempMol->componentid() == 2) {
							csvfilestream << "Xe";
						} else if (tempMol->componentid() == 3) {
							csvfilestream << "Ar";
						} else {
							csvfilestream << "C";
						}
						csvfilestream << _separator << tempMol->r(0) + tempMol->charge_d(i)[0] << _separator
									  << tempMol->r(1) + tempMol->charge_d(i)[1] << _separator
									  << tempMol->r(2) + tempMol->charge_d(i)[2] << "\n";
					}
				}
				csvfilestream.close();
			}
		}
	}
}

void CsvWriter::finish(ParticleContainer* /*particleContainer*/, DomainDecompBase* /*domainDecomp*/,
					   Domain* /*domain*/) {}
