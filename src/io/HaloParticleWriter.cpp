#include "io/HaloParticleWriter.h"

#include <sstream>
#include <string>

#include "Common.h"
#include "Domain.h"
#include "utils/Logger.h"
#include "parallel/DomainDecompBase.h"


void HaloParticleWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << std::endl;

	if(_writeFrequency == 0) {
		global_log->error() << "Write frequency must be a positive nonzero integer, but is " << _writeFrequency << std::endl;
		Simulation::exit(-1);
	}
	
	std::string HaloParticleType = "unknown";
	xmlconfig.getNodeValue("type", HaloParticleType);

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << std::endl;

	_incremental = true;
	xmlconfig.getNodeValue("incremental", _incremental);

	global_log->info() << "Incremental numbers: " << _incremental << std::endl;

	_appendTimestamp = false;
	xmlconfig.getNodeValue("appendTimestamp", _appendTimestamp);

	global_log->info() << "Append timestamp: " << _appendTimestamp << std::endl;
}

void HaloParticleWriter::init(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
                            Domain * /*domain*/) {}

void HaloParticleWriter::afterForces(ParticleContainer *particleContainer, DomainDecompBase* domainDecomp,
                               unsigned long simstep) {
	if( simstep % _writeFrequency == 0 ) {
		std::stringstream filenamestream;
		filenamestream << _outputPrefix << "-rank" << domainDecomp->getRank();

		if(_incremental) {
			/* align file numbers with preceding '0's in the required range from 0 to _numberOfTimesteps. */
			unsigned long numTimesteps = global_simulation->getNumTimesteps();
			int num_digits = (int) ceil( log( double( numTimesteps / _writeFrequency ) ) / log(10.) );
			filenamestream << "-" << aligned_number( simstep / _writeFrequency, num_digits, '0' );
		}
		if(_appendTimestamp) {
			filenamestream << "-" << gettimestring();
		}

		filenamestream << ".halos.dat";

		std::string filename = filenamestream.str();

		//global_simulation->getDomain()->writeCheckpoint(filename, particleContainer, domainDecomp, _simulation.getSimulationTime(), false);
		double rmin[3],rmax[3];
		domainDecomp->getBoundingBoxMinMax(global_simulation->getDomain(), rmin,
				rmax);
		std::ofstream checkpointfilestream;
		checkpointfilestream.open(filename.c_str());
		checkpointfilestream.precision(20);

		for (auto tempMolecule = particleContainer->iterator(ParticleIterator::ALL_CELLS);
				tempMolecule.isValid(); ++tempMolecule) {
			double r[3];
			for (int i = 0; i < 3; i++) {
				r[i] = tempMolecule->r(i);
			}
			if ((r[0] < rmin[0] or r[0] >= rmax[0])
					or (r[1] < rmin[1] or r[1] >= rmax[1])
					or (r[2] < rmin[2] or r[2] >= rmax[2])) {
				checkpointfilestream <<"cell "<< tempMolecule.getCellIndex() << ": ";
				tempMolecule->write(checkpointfilestream);
			}
		}
		checkpointfilestream.close();

	}
}

void HaloParticleWriter::finish(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
							  Domain * /*domain*/) {
}
