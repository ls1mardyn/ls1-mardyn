#include "io/ResultWriter.h"

#include "Domain.h"
#include "Simulation.h"
#include "parallel/DomainDecompBase.h"
#include "utils/Logger.h"

#include <fstream>


void ResultWriter::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	Log::global_log->info() << "[ResultWriter] Write frequency: " << _writeFrequency << std::endl;
	if (_writeFrequency == 0) {
		Log::global_log->error() << "[ResultWriter] Write frequency must be a positive nonzero integer, but is " << _writeFrequency << std::endl;
		Simulation::exit(-1);
	}

	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	Log::global_log->info() << "[ResultWriter] Output prefix: " << _outputPrefix << std::endl;

	xmlconfig.getNodeValue("writeprecision", _writePrecision);
	Log::global_log->info() << "[ResultWriter] Write precision: " << _writePrecision << std::endl;
	_writeWidth = _writePrecision + 15; // Adding a width of 15 to have enough whitespace between chars
}

void ResultWriter::init(ParticleContainer * /*particleContainer*/,
                        DomainDecompBase *domainDecomp, Domain * /*domain*/) {

	// Only main rank writes data to file
	if(domainDecomp->getRank() == 0) {
		const std::string resultfile(_outputPrefix+".res");
		std::ofstream resultStream;
		resultStream.open(resultfile.c_str(), std::ios::out);
		resultStream << std::setw(10) << "timestep"
			<< std::setw(_writeWidth) << "time"
			<< std::setw(_writeWidth) << "U_pot_avg"
			<< std::setw(_writeWidth) << "U_kin_avg"
			<< std::setw(_writeWidth) << "U_kinTrans_avg"
			<< std::setw(_writeWidth) << "U_kinRot_avg"
			<< std::setw(_writeWidth) << "p_avg"
			<< std::setw(_writeWidth) << "c_v"
			<< std::setw(_writeWidth) << "N"
			<< std::endl;
		resultStream.close();
	}
}

void ResultWriter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                           unsigned long simstep) {

	const uint64_t globalNumMolecules = domain->getglobalNumMolecules(true, particleContainer, domainDecomp);

	_uPot_acc += domain->getGlobalUpot();
	_p_acc += domain->getGlobalPressure();
	_uKinTrans_acc += domain->getGlobalUkinTrans();
	_uKinRot_acc += domain->getGlobalUkinRot();

	_numSamples++;

	if ((simstep % _writeFrequency == 0) and (simstep > 0UL)) {
		// Only main rank writes data to file
		if (domainDecomp->getRank() == 0) {
			const std::string resultfile(_outputPrefix+".res");
			std::ofstream resultStream;
			resultStream.open(resultfile.c_str(), std::ios::app);
			auto printOutput = [&](auto value) {
				resultStream << std::setw(_writeWidth) << std::scientific << std::setprecision(_writePrecision) << value;
			};
			resultStream << std::setw(10) << simstep;
				printOutput(_simulation.getSimulationTime());
				printOutput(_uPot_acc/_numSamples);
				printOutput((_uKinTrans_acc+_uKinRot_acc)/_numSamples);
				printOutput(_uKinTrans_acc/_numSamples);
				printOutput(_uKinRot_acc/_numSamples);
				printOutput(_p_acc/_numSamples);
				printOutput(domain->cv());
				printOutput(globalNumMolecules);
				resultStream << std::endl;
			resultStream.close();
		}

		// Reset values
		_numSamples = 0UL;
		_uPot_acc = 0.0F;
		_p_acc = 0.0F;
		_uKinTrans_acc = 0.0F;
		_uKinRot_acc = 0.0F;
	}
}
