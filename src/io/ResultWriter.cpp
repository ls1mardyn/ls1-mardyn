#include "io/ResultWriter.h"

#include <chrono>
#include <fstream>

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "utils/mardyn_assert.h"


void ResultWriter::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	Log::global_log->info() << "[ResultWriter] Write frequency: " << _writeFrequency << std::endl;
	if (_writeFrequency <= 0) {
		Log::global_log->error() << "[ResultWriter] Write frequency must be a positive nonzero integer, but is " << _writeFrequency << std::endl;
		mardyn_exit(123);
	}

	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	Log::global_log->info() << "[ResultWriter] Output prefix: " << _outputPrefix << std::endl;

	size_t acc_steps = 1000;
	xmlconfig.getNodeValue("accumulation_steps", acc_steps);
	_U_pot_acc = std::make_unique<Accumulator<double>>(acc_steps);
	_U_kin_acc = std::make_unique<Accumulator<double>>(acc_steps);
	_p_acc = std::make_unique<Accumulator<double>>(acc_steps);
	Log::global_log->info() << "[ResultWriter] Accumulation steps: " << acc_steps << std::endl;

	xmlconfig.getNodeValue("writeprecision", _writePrecision);
	Log::global_log->info() << "[ResultWriter] Write precision: " << _writePrecision << std::endl;
	_writeWidth = _writePrecision + 15;  // Adding a width of 15 to have enough whitespace between columns
}

void ResultWriter::init(ParticleContainer * /*particleContainer*/,
                        DomainDecompBase *domainDecomp, Domain * /*domain*/) {

	if(domainDecomp->getRank() == 0) {
		const std::string resultfile(_outputPrefix+".res");
		std::ofstream resultStream;
		resultStream.open(resultfile.c_str(), std::ios::out);
		const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
		tm unused{};
		const auto nowStr = std::put_time(localtime_r(&now, &unused), "%c");
		resultStream << "# ls1 MarDyn simulation started at " << nowStr << std::endl;
		resultStream << "# Averages are the accumulated values over " << _U_pot_acc->getWindowLength()  << " time steps."<< std::endl;
		resultStream << std::setw(10) << "# step" << std::setw(_writeWidth) << "time"
			<< std::setw(_writeWidth) << "U_pot"
			<< std::setw(_writeWidth) << "U_pot_avg"
			<< std::setw(_writeWidth) << "U_kin"
			<< std::setw(_writeWidth) << "U_kin_avg"
			<< std::setw(_writeWidth) << "p"
			<< std::setw(_writeWidth) << "p_avg"
			<< std::setw(_writeWidth) << "beta_trans"
			<< std::setw(_writeWidth) << "beta_rot"
			<< std::setw(_writeWidth) << "c_v"
			<< std::setw(_writeWidth) << "N"
			<< std::endl;
		resultStream.close();
	}
}

void ResultWriter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                           unsigned long simstep) {

	// Writing of cavities now handled by CavityWriter

	unsigned long globalNumMolecules = domain->getglobalNumMolecules(true, particleContainer, domainDecomp);
	double cv = domain->cv();
	double ekin = domain->getGlobalUkinTrans()+domain->getGlobalUkinRot();

	_U_pot_acc->addEntry(domain->getGlobalUpot());
	_U_kin_acc->addEntry(ekin);
	_p_acc->addEntry(domain->getGlobalPressure());
	if ((domainDecomp->getRank() == 0) && (simstep % _writeFrequency == 0)){
		const std::string resultfile(_outputPrefix+".res");
		std::ofstream resultStream;
		resultStream.open(resultfile.c_str(), std::ios::app);
		auto printOutput = [&](auto value) {
			resultStream << std::setw(_writeWidth) << std::scientific << std::setprecision(_writePrecision) << value;
		};
		resultStream << std::setw(10) << simstep;
		printOutput(_simulation.getSimulationTime());
		printOutput(domain->getGlobalUpot());
		printOutput(_U_pot_acc->getAverage());
		printOutput(ekin);
		printOutput(_U_kin_acc->getAverage());
		printOutput(domain->getGlobalPressure());
		printOutput(_p_acc->getAverage());
		printOutput(domain->getGlobalBetaTrans());
		printOutput(domain->getGlobalBetaRot());
		printOutput(cv);
		printOutput(globalNumMolecules);
		resultStream << std::endl;
		resultStream.close();
	}
}

void ResultWriter::finish(ParticleContainer * /*particleContainer*/,
						  DomainDecompBase *domainDecomp, Domain * /*domain*/){

	if (domainDecomp->getRank() == 0) {
		const std::string resultfile(_outputPrefix+".res");
		std::ofstream resultStream;
		resultStream.open(resultfile.c_str(), std::ios::app);
		const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
		tm unused{};
		const auto nowStr = std::put_time(localtime_r(&now, &unused), "%c");
		resultStream << "# ls1 mardyn simulation finished at " << nowStr << std::endl;
		resultStream.close();
	}
}
