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
		std::ostringstream error_message;
		error_message << "[ResultWriter] Write frequency must be a positive nonzero integer, but is " << _writeFrequency << std::endl;
		MARDYN_EXIT(error_message.str());
	}

	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	Log::global_log->info() << "[ResultWriter] Output prefix: " << _outputPrefix << std::endl;

	xmlconfig.getNodeValue("outputformat", _outputFormat);
	if (_outputFormat == "csv" || _outputFormat == "tab") {
		Log::global_log->info() << "[ResultWriter] Output format: " << _outputFormat << std::endl;
	} else {
		std::ostringstream error_message;
		error_message << "[ResultWriter] Wrong output format specified! Use \"csv\" or \"tab\" instead of " << _outputFormat << std::endl;
		MARDYN_EXIT(error_message.str());
	}
	// Set file extension to be .csv in case of "csv" and .res in case of "tab"
	_resultfilename = _outputPrefix+(_outputFormat == "csv" ? ".csv" : ".res");

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
		const std::string resultfile(_resultfilename);
		std::ofstream resultStream;
		resultStream.open(resultfile.c_str(), std::ios::out);

		auto formatOutput = [&](auto value) {
			if (_outputFormat == "tab") {
				resultStream << std::setw(_writeWidth) << std::scientific << std::setprecision(_writePrecision) << value;
			} else {  // csv
				resultStream << "," << std::scientific << std::setprecision(_writePrecision) << value;
			}
		};

		// Write header of output file in case of "tab"
		if (_outputFormat == "tab") {
			const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
			tm unused{};
			const auto nowStr = std::put_time(localtime_r(&now, &unused), "%c");
			resultStream << "# ls1 MarDyn simulation started at " << nowStr << std::endl;
			resultStream << "# Averages are the accumulated values over " << _U_pot_acc->getWindowLength()  << " time steps."<< std::endl;
		}
		
		if (_outputFormat == "tab") {
			// Do not write simstep in scientific notation
			resultStream << std::setw(_writeWidth) << "simstep";
		} else {
			// Do not write comma first
			resultStream << "simstep";
		}
		formatOutput("time");
		formatOutput("U_pot");
		formatOutput("U_pot_avg");
		formatOutput("U_kin");
		formatOutput("U_kin_avg");
		formatOutput("p");
		formatOutput("p_avg");
		formatOutput("beta_trans");
		formatOutput("beta_rot");
		formatOutput("c_v");
		formatOutput("N");
		resultStream << std::endl;
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
		const std::string resultfile(_resultfilename);
		std::ofstream resultStream;
		resultStream.open(resultfile.c_str(), std::ios::app);

		auto formatOutput = [&](auto value) {
			if (_outputFormat == "tab") {
				resultStream << std::setw(_writeWidth) << std::scientific << std::setprecision(_writePrecision) << value;
			} else {  // csv
				resultStream << "," << std::scientific << std::setprecision(_writePrecision) << value;
			}
		};

		if (_outputFormat == "tab") {
			// Do not write simstep in scientific notation
			resultStream << std::setw(_writeWidth) << simstep;
		} else {
			// Do not write comma first
			resultStream << simstep;
		}
		formatOutput(_simulation.getSimulationTime());
		formatOutput(domain->getGlobalUpot());
		formatOutput(_U_pot_acc->getAverage());
		formatOutput(ekin);
		formatOutput(_U_kin_acc->getAverage());
		formatOutput(domain->getGlobalPressure());
		formatOutput(_p_acc->getAverage());
		formatOutput(domain->getGlobalBetaTrans());
		formatOutput(domain->getGlobalBetaRot());
		formatOutput(cv);
		formatOutput(globalNumMolecules);
		resultStream << std::endl;
		resultStream.close();
	}
}
