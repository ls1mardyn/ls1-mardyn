#include "io/ResultWriter.h"

#include <fstream>
#include <iomanip>

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
		error_message << "[ResultWriter] Wrong output format specified! Use \"csv\" or \"tab\" instead of \"" << _outputFormat << "\"" << std::endl;
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

		// Write header of output file in case of "tab"
		if (_outputFormat == "tab") {
			resultStream << "# Averages are the accumulated values over " << _U_pot_acc->getWindowLength()  << " time steps."<< std::endl;
		}
		
		if (_outputFormat == "tab") {
			// Do not write simstep in scientific notation
			resultStream << std::setw(_writeWidth) << "simstep";
		} else {
			// Do not write comma first
			resultStream << "simstep";
		}
		formatOutput(resultStream, "time");
		formatOutput(resultStream, "U_pot");
		formatOutput(resultStream, "U_pot_avg");
		formatOutput(resultStream, "U_kin");
		formatOutput(resultStream, "U_kin_avg");
		formatOutput(resultStream, "p");
		formatOutput(resultStream, "p_avg");
		formatOutput(resultStream, "N");
		resultStream << std::endl;
		resultStream.close();
	}
}

void ResultWriter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                           unsigned long simstep) {

	// Writing of cavities now handled by CavityWriter

	const unsigned long globalNumMolecules = domain->getglobalNumMolecules(true, particleContainer, domainDecomp);
	const double ekin = domain->getGlobalUkinTrans()+domain->getGlobalUkinRot();

	_U_pot_acc->addEntry(domain->getGlobalUpot());
	_U_kin_acc->addEntry(ekin);
	_p_acc->addEntry(domain->getGlobalPressure());
	if ((domainDecomp->getRank() == 0) && (simstep % _writeFrequency == 0)){
		const std::string resultfile(_resultfilename);
		std::ofstream resultStream;
		resultStream.open(resultfile.c_str(), std::ios::app);

		if (_outputFormat == "tab") {
			// Do not write simstep in scientific notation
			resultStream << std::setw(_writeWidth) << simstep;
		} else {
			// Do not write comma first
			resultStream << simstep;
		}
		formatOutput(resultStream, _simulation.getSimulationTime());
		formatOutput(resultStream, domain->getGlobalUpot());
		formatOutput(resultStream, _U_pot_acc->getAverage());
		formatOutput(resultStream, ekin);
		formatOutput(resultStream, _U_kin_acc->getAverage());
		formatOutput(resultStream, domain->getGlobalPressure());
		formatOutput(resultStream, _p_acc->getAverage());
		formatOutput(resultStream, globalNumMolecules);
		resultStream << std::endl;
		resultStream.close();
	}
}

template <typename T>
void ResultWriter::formatOutput(std::ostream& os, T value) {
	if (_outputFormat == "tab") {
		os << std::setw(_writeWidth) << std::scientific << std::setprecision(_writePrecision) << value;
	} else {  // csv
		os << "," << std::scientific << std::setprecision(_writePrecision) << value;
	}
};
