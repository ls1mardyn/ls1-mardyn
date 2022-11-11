#include "io/ResultWriter.h"

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include <chrono>

using Log::global_log;

void ResultWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "[ResultWriter] Write frequency: " << _writeFrequency << endl;

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "[ResultWriter] Output prefix: " << _outputPrefix << endl;

	size_t acc_steps = 1000;
	xmlconfig.getNodeValue("accumulation_steps", acc_steps);
	_U_pot_acc = new Accumulator<double>(acc_steps);
	_p_acc = new Accumulator<double>(acc_steps);
	global_log->info() << "[ResultWriter] Accumulation steps: " << acc_steps << endl;

	_writePrecision = 5;
	xmlconfig.getNodeValue("writeprecision", _writePrecision);
	global_log->info() << "[ResultWriter] Write precision: " << _writePrecision << endl;
}

void ResultWriter::init(ParticleContainer * /*particleContainer*/,
                        DomainDecompBase *domainDecomp, Domain * /*domain*/) {

	if(domainDecomp->getRank() == 0) {
		const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
		tm unused{};
		const auto nowStr = std::put_time(localtime_r(&now, &unused), "%c");
		string resultfile(_outputPrefix+".res");
		_resultStream.open(resultfile.c_str(), std::ios::out);
		_resultStream << "# ls1 MarDyn simulation started at " << nowStr << endl;
		_resultStream << "# Averages are the accumulated values over " << _U_pot_acc->getWindowLength()  << " time steps."<< endl;
		_resultStream << std::setw(10) << "# step" << std::setw(_writePrecision+15) << "time" 
			<< std::setw(_writePrecision+15) << "U_pot"
			<< std::setw(_writePrecision+15) << "U_pot_avg"
			<< std::setw(_writePrecision+15) << "p"
			<< std::setw(_writePrecision+15) << "p_avg"
			<< std::setw(_writePrecision+15) << "beta_trans"
			<< std::setw(_writePrecision+15) << "beta_rot"
			<< std::setw(_writePrecision+15) << "c_v"
			<< std::setw(_writePrecision+15) << "N"
			<< endl;
	}
}

void ResultWriter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
                           unsigned long simstep) {

	// Writing of cavities now handled by CavityWriter

	unsigned long globalNumMolecules = domain->getglobalNumMolecules(true, particleContainer, domainDecomp);
	double cv = domain->cv();

	_U_pot_acc->addEntry(domain->getGlobalUpot());
	_p_acc->addEntry(domain->getGlobalPressure());
	if((domainDecomp->getRank() == 0) && (simstep % _writeFrequency == 0)){
		_resultStream << std::setw(10) << simstep << std::setw(_writePrecision+15) << std::scientific << std::setprecision(_writePrecision) << _simulation.getSimulationTime()
			<< std::setw(_writePrecision+15) << std::scientific << std::setprecision(_writePrecision) << domain->getGlobalUpot()
			<< std::setw(_writePrecision+15) << std::scientific << std::setprecision(_writePrecision) << _U_pot_acc->getAverage()
			<< std::setw(_writePrecision+15) << std::scientific << std::setprecision(_writePrecision) << domain->getGlobalPressure()
			<< std::setw(_writePrecision+15) << std::scientific << std::setprecision(_writePrecision) << _p_acc->getAverage()
			<< std::setw(_writePrecision+15) << std::scientific << std::setprecision(_writePrecision) << domain->getGlobalBetaTrans()
			<< std::setw(_writePrecision+15) << std::scientific << std::setprecision(_writePrecision) << domain->getGlobalBetaRot()
			<< std::setw(_writePrecision+15) << std::scientific << std::setprecision(_writePrecision) << cv
			<< std::setw(_writePrecision+15) << std::scientific << std::setprecision(_writePrecision) << globalNumMolecules
			<< endl;
	}
}

void ResultWriter::finish(ParticleContainer * /*particleContainer*/,
						  DomainDecompBase *domainDecomp, Domain * /*domain*/){

	if(domainDecomp->getRank() == 0) {
		const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
		tm unused{};
		const auto nowStr = std::put_time(localtime_r(&now, &unused), "%c");
		_resultStream << "# ls1 mardyn simulation finished at " << nowStr << endl;
		_resultStream.close();
	}
}
