#include "io/ResultWriter.h"

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "utils/Logger.h"

using Log::global_log;
using namespace std;


void ResultWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << endl;

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << endl;

	size_t acc_steps = 1000;
	xmlconfig.getNodeValue("accumulation_steps", acc_steps);
	_U_pot_acc = new Accumulator<double>(acc_steps);
	_p_acc = new Accumulator<double>(acc_steps);
	global_log->info() << "Accumulation steps: " << acc_steps << endl;
}

void ResultWriter::init(ParticleContainer * /*particleContainer*/,
                        DomainDecompBase *domainDecomp, Domain * /*domain*/) {
	// initialize result file
	string resultfile(_outputPrefix+".res");
	time_t now;
	time(&now);
	if(domainDecomp->getRank() == 0) {
		_resultStream.open(resultfile.c_str());
		_resultStream << "# ls1 MarDyn simulation started at " << ctime(&now) << endl;
		_resultStream << "# Averages are the accumulated values over " << _U_pot_acc->getWindowLength()  << " time steps."<< endl;
		_resultStream << std::setw(10) << "#step" << std::setw(28) << "time" << std::setw(28) << "U_pot" << std::setw(28) << "U_pot_avg" << std::setw(28)
			<< "p" << std::setw(28) << "p_avg" << std::setw(28) << "beta_trans" << std::setw(28) << "beta_rot" << std::setw(28) << "c_v" << std::setw(28) << "N" << "\n";
	}
}

void ResultWriter::endStep(ParticleContainer * /*particleContainer*/, DomainDecompBase *domainDecomp, Domain *domain,
                           unsigned long simstep) {

	// Writing of cavities now handled by CavityWriter

	_U_pot_acc->addEntry(domain->getGlobalUpot());
	_p_acc->addEntry(domain->getGlobalPressure());
	if((domainDecomp->getRank() == 0) && (simstep % _writeFrequency == 0)){
		_resultStream << std::setw(10) << simstep << std::setw(28) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _simulation.getSimulationTime()
			<< std::setw(28) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << domain->getGlobalUpot()
			<< std::setw(28) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _U_pot_acc->getAverage()
			<< std::setw(28) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << domain->getGlobalPressure()
			<< std::setw(28) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _p_acc->getAverage()
			<< std::setw(28) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << domain->getGlobalBetaTrans()
			<< std::setw(28) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << domain->getGlobalBetaRot()
			<< std::setw(28) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << domain->cv()
			<< std::setw(28) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << domain->getglobalNumMolecules();
		_resultStream << "\n";
	}
}

void ResultWriter::finish(ParticleContainer * /*particleContainer*/,
						  DomainDecompBase * /*domainDecomp*/, Domain * /*domain*/){
	time_t now;
	time(&now);
	_resultStream << "# ls1 mardyn simulation finished at " << ctime(&now) << endl;
	_resultStream.close();
}
