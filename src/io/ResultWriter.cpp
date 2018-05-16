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
		_resultStream << "#step\tt\t\tU_pot\tU_pot_avg\t\tp\tp_avg\t\tbeta_trans\tbeta_rot\t\tc_v\t\tN\t(N_cav*)\n";
	}
}

void ResultWriter::endStep(ParticleContainer * /*particleContainer*/, DomainDecompBase *domainDecomp, Domain *domain,
                           unsigned long simstep, list<ChemicalPotential> * /*lmu*/,
                           map<unsigned, CavityEnsemble> *mcav) {
	_U_pot_acc->addEntry(domain->getGlobalUpot());
	_p_acc->addEntry(domain->getGlobalPressure());
	if((domainDecomp->getRank() == 0) && (simstep % _writeFrequency == 0)){
		_resultStream << simstep << "\t" << _simulation.getSimulationTime()
			<< "\t\t" << domain->getGlobalUpot() << "\t" << _U_pot_acc->getAverage()
			<< "\t\t" << domain->getGlobalPressure() << "\t" << _p_acc->getAverage()
			<< "\t\t" << domain->getGlobalBetaTrans() << "\t" << domain->getGlobalBetaRot()
			<< "\t\t" << domain->cv() << "\t\t" << domain->getglobalNumMolecules();
		map<unsigned, CavityEnsemble>::iterator ceit;
		for(ceit = mcav->begin(); ceit != mcav->end(); ceit++) {
			_resultStream << "\t" << ceit->second.numCavities();
		}
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
