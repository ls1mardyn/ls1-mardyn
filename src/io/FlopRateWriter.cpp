/*
 * FlopRateWriter.cpp
 *
 *  Created on: 8 May 2017
 *      Author: tchipevn
 */

#include "FlopRateWriter.h"

#include "particleContainer/ParticleContainer.h"
#include "particleContainer/adapter/FlopCounter.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"

void FlopRateWriter::readXML(XMLfileUnits& xmlconfig) {
	std::string mode;
	xmlconfig.getNodeValue("mode", mode);
	if(mode.compare("stdout") == 0) {
		_writeToStdout = true;
		_writeToFile = false;
	} else if(mode.compare("file") == 0) {
		_writeToStdout = false;
		_writeToFile = true;
	} else if(mode.compare("both") == 0) {
		_writeToStdout = true;
		_writeToFile = true;
	} else {
		global_log->error() << "Unknown FlopRateOutputPlugin::mode. Choose \"stdout\", \"file\" or \"both\"." << endl;
	}

	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << endl;

	// TODO:
	if(_writeToFile) {
		global_log->error() << "TODO: file output not yet supported." << endl;
		global_simulation->exit(1);
	}

	if(_writeToFile) {
		_outputPrefix = "mardyn";
		xmlconfig.getNodeValue("outputprefix", _outputPrefix);
		global_log->info() << "Output prefix: " << _outputPrefix << endl;
	}
}

void FlopRateWriter::initOutput(ParticleContainer* /*particleContainer*/,
		DomainDecompBase* domainDecomp, Domain* domain) {

	if(_writeToFile != true)
		return;

	// initialize result file
	string resultfile(_outputPrefix+".res");
	time_t now;
	time(&now);
	if(domainDecomp->getRank()==0){
		_fileStream.open(resultfile.c_str());
		_fileStream << "# ls1 MarDyn simulation started at " << ctime(&now) << endl;
		_fileStream << "#step\tt\t\tFLOP-Count\tFLOP-Rate-force\t\tFLOP-Rate-loop\tefficiency(%)\t\n";
	}
}

void FlopRateWriter::doOutput(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* /*domain*/, unsigned long simstep,
		std::list<ChemicalPotential>* /*lmu*/,
		std::map<unsigned, CavityEnsemble>* /*mcav*/) {

	if (not ((_writeFrequency == 1 and simstep > 0) or simstep % _writeFrequency == 1)) {
		return;
	}

	double flops = _flopCounter->getTotalFlopCount();

	unsigned long numElapsedIterations = global_simulation->timers()->getNumElapsedIterations();
	double force_calculation_time = global_simulation->timers()->getTime("SIMULATION_FORCE_CALCULATION") / numElapsedIterations;
	global_simulation->timers()->stop("SIMULATION_LOOP");
	double loop_time = global_simulation->timers()->getTime("SIMULATION_LOOP") / numElapsedIterations ;
	global_simulation->timers()->start("SIMULATION_LOOP");

	double flop_rate_force = flops / force_calculation_time / (1024*1024);
	double flop_rate_loop = flops / loop_time / (1024*1024);

	if(_writeToStdout) {
		global_log->info() << "FlopRateWriter (simulation step " << simstep << ")" << endl
			<< "\tFLOP-Count per Iteration: " << flops << " FLOPs" << endl
			<< "\tFLOP-rate in force calculation: " << flop_rate_force << " MFLOPS" << endl
			<< "\tFLOP-rate for main loop       : " << flop_rate_loop << " MFLOPS (" << flop_rate_loop / flop_rate_force * 100. << " %)" << endl;
	}


#if 0
	if((domainDecomp->getRank() == 0) && (simstep % _writeFrequency == 0)){
		_resultStream << simstep << "\t" << _simulation.getSimulationTime()
		              << "\t\t" << domain->getGlobalUpot() << "\t" << _U_pot_acc->getAverage()
					  << "\t\t" << domain->getGlobalPressure() << "\t" << _p_acc->getAverage()
		              << "\t\t" << domain->getGlobalBetaTrans() << "\t" << domain->getGlobalBetaRot()
		              << "\t\t" << domain->cv() << "\t\t" << domain->getglobalNumMolecules();

                map<unsigned, CavityEnsemble>::iterator ceit;
                for(ceit = mcav->begin(); ceit != mcav->end(); ceit++)
                {
                   _resultStream << "\t" << ceit->second.numCavities();
                }
                _resultStream << "\n";
	}
#endif
}

void FlopRateWriter::finishOutput(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain) {

	if(_writeToFile != true)
		return;

	time_t now;
	time(&now);
	_fileStream << "# ls1 mardyn simulation finished at " << ctime(&now) << endl;
	_fileStream << "# \n# Please address your questions and suggestions to the ls1 mardyn contact point:\n# \n# E-mail: contact@ls1-mardyn.de\n# \n# Phone: +49 631 205 3227\n# University of Kaiserslautern\n# Computational Molecular Engineering\n# Erwin-Schroedinger-Str. 44\n# D-67663 Kaiserslautern, Germany\n# \n# http://www.ls1-mardyn.de/\n";

	_fileStream.close();
}

void FlopRateWriter::measureFLOPS(ParticleContainer* particleContainer, unsigned long simstep) {
	if(simstep % _writeFrequency != 0) {
		return;
	}
	if(_flopCounter == nullptr) {
		_flopCounter = new FlopCounter(global_simulation->getcutoffRadius(), global_simulation->getLJCutoff());
	}
	particleContainer->traverseCells(*_flopCounter);
}
