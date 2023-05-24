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
		global_log->error() << "Unknown FlopRateOutputPlugin::mode. Choose \"stdout\", \"file\" or \"both\"." << std::endl;
	}

	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << std::endl;

	// TODO:
	if(_writeToFile) {
		global_log->error() << "TODO: file output not yet supported." << std::endl;
		Simulation::exit(1);
	}

	if(_writeToFile) {
		_outputPrefix = "mardyn";
		xmlconfig.getNodeValue("outputprefix", _outputPrefix);
		global_log->info() << "Output prefix: " << _outputPrefix << std::endl;
	}
}

void FlopRateWriter::init(ParticleContainer * /*particleContainer*/,
                          DomainDecompBase *domainDecomp, Domain *domain) {

	if(_writeToFile != true)
		return;

	// initialize result file
	std::string resultfile(_outputPrefix+".res");
	const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	tm unused{};
	if(domainDecomp->getRank()==0){
		_fileStream.open(resultfile.c_str());
		_fileStream << "# ls1 MarDyn simulation started at " << std::put_time(localtime_r(&now, &unused), "%c") << std::endl;
		_fileStream << "#step\tt\t\tFLOP-Count\tFLOP-Rate-force\t\tFLOP-Rate-loop\tefficiency(%)\t\n";
	}
}

void FlopRateWriter::afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		unsigned long simstep){
	global_log->debug()  << "[FLOPRATEWRITER] after forces FLOPs" << std::endl;
	measureFLOPS(particleContainer, simstep);
}

void FlopRateWriter::endStep(ParticleContainer *particleContainer,
                             DomainDecompBase *domainDecomp, Domain * /*domain*/, unsigned long simstep
                             ) {

	if (not ((_writeFrequency == 1 and simstep > 0) or simstep % _writeFrequency == 1)) {
		return;
	}

	double flops = _flopCounter->getTotalFlopCount();

	unsigned long numElapsedIterations = global_simulation->timers()->getNumElapsedIterations();
	double force_calculation_time = global_simulation->timers()->getTime("SIMULATION_FORCE_CALCULATION") / numElapsedIterations;
	global_simulation->timers()->stop("SIMULATION_LOOP");
	double loop_time = global_simulation->timers()->getTime("SIMULATION_LOOP") / numElapsedIterations ;
	global_simulation->timers()->start("SIMULATION_LOOP");

	double flop_rate_force = flops / force_calculation_time;
	double flop_rate_loop = flops / loop_time;
	double percentage = flop_rate_loop / flop_rate_force * 100. ;

	// compute convenient prefixes, i.e. kilo, mega, giga, tera, peta, exa
	char prefix_flops = '0';
	double flops_normalized;
	setPrefix(flops, flops_normalized, prefix_flops);

	char prefix_flop_rate_force = '0';
	double flop_rate_force_normalized;
	setPrefix(flop_rate_force, flop_rate_force_normalized, prefix_flop_rate_force);

	char prefix_flop_rate_loop = '0';
	double flop_rate_loop_normalized;
	setPrefix(flop_rate_loop, flop_rate_loop_normalized, prefix_flop_rate_loop);

	if(_writeToStdout) {
		global_log->info() << "FlopRateWriter (simulation step " << simstep << ")" << std::endl
			<< "\tFLOP-Count per Iteration           : " << flops_normalized << " " << prefix_flops << "FLOPs" << std::endl
			<< "\tFLOP-rate in force calculation     : " << flop_rate_force_normalized << " " << prefix_flop_rate_force << "FLOP/sec" << std::endl
			<< "\tFLOP-rate for main loop            : " << flop_rate_loop_normalized << " " << prefix_flop_rate_loop << "FLOP/sec (" << percentage << " %)" << std::endl;
		_flopCounter->printStats();
	}

}

void FlopRateWriter::finish(ParticleContainer *particleContainer,
							DomainDecompBase *domainDecomp, Domain *domain) {

	delete _flopCounter;

	if(_writeToFile != true)
		return;

	const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	tm unused{};
	_fileStream << "# ls1 mardyn simulation finished at " << std::put_time(localtime_r(&now, &unused), "%c") << std::endl;
	_fileStream << "# \n# Please address your questions and suggestions to the ls1 mardyn contact point:\n# \n# E-mail: contact@ls1-mardyn.de\n# \n# Phone: +49 631 205 3227\n# University of Kaiserslautern\n# Computational Molecular Engineering\n# Erwin-Schroedinger-Str. 44\n# D-67663 Kaiserslautern, Germany\n# \n# http://www.ls1-mardyn.de/\n";

	_fileStream.close();
}

void FlopRateWriter::measureFLOPS(ParticleContainer* particleContainer, unsigned long simstep) {
	if ((simstep - 1) % _writeFrequency != 0) {
		return;
	}
	if(_flopCounter == nullptr) {
		_flopCounter = new FlopCounter(global_simulation->getcutoffRadius(), global_simulation->getLJCutoff());
	}
	particleContainer->traverseCells(*_flopCounter);
}

void FlopRateWriter::setPrefix(double f_in, double& f_out, char& prefix) const {
//	// powers of two or powers of ten? Update: resolved.

	// Powers of ten, as stated by IEC, NIST and ISO, cf. Wikipedia on Binary_prefix
	double kilo = 1e+3;
	double mega = 1e+6;
	double giga = 1e+9;
	double tera = 1e+12;
	double peta = 1e+15;
	double exa  = 1e+18;

	if (f_in >= exa) {
		f_out = f_in / exa;
		prefix = 'E';
	} else if (f_in >= peta) {
		f_out = f_in / peta;
		prefix = 'P';
	} else if (f_in >= tera) {
		f_out = f_in / tera;
		prefix = 'T';
	} else if (f_in >= giga) {
		f_out = f_in / giga;
		prefix = 'G';
	} else if (f_in >= mega) {
		f_out = f_in / mega;
		prefix = 'M';
	} else if (f_in >= kilo) {
		f_out = f_in / kilo;
		prefix = 'K';
	}
}
