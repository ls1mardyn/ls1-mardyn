/*
 * Copyright (C) 2024       High Performance Computing Center Stuttgart,
 *                          University of Stuttgart.  All rights reserved.
 *
 * This file is part of ls1-mardyn. Redistribution and use in source and
 * binary forms, with or without modification, are permitted provided that
 * the following conditions are met:
 * [Conditions]
 *
 * For full license details, see the LICENSE file in the root of this project.
 *
 */

#include "io/CsvWriter.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>

#include "Domain.h"
#include "Simulation.h"
#include "parallel/DomainDecompBase.h"
#include "utils/Logger.h"
#include "utils/mardyn_assert.h"

void CsvWriter::readXML(XMLfileUnits &xmlconfig) {
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	Log::global_log->info() << "[CSVWriter] Write frequency: " << _writeFrequency << std::endl;
	if (_writeFrequency <= 0) {
		std::string error_message =
			"[CSVWriter] Write frequency must be a positive nonzero integer, but is " + std::to_string(_writeFrequency);
		MARDYN_EXIT(error_message);
	}

	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	Log::global_log->info() << "[CSVWriter] Output prefix: " << _outputPrefix << std::endl;
	_filename = std::filesystem::path(_outputPrefix + ".csv");
	Log::global_log->info() << "[CSVWriter] Output file: " << _filename << std::endl;

	size_t acc_steps = 1000;
	xmlconfig.getNodeValue("accumulation_steps", acc_steps);
	Log::global_log->info() << "[CSVWriter] Accumulation steps: " << acc_steps << std::endl;
	_U_pot_acc = std::make_unique<Accumulator<double>>(acc_steps);
	_U_kin_acc = std::make_unique<Accumulator<double>>(acc_steps);
	_p_acc = std::make_unique<Accumulator<double>>(acc_steps);

	xmlconfig.getNodeValue("writeprecision", _writePrecision);
	Log::global_log->info() << "[CSVWriter] Write precision: " << _writePrecision << std::endl;
}

void CsvWriter::init(ParticleContainer * /*particleContainer*/, DomainDecompBase *domainDecomp, Domain * /*domain*/) {
	if (domainDecomp->getRank() == 0) {
		std::ofstream csvFile(_filename);
		const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
		tm unused{};
		const auto nowStr = std::put_time(localtime_r(&now, &unused), "%FT%T%z");  // ISO 8601 time format
		csvFile << "# ls1-mardyn simulation started at " << nowStr << std::endl;
		csvFile << "#" << std::endl;
		csvFile << "# Fields in this file are" << std::endl;
		csvFile << "# simstep: current time step" << std::endl;
		csvFile << "# N: number of particles in the system at current time step" << std::endl;
		csvFile << "# time: current simulation time" << std::endl;
		csvFile << "# U_pot: potential energy in the system at current time step" << std::endl;
		csvFile << "# U_pot_avg: average potential energy in the system over the last " << _U_pot_acc->getWindowLength()
				<< " time steps" << std::endl;
		csvFile << "# U_kin: kinetic energy in the system at current time step" << std::endl;
		csvFile << "# U_kin_avg: kinetic energy in the system over the last " << _U_kin_acc->getWindowLength()
				<< " time steps" << std::endl;
		csvFile << "# p: pressure in the system at current time step" << std::endl;
		csvFile << "# p_avg: average pressure in the system  over the last " << _p_acc->getWindowLength()
				<< " time steps" << std::endl;
		csvFile << "#" << std::endl;
		csvFile << "simstep,N,time,U_pot,U_pot_avg,U_kin,U_kin_avg,p,p_avg" << std::endl;
	}
}

void CsvWriter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
						unsigned long simstep) {
	const auto globalNumMolecules = domain->getglobalNumMolecules(true, particleContainer, domainDecomp);
	const auto U_pot = domain->getGlobalUpot();
	const auto U_kin = domain->getGlobalUkinTrans() + domain->getGlobalUkinRot();
	const auto p = domain->getGlobalPressure();

	_U_pot_acc->addEntry(U_pot);
	_U_kin_acc->addEntry(U_kin);
	_p_acc->addEntry(p);

	if ((domainDecomp->getRank() == 0) && (simstep % _writeFrequency == 0)) {
		std::ofstream csvFile(_filename, std::ios::app);

		csvFile << simstep << ",";
		csvFile << globalNumMolecules << ",";

		csvFile << std::scientific << std::setprecision(_writePrecision);
		csvFile << _simulation.getSimulationTime() << ",";
		csvFile << U_pot << ",";
		csvFile << _U_pot_acc->getAverage() << ",";
		csvFile << U_kin << ",";
		csvFile << _U_kin_acc->getAverage() << ",";
		csvFile << p << ",";
		csvFile << _p_acc->getAverage();
		csvFile << std::endl;
	}
}

void CsvWriter::finish(ParticleContainer * /*particleContainer*/, DomainDecompBase *domainDecomp, Domain * /*domain*/) {
	if (domainDecomp->getRank() == 0) {
		std::ofstream csvFile(_filename, std::ios::app);
		const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
		tm unused{};
		const auto nowStr = std::put_time(localtime_r(&now, &unused), "%FT%T%z");  // ISO 8601 time format
		csvFile << "# ls1-mardyn simulation finished at " << nowStr << std::endl;
	}
}
