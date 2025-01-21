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

#pragma once

#include <filesystem>
#include <memory>
#include <string>

#include "plugins/PluginBase.h"
#include "utils/Accumulator.h"

/** @brief Writes thermodynamic properties to a file in CSV format.
 *
 * Writes the current value and the average over a specified number of time
 * steps of values to a CSV file including a column header line as well as
 * additional information as comments in lines starting with '#'.
 *
 * The following fields are written to the file:
 * - simstep: time step
 * - N: number of particles
 * - time: simulation time
 * - U_pot: potential energy
 * - U_pot_avg: average potential energy
 * - U_kin: kinetic energy
 * - U_kin_avg: average kinetic energy
 * - p: pressure
 * - p_avg: average pressure
 */
class CsvWriter : public PluginBase {

public:
	/// @brief Read in XML configuration for CsvWriter and all its included objects.
	///
	/// The following xml object structure is handled by this method:
	// clang-format off
	/// \code{.xml}
	///	<outputplugin name="CSVWriter">
	///  	<writefrequency>INTEGER</writefrequency>          <!-- Frequency in which the output is written; Default: 1000 -->
	///		<outputprefix>STRING</outputprefix>               <!-- Prefix of the output file; Default: "results" -->
	///		<accumulation_steps>INTEGER</accumulation_steps>  <!-- Result is accumulated over the specified steps; Default: 1000 -->
	///		<writeprecision>UINTEGER</writeprecision>         <!-- Precision of output can be set here; Default: 5 -->
	///	</outputplugin>
	/// \endcode
	///
	// clang-format on
	virtual void readXML(XMLfileUnits &xmlconfig);

	void init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain);

	void endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
				 unsigned long simstep);

	void finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain);

	std::string getPluginName() { return std::string("CSVWriter"); }
	static PluginBase *createInstance() { return new CsvWriter(); }

private:
	long _writeFrequency{1000L};
	unsigned int _writePrecision{5};
	std::string _outputPrefix{"results"};
	std::filesystem::path _filename;
	std::unique_ptr<Accumulator<double>> _U_pot_acc;
	std::unique_ptr<Accumulator<double>> _U_kin_acc;
	std::unique_ptr<Accumulator<double>> _p_acc;
};
