#include "EnergyRAPL.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#ifdef ENABLE_MPI
#include <map>
#endif

#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"

EnergyRAPL::RAPLCounter::RAPLCounter(const std::string& domainBasePath) {
	// File path for reading current micro joules
	std::ostringstream microJoulesPath;
	microJoulesPath << domainBasePath << "/energy_uj";
	_microJoulesPath = microJoulesPath.str();
	// Range, i.e., maximum value of RAPL energy counter, in micro-joules
	std::ostringstream rangeMicroJoulesPath;
	rangeMicroJoulesPath << domainBasePath << "/max_energy_range_uj";
	std::ifstream rangeMicroJoulesFile(rangeMicroJoulesPath.str());
	rangeMicroJoulesFile >> _rangeMicroJoules;
	reset();
}

void EnergyRAPL::RAPLCounter::reset() {
	// Update last micro joules
	update();
	_microJoules = 0;
}

double EnergyRAPL::RAPLCounter::update() {
	long long currentMicroJoules;
	std::ifstream packageIdFile(_microJoulesPath);
	packageIdFile >> currentMicroJoules;
	long long deltaMicroJoules = currentMicroJoules - _lastMicroJoules;
	// Correct counter overflow (occurs around every 60 seconds)
	if (0 > deltaMicroJoules) {
		deltaMicroJoules += _rangeMicroJoules;
	}
	_lastMicroJoules = currentMicroJoules;
	_microJoules += deltaMicroJoules;
	return static_cast<double>(_microJoules) * 1e-6;  // Convert micro joules to joules
}

int EnergyRAPL::getNumberOfPackages() {
	int maxPackageIdx = 0;
	for (int cpuIdx = 0;; cpuIdx++) {
		std::ostringstream packageIdxPath;
		packageIdxPath << "/sys/devices/system/cpu/cpu" << cpuIdx << "/topology/physical_package_id";
		std::ifstream packageIdFile(packageIdxPath.str());
		if (!packageIdFile.good()) {
			break;  // Stop if there are no more CPUs (packages) on this machine (file does not exist)
		}
		int packageIdx = -1;  // No more packages
		packageIdFile >> packageIdx;
		if (packageIdx > maxPackageIdx) {
			maxPackageIdx = packageIdx;
		}
	}
	return maxPackageIdx + 1;
}

void EnergyRAPL::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
#ifdef ENABLE_MPI
	int processorNameLength;
	MPI_Get_processor_name(_processorName, &processorNameLength);
	MPI_Comm_rank(MPI_COMM_WORLD, &_thisRank);
#endif
	if (!_outputprefix.empty()) {
		std::ostringstream outputFilename;
		outputFilename << _outputprefix << ".tsv";
		std::ofstream outputFile(outputFilename.str().c_str());
		outputFile << "milliseconds\tsimstep\tjoules" << std::endl;
	}
	// For each package...
	const int numberOfPackages = getNumberOfPackages();
	for (int packageIdx = 0; packageIdx < numberOfPackages; packageIdx++) {
		std::ostringstream packageBasePath;
		packageBasePath << _basePathRAPL << "intel-rapl:" << packageIdx;
		// add package domain and...
		Log::global_log->info() << "[" << getPluginName() << "] Adding package domain " << packageBasePath.str()
								<< std::endl;
		_counters.push_back(RAPLCounter(packageBasePath.str()));
		// for each domain in package...
		for (int domainIdx = 0;; domainIdx++) {
			std::ostringstream domainBasePath;
			domainBasePath << packageBasePath.str() << "/intel-rapl:" << packageIdx << ":" << domainIdx;
			std::ostringstream domainNamePath;
			domainNamePath << domainBasePath.str() << "/name";
			std::string domainName;
			std::ifstream domainNameFile(domainNamePath.str());
			if (!domainNameFile.good()) {
				// (if domain exists...)
				break;  // Stop if there are no more domains for that CPU (package)
			}
			domainNameFile >> domainName;
			if (0 != domainName.compare("dram")) {
				continue;  // (and is DRAM domain)
			}
			// add domain
			Log::global_log->info() << "[" << getPluginName() << "] Adding DRAM domain " << packageBasePath.str()
									<< std::endl;
			_counters.push_back(RAPLCounter(packageBasePath.str()));
		}
	}
	_simstart = std::chrono::steady_clock::now();
	for (auto counter : _counters) {
		counter.reset();
	}
}

void EnergyRAPL::endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
						 unsigned long simstep) {
	_joules = 0;
	for (auto counter : _counters) {
		_joules += counter.update();
	}
	_simstep = simstep;
	if (0 < _writeFrequency && simstep % _writeFrequency == 0) {
		outputEnergyJoules();
	}
}

void EnergyRAPL::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	xmlconfig.getNodeValue("outputprefix", _outputprefix);
}

void EnergyRAPL::outputEnergyJoules() {
	double joules = _joules;
#ifdef ENABLE_MPI
	// Collect results from all nodes (matching of separate messages over tag)
	if (_thisRank == 0) {
		std::map<std::string, double> nodeJoules;
		nodeJoules[_processorName] = joules;  // Store result of rank 0
		int numberOfRanks;
		MPI_Comm_size(MPI_COMM_WORLD, &numberOfRanks);
		for (int otherRank = 1; otherRank < numberOfRanks; otherRank++) {
			MPI_Status status;
			char processorName[MPI_MAX_PROCESSOR_NAME];
			MPI_Recv(&processorName[0], MPI_MAX_PROCESSOR_NAME, MPI_CHAR, otherRank, /* tag */ otherRank,
					 MPI_COMM_WORLD, &status);
			MPI_Recv(&joules, 1, MPI_DOUBLE, otherRank, /* tag */ otherRank, MPI_COMM_WORLD, &status);
			nodeJoules[_processorName] = joules;
		}
		joules = 0;
		for (const auto [_, value] : nodeJoules) {
			joules += value;
		}
	} else {
		MPI_Send(_processorName, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, /* dest */ 0, /* tag */ _thisRank, MPI_COMM_WORLD);
		MPI_Send(&joules, 1, MPI_DOUBLE, /* dest */ 0, /* tag */ _thisRank, MPI_COMM_WORLD);
		return;
	}
#endif
	if (_outputprefix.empty()) {
		Log::global_log->info() << "Simstep = " << _simstep << "\tEnergy consumed = " << joules << " J" << std::endl;
	} else {
		std::ostringstream outputFilename;
		outputFilename << _outputprefix << ".tsv";
		std::ofstream outputFile(outputFilename.str().c_str(), std::ios_base::app);
		int64_t milliseconds =
			std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - _simstart).count();
		outputFile << milliseconds << "\t" << _simstep << "\t" << joules << std::endl;
	}
}

void EnergyRAPL::finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
	if (0 == _writeFrequency) {
		outputEnergyJoules();
	}
}
