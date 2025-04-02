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
	// File path for reading current micro joule
	std::ostringstream microJoulePath;
	microJoulePath << domainBasePath << "/energy_uj";
	_microJoulePath = microJoulePath.str();
	// Range, i.e., maximum value of RAPL energy counter, in micro-joule
	std::ostringstream rangeMicroJoulePath;
	rangeMicroJoulePath << domainBasePath << "/max_energy_range_uj";
	std::ifstream rangeMicroJouleFile(rangeMicroJoulePath.str());
	if(!rangeMicroJouleFile) {
		Log::global_log->error() << "Failed to open RAPL counter file " << rangeMicroJoulePath.str() << std::endl;
		exit(1);
	}
	rangeMicroJouleFile >> _rangeMicroJoule;
	reset();
}

void EnergyRAPL::RAPLCounter::reset() {
	// Update last micro joule
	update();
	_microJoule = 0;
	_joule = 0;
}

double EnergyRAPL::RAPLCounter::update() {
	long long currentMicroJoule;
	std::ifstream microJouleFile(_microJoulePath);
	if(!microJouleFile) {
		Log::global_log->error() << "Failed to open RAPL counter file " << _microJoulePath << std::endl;
		exit(1);
	}
	microJouleFile >> currentMicroJoule;
	long long deltaMicroJoule = currentMicroJoule - _lastMicroJoule;
	// Correct counter overflow (occurs around every 60 seconds)
	if (0 > deltaMicroJoule) {
		deltaMicroJoule += _rangeMicroJoule;
	}
	_lastMicroJoule = currentMicroJoule;
	_microJoule += deltaMicroJoule;
	// Avoid overflow of _microJoule
	const long long microJoulePerJoule = 1000000L;
	_joule += _microJoule / microJoulePerJoule;
	_microJoule %= microJoulePerJoule;
	return static_cast<double>(_microJoule) / microJoulePerJoule + _joule;
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

std::string EnergyRAPL::getOutputFilename() {
	if (_outputprefix.size() == 0) {
		return std::string();
	}
	std::ostringstream outputFilename;
	outputFilename << _outputprefix;
#ifdef ENABLE_MPI
	if(_per_host) {
		outputFilename << "_" << _processorName;
	}
#endif
	outputFilename << ".tsv";
	return outputFilename.str();
}

void EnergyRAPL::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
#ifdef ENABLE_MPI
	int processorNameLength;
	MPI_Get_processor_name(_processorName, &processorNameLength);
	MPI_Comm_rank(MPI_COMM_WORLD, &_thisRank);
	if (_thisRank == 0) {
		std::map<std::string, int> processorsRaplRank;
		int numberOfRanks;
		MPI_Comm_size(MPI_COMM_WORLD, &numberOfRanks);
		for (int otherRank = 1; otherRank < numberOfRanks; otherRank++) {
			char otherRankProcessorName[MPI_MAX_PROCESSOR_NAME];
			MPI_Recv(otherRankProcessorName, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, otherRank, /* tag */ otherRank, MPI_COMM_WORLD,
					 MPI_STATUS_IGNORE);
			processorsRaplRank[otherRankProcessorName] = otherRank;
		}
		processorsRaplRank[_processorName] = _thisRank;  // Root rank should measure
		std::set<int> raplRanks;                        // Only measure on one rank per processor (host)
		for (auto item : processorsRaplRank) {
			Log::global_log->debug() << "[" << getPluginName() << "] RAPL measurements for " << item.first
									 << " are performed on MPI rank " << item.second << std::endl;
			raplRanks.insert(item.second);
		}
		_thisRankShouldMeasure = raplRanks.find(_thisRank) != raplRanks.end();
		for (int i = 1; i < numberOfRanks; i++) {
			bool otherRankShouldMeasure = raplRanks.find(i) != raplRanks.end();
			MPI_Send(&otherRankShouldMeasure, 1, MPI_CXX_BOOL, /* dest */ i, /* tag */ i, MPI_COMM_WORLD);
		}
	} else {
		MPI_Send(_processorName, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, /* dest */ 0, /* tag */ _thisRank, MPI_COMM_WORLD);
		MPI_Recv(&_thisRankShouldMeasure, 1, MPI_CXX_BOOL, 0, /* tag */ _thisRank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	if (!_thisRankShouldMeasure) {
		return;  // Do not set up any counters
	}

#endif
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
	for (auto counter : _counters) {
		counter.reset();
	}
#ifdef ENABLE_MPI
	if (!_thisRankShouldMeasure || (!_per_host && _thisRank != 0)) {
		return;  // Only output on one rank per host/root rank
	}
#endif
	_simstart = std::chrono::steady_clock::now();
	if (!_outputprefix.empty()) {
		std::string outputFilename = getOutputFilename();
		std::ofstream outputFile(outputFilename.c_str());
		outputFile << "milliseconds\tsimstep\tjoule" << std::endl;
	}
}

void EnergyRAPL::endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
						 unsigned long simstep) {
	_joule = 0;
	for (auto counter : _counters) {
		_joule += counter.update();
	}
	_simstep = simstep;
	if (0 < _writeFrequency && simstep % _writeFrequency == 0) {
		outputEnergyJoule();
	}
}

void EnergyRAPL::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	xmlconfig.getNodeValue("outputprefix", _outputprefix);
	#ifdef ENABLE_MPI
	xmlconfig.getNodeValue("per-host", _per_host);
	#endif
}

void EnergyRAPL::outputEnergyJoule() {
	double joule = _joule;
#ifdef ENABLE_MPI
	if (!_per_host) {
		MPI_Reduce(&_joule, &joule, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		if (_thisRank != 0) {
			return;  // Only output on root rank
		}
	}
	else if (!_thisRankShouldMeasure) {
		return; // Only output on one rank per host
	}
#endif
	std::string jouleStr;
	std::stringstream sstream;
	sstream << std::fixed << std::setprecision(6) << joule;
	jouleStr = sstream.str();
	if (_outputprefix.empty()) {
	std::stringstream logLine;
#ifdef ENABLE_MPI
		if(_per_host) {
			logLine << "Host = " << _processorName << "\t";
		}
#endif
		 logLine << "Simstep = " << _simstep << "\tEnergy consumed = " << jouleStr << " J" << std::endl;
		 Log::global_log->info() << logLine.str() << std::flush;
	} else {
		std::string outputFilename = getOutputFilename();
		std::ofstream outputFile(outputFilename.c_str(), std::ios_base::app);
		int64_t milliseconds =
			std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - _simstart).count();
		outputFile << milliseconds << "\t" << _simstep << "\t" << jouleStr << std::endl;
	}
}

void EnergyRAPL::finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
	if (0 == _writeFrequency) {
		outputEnergyJoule();
	}
}
