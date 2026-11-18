/**
 * @file DomainDecompAnalyser.cpp
 * @author Georg von Bismarck
 * @date 18.11.2026
 */

#include "io/DomainDecompAnalyser.h"

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "Common.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"

#ifdef ENABLE_MPI
#include "parallel/DomainDecompMPIBase.h"
#endif


DomainDecompAnalyser::DomainDecompAnalyser() :
	_appendTimestamp(false), _outputPrefix("mardyn")
{}


void DomainDecompAnalyser::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	Log::global_log->info() << "DomainDecompAnalyser: Output prefix: " << _outputPrefix << std::endl;

	xmlconfig.getNodeValue("appendTimestamp", _appendTimestamp);
	Log::global_log->info() << "DomainDecompAnalyser: Append timestamp: " << _appendTimestamp << std::endl;
}


void DomainDecompAnalyser::init(ParticleContainer * /*particleContainer*/, DomainDecompBase * domainDecomp, Domain * domain) {
#ifdef ENABLE_MPI
	// store the initial DomainBox as reference for the comparisons in endStep()
	if (auto* mpiDecomp = dynamic_cast<DomainDecompMPIBase*>(domainDecomp)) {
		_rebuildStep.push_back(global_simulation->getSimulationStep());
		_rebuildDomainBox.push_back(mpiDecomp->getDomainBox(domain));
	}
#endif
}

void DomainDecompAnalyser::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain, unsigned long simstep) {
	_simstep.push_back(simstep);
	_traversalTime.push_back(domainDecomp->getLastTraversalTime());
	_averageTraversalTime.push_back(domainDecomp->getAverageLastTraversalTime());
	_numberParticles.push_back(particleContainer->getNumberOfParticles(ParticleIterator::ONLY_INNER_AND_BOUNDARY));

#ifdef ENABLE_MPI
	// store simstep and DomainBox whenever the DomainBox of this process has changed
	if (auto* mpiDecomp = dynamic_cast<DomainDecompMPIBase*>(domainDecomp)) {
		const DomainBox currentBox = mpiDecomp->getDomainBox(domain);
		if (_rebuildDomainBox.empty() || !_rebuildDomainBox.back().isNearEqual(currentBox)) {
			_rebuildStep.push_back(simstep);
			_rebuildDomainBox.push_back(currentBox);
		}
	}
#endif
}

template <typename T>
void DomainDecompAnalyser::writeVector(std::ostream& os, const std::string& name, const std::vector<T>& values) {
	os << name << ": [";
	for (std::size_t i = 0; i < values.size(); ++i) {
		if (i != 0) {
			os << ",";
		}
		os << values[i];
	}
	os << "]\n";
}

//This code is inspired by code from DecompWriter and DomainDecompMPIBase::printDecomp
void DomainDecompAnalyser::finish(ParticleContainer * /*particleContainer*/, DomainDecompBase * domainDecomp, Domain * /*domain*/) {
	std::stringstream filenamestream;
	filenamestream << _outputPrefix;
	if (_appendTimestamp) {
		filenamestream << "-" << gettimestring();
	}
	filenamestream << ".domainDecompAnalysis";
	const std::string filename = filenamestream.str();

	const int rank = domainDecomp->getRank();
	Log::global_log->info() << "DomainDecompAnalyser: Writing output to " << filename << std::endl;

	// header: general information and simsteps, written only by rank 0
	if (rank == 0) {
		// Placeholder for general information about the DomainDecompBase. It can be filled in the future
		// without changing the overall file structure (see the file structure rules in DomainDecompAnalyser.h).
		const std::string generalInfo = "";

		std::ofstream headerStream(filename);
		headerStream << "formatVersion: 1\n";
		headerStream << "[General information]\n";
		headerStream << generalInfo;
		headerStream << "[Data]\n";
		writeVector(headerStream, "simstep", _simstep);
		headerStream.close();
	}

	std::stringstream localInfo;
	localInfo << std::fixed << std::setprecision(std::numeric_limits<double>::digits10);
	localInfo << "[Process " << rank << "]\n";
	writeVector(localInfo, "traversalTime", _traversalTime);
	writeVector(localInfo, "averageTraversalTime", _averageTraversalTime);
	writeVector(localInfo, "numberParticles", _numberParticles);
	writeVector(localInfo, "rebuildStep", _rebuildStep);
#ifdef ENABLE_MPI
	localInfo << "rebuildDomainBox: " << _rebuildDomainBox << "\n";
#else
	localInfo << "rebuildDomainBox: []\n";
#endif
	const std::string localInfoStr = localInfo.str();

#ifdef ENABLE_MPI
	MPI_Comm comm = domainDecomp->getCommunicator();
	MPI_File fh;
	MPI_File_open(comm, filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_APPEND | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	uint64_t write_size = localInfoStr.size();
	uint64_t offset = 0;
	if (rank == 0) {
		MPI_Offset file_end_pos;
		MPI_File_seek(fh, 0, MPI_SEEK_END);
		MPI_File_get_position(fh, &file_end_pos);
		write_size += file_end_pos;
		MPI_Exscan(&write_size, &offset, 1, MPI_UINT64_T, MPI_SUM, comm);
		offset += file_end_pos;
	} else {
		MPI_Exscan(&write_size, &offset, 1, MPI_UINT64_T, MPI_SUM, comm);
	}
	MPI_File_write_at(fh, static_cast<MPI_Offset>(offset), localInfoStr.c_str(), static_cast<int>(localInfoStr.size()), MPI_CHAR, MPI_STATUS_IGNORE);
	MPI_File_close(&fh);
#else
	std::ofstream localStream(filename, std::ios::app);
	localStream << localInfoStr;
	localStream.close();
#endif
}
