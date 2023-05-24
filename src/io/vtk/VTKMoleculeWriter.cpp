/*
 * VTKMoleculeWriter.cpp
 *
 * @Date: 24.08.2010
 * @Author: eckhardw
 */

#include "io/vtk/VTKMoleculeWriter.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "utils/Logger.h"
#include "Domain.h"

#include <sstream>
#include <iostream>
#include <vector>


#ifdef ENABLE_MPI
#include <mpi.h>
#include <parallel/DomainDecompBase.h>
#endif


void VTKMoleculeWriter::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "VTKMoleculeWriter: Write frequency: " << _writeFrequency << std::endl;
	xmlconfig.getNodeValue("writefrequencyOffset", _writeFrequencyOffset);
	global_log->info() << "VTKMoleculeWriter: Write frequency offset: " << _writeFrequencyOffset << std::endl;
	// _writeInitialState is only relevant if the offset is != 0
	if (_writeFrequencyOffset != 0) {
		xmlconfig.getNodeValue("writeInitialState", _writeInitialState);
		global_log->info() << "VTKMoleculeWriter: Write initial state: " << _writeInitialState << std::endl;
	}
	xmlconfig.getNodeValue("outputprefix", _fileName);
	global_log->info() << "VTKMoleculeWriter: Output prefix: " << _fileName << std::endl;

	if (_writeFrequency <= 0) {
		Log::global_log->error() << "VTKMoleculeWriter: writeFrequency must be > 0!" << std::endl;
	}
}


void VTKMoleculeWriter::endStep(
        ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
        Domain * /*domain*/, unsigned long simstep) {
	// Do write only if any of these conditions is met:
	//   - We are at the beginning of the simulation and want to write the initial step.
	//   - We are beyond the offset and hit the frequency
	if (not ((simstep == 0 and _writeInitialState)
			or (simstep >= _writeFrequencyOffset and (simstep - _writeFrequencyOffset) % _writeFrequency == 0))) {
		return;
	}

	int rank = domainDecomp->getRank();

	VTKMoleculeWriterImplementation impl(rank, true);

	impl.initializeVTKFile();

	for (auto tmpMolecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
		 tmpMolecule.isValid(); ++tmpMolecule) {
		impl.plotMolecule(*tmpMolecule);
	}

	std::stringstream fileNameStream;
	fileNameStream << _fileName;

#ifdef ENABLE_MPI
	fileNameStream << "_node" << rank;

	if (rank == 0) {
		int numProcs = 0;
		MPI_CHECK( MPI_Comm_size(MPI_COMM_WORLD, &numProcs) );
		outputParallelVTKFile(numProcs,simstep, impl);
	}
#endif
	fileNameStream << "_" << simstep << ".vtu";
	impl.writeVTKFile(fileNameStream.str());
}


void VTKMoleculeWriter::outputParallelVTKFile(unsigned int numProcs, unsigned long simstep,
		VTKMoleculeWriterImplementation& impl) {

	std::vector<std::string> procFileNames;
	for (unsigned int i = 0; i < numProcs; i++) {
		std::stringstream fileNameStream;
		size_t pos = _fileName.find_last_of("/");
		if (pos != std::string::npos) {
			fileNameStream << _fileName.substr(pos+1) << "_node" << i << "_" << simstep << ".vtu";
		} else {
			fileNameStream << _fileName << "_node" << i << "_" << simstep << ".vtu";
		}
		procFileNames.push_back(fileNameStream.str());
	}

	std::stringstream fileNameStream;
	fileNameStream << _fileName << "_" << simstep << ".pvtu";
	impl.initializeParallelVTKFile(procFileNames);
	impl.writeParallelVTKFile(fileNameStream.str());
}


//! NOP
void VTKMoleculeWriter::init(ParticleContainer * /*particleContainer*/,
                             DomainDecompBase * /*domainDecomp*/, Domain * /*domain*/) {}

//! NOP
void VTKMoleculeWriter::finish(ParticleContainer * /*particleContainer*/,
							   DomainDecompBase * /*domainDecomp*/, Domain * /*domain*/) {}
