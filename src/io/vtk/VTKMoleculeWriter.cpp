/*
 * VTKMoleculeWriter.cpp
 *
 * @Date: 24.08.2010
 * @Author: eckhardw
 */

#include "io/vtk/VTKMoleculeWriter.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "Domain.h"

#include <sstream>
#include <vector>
 #ifdef PARALLEL
#include <mpi.h>
#endif

VTKMoleculeWriter::VTKMoleculeWriter(unsigned int writeFrequency, const std::string& fileName)
: _writeFrequency(writeFrequency), _fileName(fileName) {
	if (_writeFrequency <= 0) {
		Log::global_log->error() << "VTKMoleculeWriter: writeFrequency must be > 0!" << std::endl;
	}
}

VTKMoleculeWriter::~VTKMoleculeWriter() {
}


void VTKMoleculeWriter::doOutput(
		ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		Domain* domain, unsigned long simstep,
		std::list<ChemicalPotential>* lmu
) {

	if (simstep % _writeFrequency != 0) {
		return;
	}

	int rank = domain->getlocalRank();

	VTKMoleculeWriterImplementation impl(rank);

	impl.initializeVTKFile();

	Molecule* tmpMolecule = particleContainer->begin();
	while (tmpMolecule != particleContainer->end()) {
		impl.plotMolecule(*tmpMolecule);
		tmpMolecule = particleContainer->next();
	}

	std::stringstream fileNameStream;
	fileNameStream << _fileName;

#ifdef PARALLEL
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
	for (int i = 0; i < numProcs; i++) {
		std::stringstream fileNameStream;
		fileNameStream << _fileName << "_node" << i << "_" << simstep << ".vtu";
		procFileNames.push_back(fileNameStream.str());
	}

	std::stringstream fileNameStream;
	fileNameStream << _fileName << "_" << simstep << ".pvtu";
	impl.initializeParallelVTKFile(procFileNames);
	impl.writeParallelVTKFile(fileNameStream.str());
}


//! NOP
void VTKMoleculeWriter::initOutput(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain) {}

//! NOP
void VTKMoleculeWriter::finishOutput(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, Domain* domain) {}
