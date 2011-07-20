/*
 * VTKGridWriter.cpp
 *
 * @Date: 24.08.2010
 * @Author: eckhardw
 */

#include "VTKGridWriter.h"
#include "particleContainer/LinkedCells.h"
#include "io/vtk/VTKGridWriterImplementation.h"
#include "particleContainer/Cell.h"
#include "Domain.h"
#include "utils/Logger.h"

using namespace Log;

VTKGridWriter::VTKGridWriter(unsigned int writeFrequency, const std::string& fileName, const LinkedCells& container)
: _writeFrequency(writeFrequency), _fileName(fileName), _container(container), _numCells(0), _numVertices(0) {

}

VTKGridWriter::~VTKGridWriter() { }


void  VTKGridWriter::doOutput(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			Domain* domain, unsigned long simstep,
			std::list<ChemicalPotential>* lmu
	) {

#ifndef NDEBUG
	if (static_cast<const void*>(&_container) != static_cast<const void*>(particleContainer)) {
		global_log->error() << "VTKGridWriter works only with PlottableLinkCells!" << std::endl;
		exit(1);
	}
#endif

	if (simstep % _writeFrequency != 0) {
		return;
	}

	int rank = domainDecomp->getRank();

	VTKGridWriterImplementation impl(rank);
	impl.initializeVTKFile();

	setupVTKGrid();

	for (int i = 0; i < _numCells; i++) {
		getCellData(_cells[i]);
		impl.plotCell(_cells[i]);
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

	releaseVTKGrid();

}

void  VTKGridWriter::outputParallelVTKFile(unsigned int numProcs, unsigned long simstep,
			VTKGridWriterImplementation& impl) {

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


void  VTKGridWriter::initOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomposition, Domain* domain) {
#ifndef NDEBUG
	if (static_cast<const void*>(&_container) != static_cast<const void*>(particleContainer)) {
		global_log->error() << "VTKGridWriter works only with LinkCells!" << std::endl;
		exit(1);
	}
#endif
}

void VTKGridWriter::setupVTKGrid() {

	int numCellsPerDimension[3];
	for (int i = 0; i < 3; i++) {
		numCellsPerDimension[i] = (_container._cellsPerDimension[i] - 2* _container._haloWidthInNumCells[i]);
	}

	_numCells = numCellsPerDimension[2] * numCellsPerDimension[1] * numCellsPerDimension[0];
	_cells = new VTKGridCell[_numCells];
	_numVertices = (numCellsPerDimension[2]+1) * (numCellsPerDimension[1]+1) * (numCellsPerDimension[0]+1);
	_vertices = new VTKGridVertex[_numVertices];

	// assign vertices their coordinates
	for (int i = 0; i < numCellsPerDimension[2]+1; i++) {
		for (int j = 0; j < numCellsPerDimension[1]+1; j++) {
			for (int k = 0; k < numCellsPerDimension[0]+1; k++) {
				int vertexIndex = i * (numCellsPerDimension[1]+1) * (numCellsPerDimension[1]+1) + j * (numCellsPerDimension[0]+1) + k;
				assert(vertexIndex < _numVertices);
				int x = k * _container._cellLength[0];
				int y = j * _container._cellLength[1];
				int z = i * _container._cellLength[2];
				_vertices[vertexIndex].setCoordinates(x, y, z);
			}
		}
	}

	// set cell index and vertices
	for (int i = 0; i < numCellsPerDimension[2]; i++) {
		for (int j = 0; j < numCellsPerDimension[1]; j++) {
			for (int k = 0; k < numCellsPerDimension[0]; k++) {
				int cellsIndex = i * numCellsPerDimension[1] * numCellsPerDimension[1] + j * numCellsPerDimension[0] + k;
				assert(cellsIndex < _numCells);
				// calculate the indices of the inner cells, taking the halo into account
				int containerIndex = _container.cellIndexOf3DIndex(k + _container._haloWidthInNumCells[0],
						j + _container._haloWidthInNumCells[1],
						i + _container._haloWidthInNumCells[2]);
				_cells[cellsIndex].setIndex(containerIndex);

				for (int l = 0; l < 8; l++) {
					int vertexIndex =   (i + (l & 4 ? 1 : 0)) * (numCellsPerDimension[1]+1) * (numCellsPerDimension[1]+1)
						                    		 + (j + (l & 2 ? 1 : 0)) * (numCellsPerDimension[0]+1) + (k + (l & 1));
					assert(vertexIndex < _numVertices);
					_cells[cellsIndex].setVertex(l, &_vertices[vertexIndex]);
				}

			}
		}
	}
}

void VTKGridWriter::releaseVTKGrid() {
	delete[] _vertices;
	delete[] _cells;
	_numCells = 0;
	_numVertices = 0;
}

void VTKGridWriter::getCellData(VTKGridCell& cell) {
	int numberOfMolecules = _container._cells[cell.getIndex()].getMoleculeCount();
	cell.setCellData(numberOfMolecules);
}

//! NOP
void  VTKGridWriter::finishOutput(ParticleContainer* particleContainer,
			DomainDecompBase* domainDecomp, Domain* domain) { }
