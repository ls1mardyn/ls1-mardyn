/*
 * VTKGridWriter.cpp
 *
 * @Date: 24.08.2010
 * @Author: eckhardw
 */

#ifndef MARDYN_AUTOPAS

#include "VTKGridWriter.h"
#include "particleContainer/LinkedCells.h"
#include "io/vtk/VTKGridWriterImplementation.h"
#include "particleContainer/Cell.h"
#include "particleContainer/ParticleCell.h"
#include "Domain.h"
#include "utils/Logger.h"
#include "parallel/DomainDecompBase.h"


VTKGridWriter::VTKGridWriter()
: _numCells(0), _numVertices(0) {
}

VTKGridWriter::VTKGridWriter(unsigned int frequency, std::string name)
	: _writeFrequency(frequency), _fileName(name), _numCells(0), _numVertices(0) {
}

VTKGridWriter::~VTKGridWriter() { }

void VTKGridWriter::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	Log::global_log->info() << "VTKMoleculeWriter: Write frequency: " << _writeFrequency << std::endl;
	xmlconfig.getNodeValue("outputprefix", _fileName);
	Log::global_log->info() << "VTKMoleculeWriter: Output prefix: " << _fileName << std::endl;

	if (_writeFrequency <= 0) {
		Log::global_log->error() << "VTKMoleculeWriter: writeFrequency must be > 0!" << std::endl;
	}
}


void  VTKGridWriter::endStep(
        ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
        Domain * /*domain*/, unsigned long simstep
) {

	LinkedCells* container = dynamic_cast<LinkedCells*>(particleContainer);
#ifndef NDEBUG
	if (container == NULL) {
		Log::global_log->error() << "VTKGridWriter works only with plottable LinkedCells!" << std::endl;
		MARDYN_EXIT(1);
	}
#endif

	if (simstep % _writeFrequency != 0) {
		return;
	}

	int rank = domainDecomp->getRank();

	VTKGridWriterImplementation impl(rank);
	impl.initializeVTKFile();

	setupVTKGrid(particleContainer);

	for (int i = 0; i < _numCells; i++) {
		getCellData(container, _cells[i]);
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


void  VTKGridWriter::init(ParticleContainer *particleContainer,
                          DomainDecompBase * /*domainDecomposition*/, Domain * /*domain*/) {
#ifndef NDEBUG
	if (dynamic_cast<LinkedCells*>(particleContainer) == NULL) {
		Log::global_log->error() << "VTKGridWriter works only with LinkCells!" << std::endl;
		MARDYN_EXIT(1);
	}
#endif
}

void VTKGridWriter::setupVTKGrid(ParticleContainer* particleContainer) {
	LinkedCells* lc = dynamic_cast<LinkedCells*>(particleContainer);

	int numCellsPerDimension[3];
	for (int i = 0; i < 3; i++) {
		numCellsPerDimension[i] = (lc->_cellsPerDimension[i] - 2* lc->_haloWidthInNumCells[i]);
	}

	_numCells = numCellsPerDimension[2] * numCellsPerDimension[1] * numCellsPerDimension[0];
	_cells = new VTKGridCell[_numCells];
	_numVertices = (numCellsPerDimension[2]+1) * (numCellsPerDimension[1]+1) * (numCellsPerDimension[0]+1);
	_vertices = new VTKGridVertex[_numVertices];

	// assign vertices their coordinates
	for (int i = 0; i < numCellsPerDimension[2]+1; i++) {
		for (int j = 0; j < numCellsPerDimension[1]+1; j++) {
			for (int k = 0; k < numCellsPerDimension[0]+1; k++) {
				int vertexIndex = i * (numCellsPerDimension[1]+1) * (numCellsPerDimension[0]+1) + j * (numCellsPerDimension[0]+1) + k;

				mardyn_assert(vertexIndex < _numVertices);
				int x = k * lc->_cellLength[0];
				int y = j * lc->_cellLength[1];
				int z = i * lc->_cellLength[2];
				_vertices[vertexIndex].setCoordinates(x, y, z);
			}
		}
	}

	// set cell index and vertices
	for (int i = 0; i < numCellsPerDimension[2]; i++) {
		for (int j = 0; j < numCellsPerDimension[1]; j++) {
			for (int k = 0; k < numCellsPerDimension[0]; k++) {
				int cellsIndex = i * numCellsPerDimension[1] * numCellsPerDimension[0] + j * numCellsPerDimension[0] + k;
				mardyn_assert(cellsIndex < _numCells);
				// calculate the indices of the inner cells, taking the halo into account
				int containerIndex = lc->cellIndexOf3DIndex(k + lc->_haloWidthInNumCells[0],
						j + lc->_haloWidthInNumCells[1],
						i + lc->_haloWidthInNumCells[2]);
				_cells[cellsIndex].setIndex(containerIndex);

				for (int l = 0; l < 8; l++) {
					int vertexIndex =   (i + (l & 4 ? 1 : 0)) * (numCellsPerDimension[1]+1) * (numCellsPerDimension[0]+1)
						                    		 + (j + (l & 2 ? 1 : 0)) * (numCellsPerDimension[0]+1) + (k + (l & 1));
					mardyn_assert(vertexIndex < _numVertices);
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

void VTKGridWriter::getCellData(LinkedCells* container, VTKGridCell& cell) {
	int numberOfMolecules = container->_cells[cell.getIndex()].getMoleculeCount();
	cell.setCellData(numberOfMolecules, 0.0, 0);
}

//! NOP
void  VTKGridWriter::finish(ParticleContainer * /*particleContainer*/,
							DomainDecompBase * /*domainDecomp*/, Domain * /*domain*/) { }

#endif /* MARDYN_AUTOPAS */
