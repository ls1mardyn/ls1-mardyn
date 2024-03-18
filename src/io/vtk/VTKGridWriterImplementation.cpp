/*
 * VTKGridWriterImplementation.cpp
 *
 * @Date: 20.09.2010
 * @Author: eckhardw
 */

#include "io/vtk/VTKGridWriterImplementation.h"
#include "io/vtk/VTKGridCell.h"
#include "io/vtk/VTKGridVertex.h"
#include "io/vtk/vtk-punstructured.h"
#include "utils/Logger.h"

#include <fstream>


VTKGridWriterImplementation::VTKGridWriterImplementation(int rank, const std::vector<double>* processorSpeeds)
: _vtkFile(NULL), _parallelVTKFile(NULL), _numCellsPlotted(0),
  _numVerticesPlotted(0), _rank(rank), _processorSpeeds(processorSpeeds) {
}


VTKGridWriterImplementation::~VTKGridWriterImplementation() {
	if (isVTKFileInitialized()) {
		delete _vtkFile;
	}
	if (isParallelVTKFileInitialized()) {
		delete _parallelVTKFile;
	}
}

/** TODO: there's redundancy in adding the data members (the same has to be done
 * for the parallel VTK File. -> use array of structs?
 */
void VTKGridWriterImplementation::initializeVTKFile() {
	PointData pointData;
	// we don't need point data at all!?
	//DataArray_t position(type::Int32, "id", 0);
	//pointData.DataArray().push_back(position);

	CellData cellData;
	DataArray_t ranks(type::Int32, "rank", 1);
	cellData.DataArray().push_back(ranks);
	DataArray_t cells_load(type::Float32, "load", 1);
	cellData.DataArray().push_back(cells_load);
	DataArray_t cells_level(type::Int32, "level", 1);
	cellData.DataArray().push_back(cells_level);
	DataArray_t index(type::UInt32, "index", 1);
	cellData.DataArray().push_back(index);
	if (_processorSpeeds != nullptr && _processorSpeeds->size() != 0) {
		DataArray_t procSpeeds(type::Float32, "processorSpeeds", 1);
		cellData.DataArray().push_back(procSpeeds);
		DataArray_t relLoad(type::Float32, "relativeLoad", 1);
		cellData.DataArray().push_back(relLoad);
	}

	// 3 coordinates
	Points points;
	DataArray_t pointCoordinates(type::Float32, "points", 3);
	points.DataArray().push_back(pointCoordinates);

	Cells cells;
	DataArray_t cells_connectivity(type::Int32, "connectivity", 1);
	cells.DataArray().push_back(cells_connectivity);
	DataArray_t cells_offsets(type::Int32, "offsets", 1);
	cells.DataArray().push_back(cells_offsets);
	DataArray_t cells_type(type::Int32, "types", 1);
	cells.DataArray().push_back(cells_type);

	PieceUnstructuredGrid_t piece(pointData, cellData, points, cells, 0, 0);
	UnstructuredGrid_t unstructuredGrid(piece);
	_vtkFile = new VTKFile_t("UnstructuredGrid");
	_vtkFile->UnstructuredGrid(unstructuredGrid);
}


void VTKGridWriterImplementation::plotCell(VTKGridCell& cell) {
	Points::DataArray_sequence& pointsArraySequence = (*_vtkFile).UnstructuredGrid()->Piece().Points().DataArray();
	Points::DataArray_iterator points_iterator = pointsArraySequence.begin();

	VTKGridVertex* const * vertices = cell.getVertices();
	for (int i = 0; i < 8; i++) {
		if (vertices[i]->getIndex() < 0) {
			const double* coordinates = vertices[i]->getCoordinates();
			for (int j = 0; j < 3; j++) {
				double coord = coordinates[j];
				points_iterator->push_back(coord);
			}
			vertices[i]->setIndex(_numVerticesPlotted);
			_numVerticesPlotted++;
		}
	}


	Cells::DataArray_sequence& cellsArraySequence = (*_vtkFile).UnstructuredGrid()->Piece().Cells().DataArray();
	Cells::DataArray_iterator it2 = cellsArraySequence.begin();

	// connectivity
	for (int i = 0; i < 8; i++) {
		it2->push_back(vertices[i]->getIndex());
	}

	_numCellsPlotted++;
	it2++;
	// offsets;
	it2->push_back(_numCellsPlotted * 8);

	it2++;
	// types: 8 for 2d-square, 11 for cube
	it2->push_back(11);

	CellData::DataArray_sequence& cellDataArraySequence = (*_vtkFile).UnstructuredGrid()->Piece().CellData().DataArray();
	CellData::DataArray_iterator it3 = cellDataArraySequence.begin();

	it3->push_back(cell.getRank());
	it3++;
	it3->push_back(cell.getLoad());
	it3++;
	it3->push_back(cell.getLevel());
	it3++;
	it3->push_back(cell.getIndex());
	it3++;
	if (_processorSpeeds != nullptr && _processorSpeeds->size() != 0) {
		it3->push_back((*_processorSpeeds)[cell.getRank()]);
		it3++;
		it3->push_back(cell.getLoad()/(*_processorSpeeds)[cell.getRank()]);
		it3++;
	}
}


void  VTKGridWriterImplementation::writeVTKFile(const std::string& fileName) {
#ifndef NDEBUG
	if (!isVTKFileInitialized()) {
		Log::global_log->error() << "VTKMoleculeWriterImplementation::writeVTKFile(): vtkFile not initialized!" << std::endl;
		return;
	}
#endif

	(*_vtkFile).UnstructuredGrid()->Piece().NumberOfPoints(_numVerticesPlotted);
	(*_vtkFile).UnstructuredGrid()->Piece().NumberOfCells(_numCellsPlotted);
	std::ofstream file(fileName.c_str());
	VTKFile (file, *_vtkFile);
}


void VTKGridWriterImplementation::initializeParallelVTKFile(const std::vector<std::string>& fileNames) {
	// init parallel file
		PPointData p_pointData;
		DataArray_t p_moleculeId(type::Int32, "id", 0);

		PCellData p_cellData;
		DataArray_t p_numberOfMolecules(type::Int32, "numberOfMolecules", 1);
		p_cellData.PDataArray().push_back( p_numberOfMolecules);
		DataArray_t p_load(type::Float32, "load", 1);
		p_cellData.PDataArray().push_back( p_load);
		DataArray_t p_cells_level(type::Int32, "level", 1);
		p_cellData.PDataArray().push_back(p_cells_level);
		DataArray_t p_node_rank(type::Int32, "node-rank", 1);
		p_cellData.PDataArray().push_back(p_node_rank);
		DataArray_t index(type::UInt32, "index", 1);
		p_cellData.PDataArray().push_back(index);

		// 3 coordinates
		PPoints p_points;
		DataArray_t p_pointCoordinates(type::Float32, "points", 3);
		p_points.PDataArray().push_back(p_pointCoordinates);

		PCells p_cells;
		DataArray_t p_cells_connectivity(type::Int32, "connectivity", 1);
		p_cells.PDataArray().push_back(p_cells_connectivity);
		DataArray_t p_cells_offsets(type::Int32, "offsets", 1);
		p_cells.PDataArray().push_back(p_cells_offsets);
		DataArray_t p_cells_type(type::Int32, "types", 1);
		p_cells.PDataArray().push_back(p_cells_type);


		PUnstructuredGrid_t p_unstructuredGrid(p_pointData, p_cellData, p_points, p_cells);
		for (unsigned int i = 0; i < fileNames.size(); i++) {
			Piece p_piece(fileNames[i]);
			p_unstructuredGrid.Piece().push_back(p_piece);
		}

		_parallelVTKFile = new VTKFile_t("PUnstructuredGrid");
		_parallelVTKFile->PUnstructuredGrid(p_unstructuredGrid);
}


void VTKGridWriterImplementation::writeParallelVTKFile(const std::string& fileName) {
#ifndef NDEBUG
	if (!isParallelVTKFileInitialized()) {
		Log::global_log->error() << "VTKMoleculeWriterImplementation::writeParallelVTKFile(): parallelVTKFile not initialized!" << std::endl;
		return;
	}
#endif
	std::ofstream file(fileName.c_str());
	VTKFile (file, *_parallelVTKFile);
}


bool VTKGridWriterImplementation::isVTKFileInitialized() {
	return _vtkFile != NULL;
}

bool VTKGridWriterImplementation::isParallelVTKFileInitialized() {
	return _parallelVTKFile != NULL;
}

unsigned int VTKGridWriterImplementation::getNumCellsPlotted() {
	return _numCellsPlotted;
}

unsigned int VTKGridWriterImplementation::getNumVerticesPlotted() {
	return _numVerticesPlotted;
}
