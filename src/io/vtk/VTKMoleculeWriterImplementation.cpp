/*
 * VTKMoleculeWriterImplementation.cpp
 *
 * @Date: 25.08.2010
 * @Author: eckhardw
 */

#include "io/vtk/VTKMoleculeWriterImplementation.h"
#include "molecules/Molecule.h"
#include "io/vtk/vtk-punstructured.h"
#include "utils/Logger.h"

#include <fstream>


VTKMoleculeWriterImplementation::VTKMoleculeWriterImplementation(int rank, bool plotCenters)
: _vtkFile(NULL), _parallelVTKFile(NULL), _numMoleculesPlotted(0), _rank(rank), _plotCenters(plotCenters) {
}

VTKMoleculeWriterImplementation::~VTKMoleculeWriterImplementation() {
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
void VTKMoleculeWriterImplementation::initializeVTKFile() {

	// Add the data we want to output for each molecule.
	// The iterator over PointData traverses the DataArrays just in the order
	// in which we add them here.
	PointData pointData;
	DataArray_t moleculeId(type::UInt64, "id", 1);
	pointData.DataArray().push_back(moleculeId);
	DataArray_t componentId(type::Int32, "component-id", 1);
	pointData.DataArray().push_back(componentId);
	DataArray_t node_rank(type::Int32, "node-rank", 1);
	pointData.DataArray().push_back(node_rank);
	DataArray_t forces(type::Float32, "forces", 3);
	pointData.DataArray().push_back(forces);

	if (_plotCenters) {
		DataArray_t centerId(type::Float32, "center-id", 1);
		pointData.DataArray().push_back(centerId);
		DataArray_t centerType(type::UInt8, "center-type", 1);
		pointData.DataArray().push_back(centerType);
	}

	CellData cellData; // we don't have cell data => leave it empty

	// 3 coordinates
	Points points;
	DataArray_t pointCoordinates(type::Float32, "points", 3);
	points.DataArray().push_back(pointCoordinates);

	Cells cells; // we don't have cells, => leave it empty
	// for some reasons, we have to add a dummy entry for paraview
	DataArray_t cells_data(type::Float32, "types", 0);
	cells.DataArray().push_back(cells_data);

	PieceUnstructuredGrid_t piece(pointData, cellData, points, cells, 0, 0);
	UnstructuredGrid_t unstructuredGrid(piece);
	_vtkFile = new VTKFile_t("UnstructuredGrid");
	_vtkFile->UnstructuredGrid(unstructuredGrid);
}


void VTKMoleculeWriterImplementation::plotMolecule(Molecule& molecule) {

#ifndef NDEBUG
	if (!isVTKFileInitialized()) {
		global_log->error() << "VTKMoleculeWriterImplementation::plotMolecule(): vtkFile not initialized!" << std::endl;
		return;
	}
#endif

	if (_plotCenters) {
		int centerID = 0;
		for (size_t i = 0; i < molecule.numLJcenters(); i++) {
			plotCenter(molecule, centerID, LJ);
			centerID++;
		}
		for (size_t i = 0; i < molecule.numCharges(); i++) {
			plotCenter(molecule, centerID, Charge);
			centerID++;
		}
		for (size_t i = 0; i < molecule.numDipoles(); i++) {
			plotCenter(molecule, centerID, Dipole);
			centerID++;
		}
		for (size_t i = 0; i < molecule.numQuadrupoles(); i++) {
			plotCenter(molecule, centerID, Quadrupole);
			centerID++;
		}
	} else {
		PointData::DataArray_sequence& pointDataArraySequence =
				(*_vtkFile).UnstructuredGrid()->Piece().PointData().DataArray();
		PointData::DataArray_iterator data_iterator = pointDataArraySequence.begin();
		// id
		data_iterator->push_back(molecule.getID());
		data_iterator++;
		// componentID
		data_iterator->push_back(molecule.componentid());
		data_iterator++;
		// mpi-node rank
		data_iterator->push_back(_rank);
		data_iterator++;
		// forces
		data_iterator->push_back(molecule.F(0));
		data_iterator->push_back(molecule.F(1));
		data_iterator->push_back(molecule.F(2));

		// Coordinates
		Points::DataArray_sequence& pointsArraySequence = (*_vtkFile).UnstructuredGrid()->Piece().Points().DataArray();
		Points::DataArray_iterator coordinates_iterator = pointsArraySequence.begin();
		coordinates_iterator->push_back(molecule.r(0));
		coordinates_iterator->push_back(molecule.r(1));
		coordinates_iterator->push_back(molecule.r(2));
		_numMoleculesPlotted++;
	}
}


void VTKMoleculeWriterImplementation::plotCenter(Molecule& molecule, int centerID, CenterType centerType) {
	PointData::DataArray_sequence& pointDataArraySequence = (*_vtkFile).UnstructuredGrid()->Piece().PointData().DataArray();
	PointData::DataArray_iterator data_iterator = pointDataArraySequence.begin();

	data_iterator->push_back(molecule.getID());
	data_iterator++;

	data_iterator->push_back(molecule.componentid());
	data_iterator++;

	data_iterator->push_back(_rank);
	data_iterator++;

	const std::array<double, 3> center_force = molecule.site_F(centerID);
	data_iterator->push_back(center_force[0]);
	data_iterator->push_back(center_force[1]);
	data_iterator->push_back(center_force[2]);
	data_iterator++;

	data_iterator->push_back(centerID);
	data_iterator++;
	data_iterator->push_back(centerType);

	Points::DataArray_sequence& pointsArraySequence = (*_vtkFile).UnstructuredGrid()->Piece().Points().DataArray();
	Points::DataArray_iterator coordinates_iterator = pointsArraySequence.begin();

	const std::array<double, 3> curr_center = molecule.site_d_abs(centerID);
	coordinates_iterator->push_back(curr_center[0]);
	coordinates_iterator->push_back(curr_center[1]);
	coordinates_iterator->push_back(curr_center[2]);
	_numMoleculesPlotted++;
}

void  VTKMoleculeWriterImplementation::writeVTKFile(const std::string& fileName) {
#ifndef NDEBUG
	if (!isVTKFileInitialized()) {
		global_log->error() << "VTKMoleculeWriterImplementation::writeVTKFile(): vtkFile not initialized!" << std::endl;
		return;
	}
#endif

	(*_vtkFile).UnstructuredGrid()->Piece().NumberOfPoints(_numMoleculesPlotted); // sets the number of points
	std::ofstream file(fileName.c_str());
	VTKFile (file, *_vtkFile); //actually writes the file
}

void VTKMoleculeWriterImplementation::initializeParallelVTKFile(const std::vector<std::string>& fileNames) {
	// init parallel file
	PPointData p_pointData;
	DataArray_t p_moleculeId(type::Float32, "id", 1);
	p_pointData.PDataArray().push_back(p_moleculeId);
	DataArray_t p_componentId(type::Float32, "component-id", 1);
	p_pointData.PDataArray().push_back(p_componentId);
	DataArray_t p_node_rank(type::Int32, "node-rank", 1);
	p_pointData.PDataArray().push_back(p_node_rank);
	DataArray_t p_forces(type::Float32, "forces", 3);
	p_pointData.PDataArray().push_back(p_forces);

	if (_plotCenters) {
		DataArray_t p_centerId(type::Float32, "center-id", 1);
		p_pointData.PDataArray().push_back(p_centerId);
		DataArray_t p_centerType(type::UInt8, "center-type", 1);
		p_pointData.PDataArray().push_back(p_centerType);
	}

	PCellData p_cellData; // we don't have cell data => leave it empty

	// 3 coordinates
	PPoints p_points;
	DataArray_t p_pointCoordinates(type::Float32, "points", 3);
	p_points.PDataArray().push_back(p_pointCoordinates);

	PCells p_cells; // we don't have cells, => leave it empty
	// for some reasons, we have to add a dummy entry for paraview
	DataArray_t p_cells_data(type::Float32, "types", 0);
	p_cells.PDataArray().push_back(p_cells_data);

	PUnstructuredGrid_t p_unstructuredGrid(p_pointData, p_cellData, p_points, p_cells);
	for (size_t i = 0; i < fileNames.size(); i++) {
		Piece p_piece(fileNames[i]);
		p_unstructuredGrid.Piece().push_back(p_piece);
	}

	_parallelVTKFile = new VTKFile_t("PUnstructuredGrid");
	_parallelVTKFile->PUnstructuredGrid(p_unstructuredGrid);
}


void VTKMoleculeWriterImplementation::writeParallelVTKFile(const std::string& fileName) {
#ifndef NDEBUG
	if (!isParallelVTKFileInitialized()) {
		global_log->error() << "VTKMoleculeWriterImplementation::writeParallelVTKFile(): parallelVTKFile not initialized!" << std::endl;
		return;
	}
#endif
	std::ofstream file(fileName.c_str());
	VTKFile (file, *_parallelVTKFile);
}


bool VTKMoleculeWriterImplementation::isVTKFileInitialized() {
	return _vtkFile != NULL;
}

bool VTKMoleculeWriterImplementation::isParallelVTKFileInitialized() {
	return _parallelVTKFile != NULL;
}

unsigned int VTKMoleculeWriterImplementation::getNumMoleculesPlotted() {
	return _numMoleculesPlotted;
}
