/*
 * VTKGridWriterImplementationTest.cpp
 *
 * @Date: 20.09.2010
 * @Author: eckhardw
 */

#include "io/vtk/tests/VTKGridWriterImplementationTest.h"
#include "io/vtk/VTKGridWriterImplementation.h"
#include "io/vtk/VTKGridCell.h"
#include "io/vtk/VTKGridVertex.h"
#include "utils/FileUtils.h"

#include <vector>

#ifdef ENABLE_MPI
#include <parallel/MPI_TIMED/mpi_timed.h>
#endif

TEST_SUITE_REGISTRATION(VTKGridWriterImplementationTest);

VTKGridWriterImplementationTest::VTKGridWriterImplementationTest() {
}

VTKGridWriterImplementationTest::~VTKGridWriterImplementationTest() {
}

void VTKGridWriterImplementationTest::testInitialization() {
	VTKGridWriterImplementation writer(0);

	ASSERT_EQUAL(writer.isVTKFileInitialized(), false);
	writer.initializeVTKFile();
	ASSERT_EQUAL(writer.isVTKFileInitialized(), true);

	ASSERT_EQUAL(writer.isParallelVTKFileInitialized(), false);
	std::vector<std::string> v;
	writer.initializeParallelVTKFile(v);
	ASSERT_EQUAL(writer.isParallelVTKFileInitialized(), true);
}

void VTKGridWriterImplementationTest::testWriteVTKFile() {

#ifdef ENABLE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank != 0) {
		return;
	}
#endif
	VTKGridWriterImplementation writer(0);

	VTKGridCell cell_1;
	cell_1.setIndex(0);
	VTKGridCell cell_2;
	cell_2.setIndex(1);

	std::vector<VTKGridVertex> vertices;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 3; k++) {
				VTKGridVertex vertex;
				vertex.setCoordinates(1. * i, 1. * j, 1. * k);
				vertices.push_back(vertex);
			}
		}
	}

	cell_1.setVertex(0, &vertices[0]);
	cell_1.setVertex(1, &vertices[1]);
	cell_1.setVertex(2, &vertices[6]);
	cell_1.setVertex(3, &vertices[7]);
	cell_1.setVertex(4, &vertices[3]);
	cell_1.setVertex(5, &vertices[4]);
	cell_1.setVertex(6, &vertices[9]);
	cell_1.setVertex(7, &vertices[10]);
	cell_1.setCellData(14, -1.0, 0);

	cell_2.setVertex(0, &vertices[1]);
	cell_2.setVertex(1, &vertices[2]);
	cell_2.setVertex(2, &vertices[7]);
	cell_2.setVertex(3, &vertices[8]);
	cell_2.setVertex(4, &vertices[4]);
	cell_2.setVertex(5, &vertices[5]);
	cell_2.setVertex(6, &vertices[10]);
	cell_2.setVertex(7, &vertices[11]);
	cell_2.setCellData(4, -1.0, 0);

	std::vector<std::string> v;
	writer.initializeVTKFile();
	writer.initializeParallelVTKFile(v);

	writer.plotCell(cell_1);
	ASSERT_EQUAL(writer.getNumCellsPlotted(), 1u);
	ASSERT_EQUAL(writer.getNumVerticesPlotted(), 8u);

	writer.plotCell(cell_2);
	ASSERT_EQUAL(writer.getNumCellsPlotted(), 2u);
	ASSERT_EQUAL(writer.getNumVerticesPlotted(), 12u);

	writer.writeVTKFile("VTKGridWriterImplementation00.vtu");
	ASSERT_TRUE(fileExists("VTKGridWriterImplementation00.vtu"));
	writer.writeParallelVTKFile("VTKGridWriterImplementation00.pvtu");
	ASSERT_TRUE(fileExists("VTKGridWriterImplementation00.pvtu"));

	removeFile("VTKGridWriterImplementation00.vtu");
	removeFile("VTKGridWriterImplementation00.pvtu");
}
