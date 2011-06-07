/*
 * VTKGridWriterTest.cpp
 *
 * @Date: 21.09.2010
 * @Author: eckhardw
 */

#include "io/vtk/tests/VTKGridWriterTest.h"
#include "io/vtk/VTKGridWriter.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/LinkedCells.h"
#include "utils/FileUtils.h"
#include "particleContainer/tests/ParticleContainerFactory.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

typedef ParticleContainerFactory Factory;

TEST_SUITE_REGISTRATION(VTKGridWriterTest);

VTKGridWriterTest::VTKGridWriterTest() {
}

VTKGridWriterTest::~VTKGridWriterTest() {
}

void VTKGridWriterTest::testEmptyGrid() {
	ParticleContainer* container = Factory::createEmptyParticleContainer(Factory::LinkedCell);
	LinkedCells* linkedCells = dynamic_cast<LinkedCells*> (container);
	assert(linkedCells != NULL);

	VTKGridWriter writer(2, "VTKGridWriterTest", *linkedCells);

#ifdef ENABLE_MPI
	// in the parallel case we check only that the right files are written.
	// Their content should be right, if the sequential tests pass.
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	Domain domain(rank, NULL);
	writer.doOutput(container, NULL, &domain, 1, NULL);
	writer.doOutput(container, NULL, &domain, 2, NULL);

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
		ASSERT_TRUE_MSG("Check that files are written in the right interval.", !fileExists("VTKGridWriterTest_1.pvtu"));
		ASSERT_TRUE_MSG("Check that files are written in the right interval.", fileExists("VTKGridWriterTest_2.pvtu"));
		removeFile("VTKGridWriterTest_2.pvtu");

		int numProcs = 0;
		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
		for (int i = 0; i < numProcs; i++) {
			std::stringstream str;
			str << "VTKGridWriterTest_node" << i << "_2.vtu";
			std::stringstream msg;
			msg << "Check that parallel files are written in the right interval. (File " << str.str() <<")";
			ASSERT_TRUE_MSG(msg.str(), fileExists(str.str().c_str()));
			removeFile(str.str().c_str());
		}
	}
#else
	Domain domain(0, NULL);
	writer.doOutput(container, NULL, &domain, 1, NULL);
	ASSERT_TRUE_MSG("Check that files are written in the right interval.", !fileExists("VTKGridWriterTest_1.vtu"));

	writer.doOutput(container, NULL, &domain, 2, NULL);
	ASSERT_TRUE_MSG("Check that files are written in the right interval.", fileExists("VTKGridWriterTest_2.vtu"));

	removeFile("VTKGridWriterTest_2.vtu");
#endif


}
