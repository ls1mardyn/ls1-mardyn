/*
 * VTKGridWriterTest.cpp
 *
 * @Date: 21.09.2010
 * @Author: eckhardw
 */

#ifndef MARDYN_AUTOPAS

#include "io/vtk/tests/VTKGridWriterTest.h"
#include "io/vtk/VTKGridWriter.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/LinkedCells.h"
#include "parallel/DomainDecompBase.h"
#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#endif
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
	VTKGridWriter writer(2, "VTKGridWriterTest");

#ifdef ENABLE_MPI
	// in the parallel case we check only that the right files are written.
	// Their content should be right, if the sequential tests pass.
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	DomainDecomposition domainDecomposition;
	Domain domain(rank);
	writer.endStep(container, &domainDecomposition, &domain, 1);
	writer.endStep(container, &domainDecomposition, &domain, 2);

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
	Domain domain(0);
	DomainDecompBase dummy;
    writer.endStep(container, &dummy, &domain, 1);
	ASSERT_TRUE_MSG("Check that files are written in the right interval.", !fileExists("VTKGridWriterTest_1.vtu"));

    writer.endStep(container, &dummy, &domain, 2);
	ASSERT_TRUE_MSG("Check that files are written in the right interval.", fileExists("VTKGridWriterTest_2.vtu"));

	removeFile("VTKGridWriterTest_2.vtu");
#endif
}


#endif
