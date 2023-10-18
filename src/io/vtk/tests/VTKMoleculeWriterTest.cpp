/*
 * VTKMoleculeWriterTest.cpp
 *
 * @Date: 26.08.2010
 * @Author: eckhardw
 */

#include "VTKMoleculeWriterTest.h"
#ifndef MARDYN_AUTOPAS
#include "particleContainer/LinkedCells.h"
#else
#include "particleContainer/AutoPasContainer.h"
#endif
#include "parallel/DomainDecompBase.h"
#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#endif
#include "utils/FileUtils.h"
#include "utils/Logger.h"
#include "Domain.h"

#include "io/vtk/vtk-punstructured.h"

#include <sstream>


TEST_SUITE_REGISTRATION(VTKMoleculeWriterTest);

VTKMoleculeWriterTest::VTKMoleculeWriterTest() {

}

VTKMoleculeWriterTest::~VTKMoleculeWriterTest() {
}


void VTKMoleculeWriterTest::testDoOutput() {
	double boundings_min[] = {-1., -1., -1. };
	double boundings_max[] = {10.0, 10.0, 10.0 };
	double cutoff = 1.;
#ifndef MARDYN_AUTOPAS
	LinkedCells container(boundings_min, boundings_max, cutoff);
#else
	AutoPasContainer container(cutoff);
	container.rebuild(boundings_min, boundings_max);
#endif

	std::vector<Component> components;
	Component dummyComponent(0);
	dummyComponent.addLJcenter(0,0,0,0,0,0,0,false);
	components.push_back(dummyComponent);

	Molecule dummyMolecule1(0,&components[0],2.20,3.0,3.4,0,0,0,0,0,0,1,0,0,0);
	container.addParticle(dummyMolecule1);

	Molecule dummyMolecule2(0,&components[0], 1.0,0,0,0,0,0,0,0,0,1,0,0,0);
	container.addParticle(dummyMolecule2);

	Molecule dummyMolecule3(0,&components[0], 0,1.0,0,0,0,0,0,0,0,1,0,0,0);
	container.addParticle(dummyMolecule3);

	Molecule dummyMolecule4(0,&components[0], 1.0,1.5,0.234,0,0,0,0,0,0,1,0,0,0);
	container.addParticle(dummyMolecule4);

	container.updateMoleculeCaches();

	VTKMoleculeWriter writer(2, "VTKMoleculeWriterTest");

#ifdef ENABLE_MPI
	// in the parallel case we check only that the right files are written.
	// Their content should be right, if the sequential tests pass.
	int rank = 0;
	MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &rank) );
	DomainDecomposition domainDecomposition;
	Domain domain(rank);
	writer.endStep(&container, &domainDecomposition, &domain, 1);
	writer.endStep(&container, &domainDecomposition, &domain, 2);

	MPI_CHECK( MPI_Barrier(MPI_COMM_WORLD) );
	if (rank == 0) {
		ASSERT_TRUE_MSG("Check that files are written in the right interval.", !fileExists("VTKMoleculeWriterTest_1.pvtu"));
		ASSERT_TRUE_MSG("Check that files are written in the right interval.", fileExists("VTKMoleculeWriterTest_2.pvtu"));
		removeFile("VTKMoleculeWriterTest_2.pvtu");

		int numProcs = 0;
		MPI_CHECK( MPI_Comm_size(MPI_COMM_WORLD, &numProcs) );
		for (int i = 0; i < numProcs; i++) {
			std::stringstream str;
			str << "VTKMoleculeWriterTest_node" << i << "_2.vtu";
			std::stringstream msg;
			msg << "Check that parallel files are written in the right interval. (File " << str.str() <<")";
			ASSERT_TRUE_MSG(msg.str(), fileExists(str.str().c_str()));
			removeFile(str.str().c_str());
		}
	}
#else
	Domain domain(0);
	DomainDecompBase dummy;
    writer.endStep(&container, &dummy, &domain, 1);
	ASSERT_TRUE_MSG("Check that files are written in the right interval.", !fileExists("VTKMoleculeWriterTest_1.vtu"));

    writer.endStep(&container, &dummy, &domain, 2);
	ASSERT_TRUE_MSG("Check that files are written in the right interval.", fileExists("VTKMoleculeWriterTest_2.vtu"));

	try {
		std::unique_ptr<VTKFile_t> vtkFile(VTKFile ("VTKMoleculeWriterTest_2.vtu", xml_schema::flags::dont_validate));
		ASSERT_EQUAL( 4, (int) vtkFile->UnstructuredGrid()->Piece().NumberOfPoints());


		Points::DataArray_sequence& pointsArraySequence = vtkFile->UnstructuredGrid()->Piece().Points().DataArray();
		Points::DataArray_iterator coordinates_iterator = pointsArraySequence.begin();

		ASSERT_EQUAL_MSG("check number of point coordinates", 12, (int) coordinates_iterator->size());
		ASSERT_DOUBLES_EQUAL_MSG("check point coordinates", 2.2, coordinates_iterator->at(9),1e-7);
		ASSERT_DOUBLES_EQUAL_MSG("check point coordinates", 3., coordinates_iterator->at(10),1e-7);
		ASSERT_DOUBLES_EQUAL_MSG("check point coordinates", 3.4, coordinates_iterator->at(11),1e-7);
	} catch (const xml_schema::exception& e) {
		// todo printing the exception works only for ostream -> find some way to use our logger
		std::cerr << e << std::endl;
		ASSERT_FAIL("Could not parse .vtu-file!");
	}

	removeFile("VTKMoleculeWriterTest_2.vtu");
#endif
}
