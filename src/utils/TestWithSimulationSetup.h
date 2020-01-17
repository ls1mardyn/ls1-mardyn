/*
 * TestWithSimulationSetup.h
 *
 * @Date: 24.02.2011
 * @Author: eckhardw
 */

#ifndef TESTWITHSIMULATIONSETUP_H_
#define TESTWITHSIMULATIONSETUP_H_

#ifdef UNIT_TESTS

#include "utils/Testing.h"

namespace utils {
	class TestWithSimulationSetup;
}


class Domain;
class DomainDecompBase;
class ParticleContainer;


/**
 * A TestWithSimulationSetup has the ability to setup every thing needed for a
 * simulation, that is mainly
 *
 * - rank
 * - domain
 * - domainDecomposition
 * - particleContainer
 *
 * With those members set up it should be possible to test all the single aspects
 * of a molecular dynamics simulation.
 *
 * Those fields are setup from scratch before the execution of every test method, so
 * every test has a clean test bed.
 * The member objects are automatically destroyed after a test method has been executed.
 *
 * @note Testing the sequential algorithm independent of the number of MPI-processes
 *       Mardyn has been started with, should be possible if you replace the
 *       domainDecomposition with a dummyDomainDecomposition before the particleContainer
 *       is initialized.
 *
 * @note The components (i.e. molecule types) specified in the input file used
 *       in initializeFromFile() can be accessed via _domain->getComponents();
 */
class utils::TestWithSimulationSetup : public utils::Test {

public:

	TestWithSimulationSetup();

	virtual ~TestWithSimulationSetup();

	virtual void setUp();

	virtual void tearDown();


protected:
	/**
	 * Initialize a particle container from the given phase specification file.
	 * The domain class is setup as far as the input file is concerned, and for
	 * the initialization the _domain and _domainDecomposition of this test case
	 * are used.
	 *
	 * @param type the concrete implementation of the particleContainer to be created
	 * @param fileName the name of the file (relative to the testDataDirectory) which
	 *                 is used for initialization, or the prefix of the filename for binary files.
	 * @param cutoff the value used for all kinds of cutoff radii
	 * @param binary specifies that the file pointed to is a binary file. If true the header of this file has to have the
	 * extension ".header.xml" and the data file ".dat".
	 * @note The caller is responsible for deleting the particle container.
	 * @see ParticleContainerFactory::createInitializedParticleContainer()
	 */
	ParticleContainer* initializeFromFile(ParticleContainerFactory::Type type, const std::string& fileName,
										  double cutoff, bool binary = false);

	int _rank;

	Domain* _domain;

	DomainDecompBase* _domainDecomposition;


};

#endif /* UNIT_TESTS */

#endif /* TESTWITHSIMULATIONSETUP_H_ */
