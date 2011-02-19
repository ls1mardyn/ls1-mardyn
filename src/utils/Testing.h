/*
 * Testing.h
 *
 * This file is used to encapsulate the testframework used (especially, I want to
 * be able to switch from CppUnit, which should be default, to the one used within
 * the peano project).
 *
 * @todo move the class Test to a seperate file!?
 *
 * @Date: 26.07.2010
 * @Author: Wolfgang Eckhardt
 *
 */

#ifndef TESTING_H_
#define TESTING_H_

#include "particleContainer/tests/ParticleContainerFactory.h"
#include <string>

class Domain;
class DomainDecompBase;
class ParticleContainer;

//! execute unit tests
//! @return false if no errors occured, true otherwise
bool runTests();

//! delegate to Test::setTestDataDirectory
void setTestDataDirectory(std::string& testDataDirectory);


#ifdef UNIT_TESTS

#define USE_CPPUNIT

#ifdef USE_CPPUNIT

#include <cppunit/extensions/HelperMacros.h>

namespace utils {
	typedef CppUnit::TestFixture TestBaseClass;
}


#define TEST_SUITE_REGISTRATION(name) CPPUNIT_TEST_SUITE_REGISTRATION(name)

#define TEST_SUITE(name) CPPUNIT_TEST_SUITE(name)
#define TEST_METHOD(name) CPPUNIT_TEST(name)
#define TEST_SUITE_END CPPUNIT_TEST_SUITE_END()

/***************  Assertion Macros  ****************/

/**
 * Assert that "expression" evaluates to true.
 * @see CPPUNIT_ASSERT
 */
#define ASSERT_TRUE(expression) CPPUNIT_ASSERT(expression)

/**
 * Assert that "expression" evaluates to true. If it is not, message msg is being output additionally to
 * the diagnostics of the assertion.
 * @see CPPUNIT_ASSERT_MESSAGE
 */
#define ASSERT_TRUE_MSG(msg, expression) CPPUNIT_ASSERT_MESSAGE(msg, expression)

#define ASSERT_FAIL(message) CPPUNIT_FAIL( message )

/**
 * Assert that a == b
 * @see CPPUNIT_ASSERT_EQUAL
 */
#define ASSERT_EQUAL(a,b) CPPUNIT_ASSERT_EQUAL(a,b)

/**
 * Assert that a == b. If it is not, message msg is being output additionally to
 * the diagnostics of the assertion.
 * @see CPPUNIT_ASSERT_EQUAL_MESSAGE
 */
#define ASSERT_EQUAL_MSG(msg, a, b) CPPUNIT_ASSERT_EQUAL_MESSAGE(msg, a, b)

/**
 * Assert that two double are equals given a tolerance.
 * @see CPPUNIT_ASSERT_DOUBLES_EQUAL
 */
#define ASSERT_DOUBLES_EQUAL(expected,actual,delta) CPPUNIT_ASSERT_DOUBLES_EQUAL(expected,actual,delta)

/**
 * Assert that two double are equals given a tolerance. If they are not,
 * message msg is being output additionally to the diagnostics of the assertion.
 * @see CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE
 */
#define ASSERT_DOUBLES_EQUAL_MSG(expected,actual,delta) CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(expected,actual,delta)

#else /***********  DEFINITIONS TO USE THE PEANO Test-Package    ****/



#include <tarch/tests/TestCase.h>
#include <tarch/tests/TestCaseFactory.h>

namespace utils {
	typedef tarch::tests::TestCase TestBaseClass;
}


#define TEST_SUITE_REGISTRATION(name) registerTest(name)

#define TEST_SUITE(name) void run() {
#define TEST_METHOD(name) testMethod(name)
#define TEST_SUITE_END() }

/***************  Assertion Macros  ****************/

/**
 * Assert that "expression" evaluates to true.
 * @see CPPUNIT_ASSERT
 */
#define ASSERT_TRUE(expression) validate(expression)

/**
 * Assert that "expression" evaluates to true. If it is not, message msg is being output additionally to
 * the diagnostics of the assertion.
 * @see CPPUNIT_ASSERT_MESSAGE
 */
#define ASSERT_TRUE_MSG(msg, expression) validateWithMessage(expression, msg)

#define ASSERT_FAIL(message) validateWithMessage(false, msg)

/**
 * Assert that a == b
 * @see CPPUNIT_ASSERT_EQUAL
 */
#define ASSERT_EQUAL(a,b) validateEquals(a, b)

/**
 * Assert that a == b. If it is not, message msg is being output additionally to
 * the diagnostics of the assertion.
 * @see CPPUNIT_ASSERT_EQUAL_MESSAGE
 */
#define ASSERT_EQUAL_MSG(msg, a, b) validateEqualsWithMessage(a, b, msg)

/**
 * Assert that two double are equals given a tolerance.
 * @see CPPUNIT_ASSERT_DOUBLES_EQUAL
 */
#define ASSERT_DOUBLES_EQUAL(expected,actual,delta) validateNumericalEqualsWithEps(actual, expected, delta)

/**
 * Assert that two double are equals given a tolerance. If they are not,
 * message msg is being output additionally to the diagnostics of the assertion.
 * @see CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE
 */
#define ASSERT_DOUBLES_EQUAL_MSG(expected,actual,delta) validateNumericalEqualsWithEps(actual, expected, delta)


#endif /* USE_CPPUNIT */


namespace utils {

/**
 * Now for every Test class it's rank, domain and domainDecomposition are setup
 * from scratch before the execution of every test method. The objects are also
 * automatically destroyed after the testcase has been executed.
 *
 * All the files which might be required by tests have to reside in one folder which
 * can be set via setTestDataDirectory().
 * If particle containers are created from input files, it is sufficient to just
 * give the name of the inputfile. The complete path is also determined by this
 * class.
 *
 * @note Testing the sequential algorithm indepentend of the number of MPI-processes
 *       Mardyn has been started with, should be possible if you replace the
 *       domainDecomposition with a dummyDomainDecomposition before the particleContainer
 *       is initialized.
 *
 * @todo Probably not every test needs the basic simulation classes like domain,
 *       domainDecomposition, particleContainer, etc... So should we move that into
 *       a seperate subclass?
 *
 */
class Test : public utils::TestBaseClass {

public:

	Test();
	~Test();


	virtual void setUp();

	virtual void tearDown();

	static void setTestDataDirectory(std::string& testDataDirectory);

protected:

	/**
	 * Initialize a particle container from the given phase specification file.
	 * The domain class is setup as far as the input file is concerned, and for
	 * the initialization the _domain and _domainDecomposition of this test case
	 * are used.
	 *
	 * @note The caller is responsible for deleting the particle container.
	 *
	 * @see ParticleContainerFactory::createInitializedParticleContainer()
	 */
	ParticleContainer* initializeFromFile(ParticleContainerFactory::type type, const char* fileName, double cutoff);

	int _rank;

	Domain* _domain;

	DomainDecompBase* _domainDecomposition;

private:

	std::string getTestDataFilename(const std::string& file);

	static std::string testDataDirectory;

};

}


#endif /* UNIT_TESTS */



#endif /* TESTING_H_ */
