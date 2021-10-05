/*
 * Adios2IOTest.h
 *
 * Check whether an ADIOS2 checkpoint can be successfully written and read again.
 *
 *  Created on: July 2021
 *      Author: Heinen, Gralka, Rau, Homes
 */
#pragma once

#ifdef ENABLE_ADIOS2

#include "utils/TestWithSimulationSetup.h"
#include "molecules/Component.h"
#include <vector>
#include <array>

class Adios2IOTest : public utils::TestWithSimulationSetup {

	// declare testsuite inputFileTest
	TEST_SUITE(Adios2IOTest);

	// add a method which perform test
	TEST_METHOD(testWriteCheckpoint);

	// add a method which perform test
	TEST_METHOD(testReadCheckpoint);

	// end suite declaration
	TEST_SUITE_END();

public:
	Adios2IOTest() = default;
	virtual ~Adios2IOTest() = default;

	void initParticles();
	
	void testWriteCheckpoint();

	void testReadCheckpoint();

private:

	std::vector<unsigned long> _ids;
	std::vector<std::array<double,3>> _positions;
	std::vector<std::array<double,3>> _velocities;
	std::vector<std::array<double,3>> _Ds;
	std::vector<std::array<double,4>> _quaternions;
	std::array<double,3> _box_lower;
	std::array<double,3> _box_upper;

	std::vector<Component> _comp;

	std::shared_ptr<ParticleContainer> _inputPatricleContainer;
	std::shared_ptr<DomainDecompBase> _inputDomainDecomp;
	std::shared_ptr<Domain> _inputDomain;

	double _cutoff = 2.5;
	std::string _filename = "adios2restart.bp";
	
};

#endif