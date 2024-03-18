/*
 * CollectiveCommunicationTest.cpp
 *
 *  Created on: May 3, 2017
 *      Author: seckler
 */

#include "CollectiveCommunicationTest.h"

#include "parallel/CollectiveCommBase.h"
#include "parallel/CollectiveCommunication.h"
#include "parallel/CollectiveCommunicationNonBlocking.h"

#include <sstream>
#include <cmath>

TEST_SUITE_REGISTRATION(CollectiveCommunicationTest);


CollectiveCommunicationTest::CollectiveCommunicationTest() :
		_rank(0) {
	MPI_Comm_size(MPI_COMM_WORLD, &_commSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
}

CollectiveCommunicationTest::~CollectiveCommunicationTest() {
}

void CollectiveCommunicationTest::testCollectiveCommBase() {

}

void CollectiveCommunicationTest::testCollectiveCommunication() {
	CollectiveCommunication collComm;

	testSingleIteration(collComm);

}

void CollectiveCommunicationTest::testCollectiveCommunicationNonBlocking() {
#if MPI_VERSION >= 3
	CollectiveCommunicationNonBlocking collComm;

	testSingleIteration(collComm);

	collComm.init(MPI_COMM_WORLD, 1, 1);
	collComm.appendDouble(1.);
	collComm.allreduceSumAllowPrevious();
	double val1 = collComm.getDouble();
	collComm.finalize();
	ASSERT_DOUBLES_EQUAL(1. * _commSize, val1, 1e-8);

	collComm.init(MPI_COMM_WORLD, 1, 1);
	collComm.appendDouble(2.);
	collComm.allreduceSumAllowPrevious();
	val1 = collComm.getDouble();
	collComm.finalize();
	ASSERT_DOUBLES_EQUAL(1. * _commSize, val1, 1e-8);

	collComm.init(MPI_COMM_WORLD, 1, 1);
	collComm.appendDouble(3.);
	collComm.allreduceSumAllowPrevious();
	val1 = collComm.getDouble();
	collComm.finalize();
	ASSERT_DOUBLES_EQUAL(2. * _commSize, val1, 1e-8);

	collComm.init(MPI_COMM_WORLD, 1, 1);
	collComm.appendDouble(4.);
	collComm.allreduceSumAllowPrevious();
	val1 = collComm.getDouble();
	collComm.finalize();
	ASSERT_DOUBLES_EQUAL(3. * _commSize, val1, 1e-8);

	collComm.init(MPI_COMM_WORLD, 1, 2);
	collComm.appendDouble(-1.);
	collComm.allreduceSumAllowPrevious();
	val1 = collComm.getDouble();
	collComm.finalize();
	ASSERT_DOUBLES_EQUAL(-1. * _commSize, val1, 1e-8);

	collComm.init(MPI_COMM_WORLD, 1, 3);
	collComm.appendDouble(-1.);
	collComm.allreduceSumAllowPrevious();
	val1 = collComm.getDouble();
	collComm.finalize();
	ASSERT_DOUBLES_EQUAL(-1. * _commSize, val1, 1e-8);

	collComm.init(MPI_COMM_WORLD, 1, 2);
	collComm.appendDouble(-2.);
	collComm.allreduceSumAllowPrevious();
	val1 = collComm.getDouble();
	collComm.finalize();
	ASSERT_DOUBLES_EQUAL(-1. * _commSize, val1, 1e-8);

	collComm.init(MPI_COMM_WORLD, 1, 1);
	collComm.appendDouble(5.);
	collComm.allreduceSumAllowPrevious();
	val1 = collComm.getDouble();
	collComm.finalize();
	ASSERT_DOUBLES_EQUAL(4. * _commSize, val1, 1e-8);
#else
#pragma message "CollectiveCommunicationNonBlocking not supported and not tested. MPI_Version is too old."
#endif
}

void CollectiveCommunicationTest::testSingleIteration(CollectiveCommunicationInterface& collComm) {
	// allreduce with double
	collComm.init(MPI_COMM_WORLD, 1);
	collComm.appendDouble(2.);
	collComm.allreduceSum();
	double val1 = collComm.getDouble();
	collComm.finalize();
	ASSERT_DOUBLES_EQUAL(2. * _commSize, val1, 1e-8);

	// broadcast with int
	collComm.init(MPI_COMM_WORLD, 1);
	collComm.appendInt(_rank);
	collComm.broadcast(0);
	int val2 = collComm.getInt();
	collComm.finalize();
	mardyn_assert(val2 == 0);

	// scansum with unsigned long
	collComm.init(MPI_COMM_WORLD, 1);
	collComm.appendUnsLong(_rank);
	collComm.scanSum();
	unsigned long val3 = collComm.getUnsLong();
	collComm.finalize();
	unsigned long shouldBe3 = 0;
	for (int i = 0; i <= _rank; i++) {  // (we could indeed do this better)
		shouldBe3 += i;
	}
	mardyn_assert(val3 == shouldBe3);

	// allreduce with double, float, double
	collComm.init(MPI_COMM_WORLD, 3);
	collComm.appendDouble(1.);
	collComm.appendFloat(2.f);
	collComm.appendDouble(3.);
	collComm.allreduceSum();
	double val4 = collComm.getDouble();
	float val5 = collComm.getFloat();
	double val6 = collComm.getDouble();
	collComm.finalize();
	ASSERT_DOUBLES_EQUAL(1. * _commSize, val4, 1e-8);
	ASSERT_DOUBLES_EQUAL(2.f * _commSize, val5, 1e-8);
	ASSERT_DOUBLES_EQUAL(3. * _commSize, val6, 1e-8);
}
