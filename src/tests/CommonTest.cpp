/*
 * CommonTest.cpp
 *
 * @Date: 21.05.2010
 * @Author: Wolfgang Eckhardt
 */

#include "CommonTest.h"
#include "Common.h"
#include "utils/Logger.h"

#include <iostream>

TEST_SUITE_REGISTRATION(CommonTest);

CommonTest::CommonTest() {
}

CommonTest::~CommonTest() {
}

void CommonTest::testGetTimeString() {
	std::string time = gettimestring();
	ASSERT_EQUAL((size_t)13, time.size());
	ASSERT_EQUAL(time[6], 'T');
}


void CommonTest::testAlignedNumber() {
//	ASSERT_EQUAL_MSG("One should be zero!", 1, 0);

	std::string result = aligned_number(123, 7, '.');
	std::string expected("....123");
	ASSERT_EQUAL_MSG("Align number 123 to 7 digits, filling char is .", expected, result);

	result = aligned_number(-2, 6, ' ');
	expected = "    -2";
	ASSERT_EQUAL_MSG("Align number -2 to 6 digits, filling char is ' '(space)", expected, result);
}

void CommonTest::testCalculateDistances() {
	typedef float fp_type;

	fp_type* valuesA[3];
	fp_type* valuesB[3];
	for (int i = 0; i < 3; i++) {
		valuesA[i] = new fp_type[3];
		valuesB[i] = new fp_type[1];
		valuesB[i][0] = 4.;
		for (int j = 0; j < 3; j++) {
			valuesA[i][j] = j+1;
		}
	}

	for (int i = 0; i < 3; i++) {
		std::cout << "A Points: " << valuesA[0][i] << "," << valuesA[1][i] << "," << valuesA[2][i] << std::endl;
	}
	for (int i = 0; i < 1; i++) {
			std::cout << "B Points: " << valuesB[0][i] << "," << valuesB[1][i] << "," << valuesB[2][i] << std::endl;
	}


	fp_type** distances;
	distances = new fp_type*[3];
	for (int i = 0; i < 3; i++) {
		distances[i] = new fp_type[1];
	}
	fp_type** distanceVectors[3];
	for (int i = 0; i < 3; i++) {
		distanceVectors[i] = new fp_type*[3];
		for (int j = 0; j < 4; j++) {
			distanceVectors[i][j] = new fp_type[1];
		}
	}

	calculateDistances(valuesA, valuesB, 3, 1, distances, distanceVectors);

	ASSERT_DOUBLES_EQUAL(distanceVectors[0][0][0], -3, 0.000001);
	ASSERT_DOUBLES_EQUAL(distanceVectors[1][0][0], -3, 0.000001);
	ASSERT_DOUBLES_EQUAL(distanceVectors[2][0][0], -3, 0.000001);
	ASSERT_DOUBLES_EQUAL(distances[0][0], 27, 0.000001);

	ASSERT_DOUBLES_EQUAL(distanceVectors[0][1][0], -2, 0.000001);
	ASSERT_DOUBLES_EQUAL(distanceVectors[1][1][0], -2, 0.000001);
	ASSERT_DOUBLES_EQUAL(distanceVectors[2][1][0], -2, 0.000001);
	ASSERT_DOUBLES_EQUAL(distances[1][0], 12, 0.000001);

	ASSERT_DOUBLES_EQUAL(distanceVectors[0][2][0], -1, 0.000001);
	ASSERT_DOUBLES_EQUAL(distanceVectors[1][2][0], -1, 0.000001);
	ASSERT_DOUBLES_EQUAL(distanceVectors[2][2][0], -1, 0.000001);
	ASSERT_DOUBLES_EQUAL(distances[2][0], 3, 0.000001);
}
