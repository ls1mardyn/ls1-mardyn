/**
 * @file DomainBoxTest.h
 * @author Georg von Bismarck
 * @date 17.11.2026
 */

#pragma once

#include "utils/Testing.h"

class DomainBox;

/**
 * Unit tests for DomainBox.
 */
class DomainBoxTest : public utils::Test {

	TEST_SUITE(DomainBoxTest);
	TEST_METHOD(testConstructors);
	TEST_METHOD(testFromArrayToArray);
	TEST_METHOD(testFromFlatArray);
	TEST_METHOD(testCreateAsRegularGrid);
	TEST_METHOD(testDimensionsAndVolume);
	TEST_METHOD(testIsValid);
	TEST_METHOD(testIntersection);
	TEST_METHOD(testIsNearEqual);
	TEST_METHOD(testLatchToGridSize);
	TEST_SUITE_END();

public:
	DomainBoxTest() = default;

	~DomainBoxTest() override = default;

	/**
	 * Default constructor, constructor from two corner points and constructor from six raw values.
	 */
	void testConstructors();

	/**
	 * fromArray() and toArray()
	 */
	void testFromArrayToArray();

	/**
	 * fromFlatArray(): The error case (size not a multiple of 6) is not tested, because MARDYN_EXIT terminates the process.
	 */
	void testFromFlatArray();

	/**
	 * createAsRegularGrid()
	 */
	void testCreateAsRegularGrid();

	/**
	 * width(), height(), depth() and volume()
	 */
	void testDimensionsAndVolume();

	/**
	 * isValid()
	 */
	void testIsValid();

	/**
	 * intersection() and intersectionVolume()
	 * intersectionVolume() is checked against intersection()->volume(); both have to be symmetric.
	 */
	void testIntersection();

	/**
	 * isNearEqual()
	 */
	void testIsNearEqual();

	/**
	 * latchToGridSize()
	 */
	void testLatchToGridSize();

private:
	/**
	 * Asserts that all six coordinates of two boxes are equal within delta.
	 * @param expected The expected box
	 * @param actual The box to check
	 * @param delta Allowed absolute deviation per coordinate
	 */
	void assertBoxEquals(const DomainBox& expected, const DomainBox& actual, double delta = 1e-12);
};
