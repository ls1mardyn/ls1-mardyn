/*
 * VTKGridWriterTest.h
 *
 * @Date: 21.09.2010
 * @Author: eckhardw
 */

#ifndef VTKGRIDWRITERTEST_H_
#define VTKGRIDWRITERTEST_H_

#include "utils/Testing.h"

class VTKGridWriterTest : public utils::Test {

	TEST_SUITE(VTKGridWriterTest);
	TEST_METHOD(testEmptyGrid);
	TEST_SUITE_END;

public:

	VTKGridWriterTest();

	virtual ~VTKGridWriterTest();

	void testEmptyGrid();
};

#endif /* VTKGRIDWRITERTEST_H_ */
