/*
 * VTKGridWriterImplementationTest.h
 *
 * @Date: 20.09.2010
 * @Author: eckhardw
 */

#ifndef VTKGRIDWRITERIMPLEMENTATIONTEST_H_
#define VTKGRIDWRITERIMPLEMENTATIONTEST_H_

#include "utils/Testing.h"

class VTKGridWriterImplementationTest : public utils::Test {

	TEST_SUITE(VTKGridWriterImplementationTest);
	TEST_METHOD(testInitialization);
	TEST_METHOD(testWriteVTKFile);
	TEST_SUITE_END();

public:

	VTKGridWriterImplementationTest();

	virtual ~VTKGridWriterImplementationTest();

	void testInitialization();

	void testWriteVTKFile();
};

#endif /* VTKGRIDWRITERIMPLEMENTATIONTEST_H_ */
