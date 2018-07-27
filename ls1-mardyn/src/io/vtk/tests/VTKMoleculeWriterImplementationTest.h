/*
 * VTKMoleculeWriterImplementationTest.h
 *
 * @Date: 25.08.2010
 * @Author: eckhardw
 */

#ifndef VTKMOLECULEWRITERIMPLEMENTATIONTEST_H_
#define VTKMOLECULEWRITERIMPLEMENTATIONTEST_H_

#include "utils/Testing.h"

class VTKMoleculeWriterImplementationTest : public utils::Test {

	TEST_SUITE(VTKMoleculeWriterImplementationTest);
	TEST_METHOD(testInitialization);
	TEST_METHOD(testWriteVTKFile);
	TEST_SUITE_END();

public:

	VTKMoleculeWriterImplementationTest();

	virtual ~VTKMoleculeWriterImplementationTest();

	void testInitialization();

	void testWriteVTKFile();
};

#endif /* VTKMOLECULEWRITERIMPLEMENTATIONTEST_H_ */
