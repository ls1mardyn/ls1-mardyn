/*
 * VTKMoleculeWriterTest.h
 *
 * @Date: 26.08.2010
 * @Author: eckhardw
 */

#ifndef VTKMOLECULEWRITERTEST_H_
#define VTKMOLECULEWRITERTEST_H_

#include "utils/Testing.h"
#include "io/vtk/VTKMoleculeWriter.h"

class VTKMoleculeWriterTest : public utils::Test {

	TEST_SUITE(VTKMoleculeWriterTest);
	TEST_METHOD(testDoOutput);
	TEST_SUITE_END;

public:
	VTKMoleculeWriterTest();

	virtual ~VTKMoleculeWriterTest();

	void testDoOutput();
};

#endif /* VTKMOLECULEWRITERTEST_H_ */
