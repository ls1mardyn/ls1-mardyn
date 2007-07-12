#include "TestCaseRunner.h"
#include "utils/Watch.h"
#include "utils/TestCaseCollection.h"
#include "tests/MoleculeTest.h"
#include "tests/datastructures/CellTest.h"
#include <iostream>
using namespace std;


utils::Log TestCaseRunner::_log("TestCaseRunner");

TestCaseRunner::TestCaseRunner():
  _numberOfErrors(0) {
}

TestCaseRunner::~TestCaseRunner() {
}


void TestCaseRunner::runMoleculeTestCases() {
  utils::TestCaseCollection MoleculeTestCases( "Molecule" );
  MoleculeTest moleculeTest;

  MoleculeTestCases.addTestCase( &moleculeTest );

  MoleculeTestCases.run();

  _numberOfErrors += MoleculeTestCases.getNumberOfErrors();
}


void TestCaseRunner::runDatastructureTestCases() {
  utils::TestCaseCollection DatastructureTestCases( "Datastructures" );
  CellTest cellTest;

  DatastructureTestCases.addTestCase( &cellTest );

  DatastructureTestCases.run();

  _numberOfErrors += DatastructureTestCases.getNumberOfErrors();
}



void TestCaseRunner::run() {
  _log.info( "run(...)", "running test cases" );
  utils::Watch watch("TestCaseRunner", "run(...)");

  _numberOfErrors = 0;

  runMoleculeTestCases();
  runDatastructureTestCases();

  if ( _numberOfErrors == 0) {
    _log.info( "run(...)", "automatic tests successful" );
  }
  else {
    _log.error( "run(...)", "automatic tests failed" );
  }
}
