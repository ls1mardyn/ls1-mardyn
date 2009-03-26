#include "utils/TestCaseCollection.h"

#include <iostream>

utils::Log utils::TestCaseCollection::_log("utils::TestCaseCollection");

utils::TestCaseCollection::TestCaseCollection(const std::string& testCaseCollectionName):
  TestCase::TestCase( testCaseCollectionName ) {}

utils::TestCaseCollection::~TestCaseCollection() {}

void utils::TestCaseCollection::run() {
  std::string logInformation = "running test case collection \"" + _testCaseName + "\" ";
  for (std::list<utils::TestCase*>::iterator p = _testCases.begin(); p!=_testCases.end(); p++ ) {
    utils::TestCase* currentTestCase = *p;
    currentTestCase->run();
    logInformation += ".";
    _errors += currentTestCase->getNumberOfErrors();
  }
  if (_errors==0) {
    logInformation += " ok";
  }
  else {
    logInformation += " failed";
  }
  _log.info("run()",logInformation );
}

void utils::TestCaseCollection::addTestCase( utils::TestCase* testCase ) {
  _testCases.push_back(testCase);
}
