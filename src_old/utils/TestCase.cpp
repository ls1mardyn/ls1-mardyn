#include "utils/TestCase.h"

utils::TestCase::TestCase( const std::string& testCaseName ):
  _testCaseName(testCaseName),
  _errors(0) {
  setAssertionOutputFormat;  	
}
    
utils::TestCase::~TestCase() {}
    
int utils::TestCase::getNumberOfErrors() const {
  return _errors;
}
