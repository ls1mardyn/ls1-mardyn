#ifndef _TESTCASERUNNER_H_
#define _TESTCASERUNNER_H_

#include "utils/Log.h"


//! Run all Tests
class TestCaseRunner {
  public:
    //! Constructor
    TestCaseRunner();

    //! Destructor
    virtual ~TestCaseRunner();

    //! Run the tests
    void run();

  private:
    //! Logging interface
    static utils::Log _log;

    //! total number of errors for all Tests
    int _numberOfErrors;

    //! Test all Methods of the class Molecule  
    void runMoleculeTestCases();
    
    //! Test all Methods of the datastructure classes  
    void runDatastructureTestCases();
    
};

#endif
