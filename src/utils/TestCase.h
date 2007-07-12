#ifndef _UTILS_TESTCASE_H_
#define _UTILS_TESTCASE_H_

#include <string>

#include "utils/Globals.h"

#define validate(booleanExpr,testCaseMethodName) if (!(booleanExpr)) { \
    _errors++; \
    std::cerr << "failure in test case " << _testCaseName << "::" << testCaseMethodName << ":" << std::endl \
              << "  boolean test failed " << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #booleanExpr << std::endl; \
  }

#define validateEquals(actualValue, validValue,testCaseMethodName) if (!(actualValue == validValue)) { \
    _errors++; \
    std::cerr << "failure in test case " << _testCaseName << "::" << testCaseMethodName << ":" << std::endl \
              << "  equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl; \
  }

#define validateEqualsWithParams1(actualValue, validValue,testCaseMethodName, param0) if (!(actualValue == validValue)) { \
    _errors++; \
    std::cerr << "failure in test case " << _testCaseName << "::" << testCaseMethodName << ":" << std::endl \
              << "  equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << std::endl; \
  }

#define validateEqualsWithParams2(actualValue, validValue,testCaseMethodName, param0, param1) if (!(actualValue == validValue)) { \
    _errors++; \
    std::cerr << "failure in test case " << _testCaseName << "::" << testCaseMethodName << ":" << std::endl \
              << "  equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << ", parameter " << #param1 << "=" << param1 << std::endl; \
  }


#define validateEqualsWithParams3(actualValue, validValue,testCaseMethodName, param0, param1, param2) if (!(actualValue == validValue)) { \
    _errors++; \
    std::cerr << "failure in test case " << _testCaseName << "::" << testCaseMethodName << ":" << std::endl \
              << "  equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << ", parameter " << #param1 << "=" << param1 << ", parameter " << #param2 << "=" << param2 << std::endl; \
  }

#define validateEqualsWithParams4(actualValue, validValue,testCaseMethodName, param0, param1, param2, param3) if (!(actualValue == validValue)) { \
    _errors++; \
    std::cerr << "failure in test case " << _testCaseName << "::" << testCaseMethodName << ":" << std::endl \
              << "  equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << ", parameter " << #param1 << "=" << param1 << ", parameter " << #param2 << "=" << param2 << ", parameter " << #param3 << "=" << param3 << std::endl; \
  }


#define validateNotEqual(actualValue, validValue,testCaseMethodName) if ((actualValue == validValue)) { \
    _errors++; \
    std::cerr << "failure in test case " << _testCaseName << "::" << testCaseMethodName << ":" << std::endl \
              << "  inequality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "!="  << #validValue << std::endl; \
  }

#define validateNumericalEquals(actualValue, validValue,testCaseMethodName) if (!equals(actualValue,validValue)) { \
    _errors++; \
    std::cerr << "failure in test case " << _testCaseName << "::" << testCaseMethodName << ":" << std::endl \
              << "  numerical equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  difference:  " << fabs((actualValue)-(validValue)) << std::endl; \
  }
  
#define validateNumericalEqualsWithEps(actualValue, validValue, testCaseMethodName, eps) \
  if (  (fabs((actualValue) - (validValue)) > eps ) ) { \
    _errors++; \
    std::cerr << "failure in test case " << _testCaseName << "::" << testCaseMethodName << ":" << std::endl \
              << "  numerical equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  difference:  " << fabs((actualValue)-(validValue)) << std::endl; \
  }  

#define validateNumericalEqualsWithParams1(actualValue, validValue,testCaseMethodName, param0) if (!equals(actualValue,validValue)) { \
    _errors++; \
    std::cerr << "failure in test case " << _testCaseName << "::" << testCaseMethodName << ":" << std::endl \
              << "  numerical equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << std::endl \
              << "  difference:  " << fabs((actualValue)-(validValue)) << std::endl; \
  }

#define validateNumericalEqualsWithParams2(actualValue, validValue,testCaseMethodName, param0, param1) if (!equals(actualValue,validValue)) { \
    _errors++; \
    std::cerr << "failure in test case " << _testCaseName << "::" << testCaseMethodName << ":" << std::endl \
              << "  numerical equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << ", parameter " << #param1 << "=" << param1 << std::endl \
              << "  difference:  " << fabs((actualValue)-(validValue)) << std::endl; \
  }

#define validateNumericalEqualsWithParams3(actualValue, validValue,testCaseMethodName, param0, param1, param2) if (!equals(actualValue,validValue)) { \
    _errors++; \
    std::cerr << "failure in test case " << _testCaseName << "::" << testCaseMethodName << ":" << std::endl \
              << "  numerical equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << ", parameter " << #param1 << "=" << param1 << ", parameter " << #param2 << "=" << param2 << std::endl \
              << "  difference:  " << fabs((actualValue)-(validValue)) << std::endl; \
  }

#define validateNumericalEqualsWithParams4(actualValue, validValue,testCaseMethodName, param0, param1, param2, param3) if (!equals(actualValue,validValue)) { \
    _errors++; \
    std::cerr << "failure in test case " << _testCaseName << "::" << testCaseMethodName << ":" << std::endl \
              << "  numerical equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << ", parameter " << #param1 << "=" << param1 << ", parameter " << #param2 << "=" << param2 << ", parameter " << #param3 << "=" << param3 << std::endl \
              << "  difference:  " << fabs((actualValue)-(validValue)) << std::endl; \
  }

#define validateNumericalEqualsWithParams5(actualValue, validValue,testCaseMethodName, param0, param1, param2, param3, param4) if (!equals(actualValue,validValue)) { \
    _errors++; \
    std::cerr << "failure in test case " << _testCaseName << "::" << testCaseMethodName << ":" << std::endl \
              << "  numerical equality test failed: " << actualValue << " instead of " << validValue << std::endl \
              << "  file: " << __FILE__ << " \t line: " << __LINE__ << std::endl \
              << "  statement: " << #actualValue << "=="  << #validValue << std::endl \
              << "  parameter " << #param0 << "=" << param0 << ", parameter " << #param1 << "=" << param1 << ", parameter " << #param2 << "=" << param2 << ", parameter " << #param3 << "=" << param3<< ", parameter " << #param4 << "=" << param4 << std::endl \
              << "  difference:  " << fabs((actualValue)-(validValue)) << std::endl; \
  }

namespace utils {
  class TestCase;
}

/**
 * Represents one test case. Every test case should be a subclass of this
 * class implementing run(). Furthermore subtypes should use the macros given
 * above to check any assumption. All the test cases are managed by the
 * TestCaseCollection.
 *
 * @see TestCaseCollection
 *
 * @author Tobias Weinzierl
 * @version $Revision: 1.16 $
 */
class utils::TestCase {
  protected:
    /**
     * Name of the test case.
     */
    std::string _testCaseName;

    /**
     * Error counter. This counter is increased on errors iff you use the
     * validate macros within your test cases.
     */
    int _errors;

  public:
    /**
     * Constructor.
     *
     * @param testCaseName Name of the test case. If a class is tested, this
     *                     string should equal the class name added the
     *                     namespace.
     */
    TestCase( const std::string& testCaseName );

    /**
     * Destructor.
     */
    virtual ~TestCase();

    /**
     * @return Number of errors.
     */
    int getNumberOfErrors() const;

    /**
     * This routine is triggered by the TestCaseCollection
     */
    void virtual run() = 0;
};

#endif
