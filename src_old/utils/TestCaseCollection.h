#ifndef _UTILS_TESTCASECOLLECTION_H_
#define _UTILS_TESTCASECOLLECTION_H_

#include "utils/TestCase.h"
#include "utils/Log.h"

#include <list>

namespace utils {
  class TestCaseCollection;
}


/**
 * Contains a sequence of tests which have to be executed sequentially. The
 * object manages a sequence of references. The user is responsible for
 * creating and destroying the test case instances.
 *
 * @author Tobias Weinzierl
 * @version $Revision: 1.6 $
 */
class utils::TestCaseCollection: public utils::TestCase {
  private:
    /**
     * Sequence of test cases that are executes on a run() call. The class is
     * not responsible for destroying them.
     */
    std::list<utils::TestCase*> _testCases;

    /**
     * Log interface the class writes to.
     */
    static utils::Log _log;
  public:
    /**
     * Creates a test case collection.
     *
     * @param testCaseCollectionName Name of the test case collection.
     */
    TestCaseCollection(const std::string& testCaseCollectionName);

    /**
     * Destructor
     */
    virtual ~TestCaseCollection();

    /**
     * Runs all test cases assigned.
     */
    virtual void run();

    /**
     * Adds a new test case.
     */
    void addTestCase( utils::TestCase* testCase );
};


#endif
