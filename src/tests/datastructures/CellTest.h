#ifndef CELLTEST_H_
#define CELLTEST_H_

#include "utils/TestCase.h"
#include "utils/Log.h"

//! Test for Class Cell
class CellTest : public utils::TestCase {
  public:
    //! Constructor
    CellTest();

    //! Destructor
    virtual ~CellTest();
  
    //! Run the tests
    virtual void run();

  private:
    //! Logging interface
    static utils::Log _log;

    //! Test everything
    //! - get the list of Molecule pointers from a empty cell
    //! - remove all Molecules from a empty cell
    //! - insert some molecules into a cell
    //! - get the list of Molecule pointers from the filled cell
    //! - remove all molecules from the filled cell
    void allTests();
    

};

#endif /*CELLTEST_H_*/
