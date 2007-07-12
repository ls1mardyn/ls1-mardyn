#ifndef MOLECULETESTS_H_
#define MOLECULETESTS_H_

#include "utils/TestCase.h"
#include "utils/Log.h"

/* namespace fluid { */
/*   class CalculatePressureGradientTest; */
/* } */

//! Test for Class Molecule
class MoleculeTest : public utils::TestCase {
  public:
    //! Constructor
    MoleculeTest();

    //! Destructor
    virtual ~MoleculeTest();
  
    //! Run the tests
    virtual void run();

  private:
    //! Logging interface
    static utils::Log _log;

    //! Test the constructor,  all get and all set methods
    void get_and_set_Test();
    
    //! Test the distance calculation
    void distance_calc_Test();

    //! Test the force calculation
    void force_calc_Test();

};

#endif /*MOLECULETESTS_H_*/
