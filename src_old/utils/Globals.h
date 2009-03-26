/**
 * This file defines some macros and definitions used by most of the files of
 * the project. Among the definition of assertions and the global constants
 * (DIMENSIONS and NUMBER_OF_VERTICES_PER_ELEMENT) are the numerical comparison
 * operations.
 *
 * @version $Revision: 1.30 $
 * @author Tobias Weinzierl
 */
#ifndef _UTILS_GLOBALS_H_
#define _UTILS_GLOBALS_H_

  /**
   * The global header file provides a numerical equals operator that takes
   * the machine precision into account. Therefore the abs operator of the math
   * header is required.
   */
  #include <math.h>

  /**
   * Required for the inexact (numerical) comparisons.
   */
  #include <limits>

  /**
   * The IO stream operations are required for the assertion error messages.
   */
  #include <iostream>

  /**
   * Number of dimensions the code is working on. Furthermore the constant
   * NUMBER_OF_VERTICES_PER_ELEMENT is set. The value of this constant is
   * $f 2^d $f.
   */
  #ifdef Dim2
    #define DIMENSIONS 2
    #define TWO_POWER_D   2*2
    #define THREE_POWER_D 3*3
    #define FOUR_POWER_D  4*4
    #define NUMBER_OF_VERTICES_PER_ELEMENT TWO_POWER_D
  #elif Dim3
    #define DIMENSIONS 3
    #define TWO_POWER_D   2*2*2
    #define THREE_POWER_D 3*3*3
    #define FOUR_POWER_D  4*4*4
    #define NUMBER_OF_VERTICES_PER_ELEMENT TWO_POWER_D
  #elif Dim4
    #define DIMENSIONS 4
    #define TWO_POWER_D   2*2*2*2
    #define THREE_POWER_D 3*3*3*3
    #define FOUR_POWER_D  4*4*4*4
    #define NUMBER_OF_VERTICES_PER_ELEMENT TWO_POWER_D
  #endif


  #define NUMERICAL_ZERO_DIFFERENCE 1.0e-14
  #define equals(lhs,rhs) (fabs((rhs) - (lhs)) <= NUMERICAL_ZERO_DIFFERENCE )
  #define smaller(lhs,rhs) ((lhs - rhs) < -NUMERICAL_ZERO_DIFFERENCE )
  #define greater(lhs,rhs) ((lhs - rhs) >  NUMERICAL_ZERO_DIFFERENCE )



  /**
   * This vector type is a double vector of length DIMENSIONS.
   */
  //typedef blitz::TinyVector<double,DIMENSIONS> Vector;
  /**
   * This vector type is a integer vector of length DIMENSIONS.
   */
  //typedef blitz::TinyVector<int,DIMENSIONS> IntVector;

  /**
   * Define the assert macro. An assertion is given a boolean expression. If
   * the expression isn't true, the program immediatly quits giving the user
   * filename and line of the assertion failed. Assertions should be used to
   * verify preconditions and invariants but may not be used to validate any
   * arguments given by users. If the assertion fails, the program quits with
   * error code ASSERTION_EXIT_CODE.
   *
   * Whenever possible one should use the assertion macro assertMsg instead of
   * the pure assert. This operation is given an additional message
   * describing what the assertion does verify. This enables the programmer to
   * identify failures immediately.
   */
  #define ASSERTION_EXIT_CODE -1

  #define setAssertionOutputFormat { \
    std::cerr.setf( std::ios_base::scientific, std::ios_base::floatfield ); \
    std::cerr.precision(20); \
  }

  #ifdef Asserts
    #define assertion(expr) if (!(expr)) { \
      setAssertionOutputFormat; \
      std::cerr << "assertion in file " << __FILE__ << ", line " << __LINE__ << " failed: " << #expr <<  std::endl; \
      std::cerr.flush(); \
      exit(ASSERTION_EXIT_CODE); \
    }
    #define assertion1(expr,param) if (!(expr)) { \
      setAssertionOutputFormat; \
      std::cerr << "assertion in file " << __FILE__ << ", line " << __LINE__ << " failed: " << #expr <<  std::endl; \
      std::cerr << "parameter " << #param << ": " << param << std::endl; \
      std::cerr.flush(); \
      exit(ASSERTION_EXIT_CODE); \
    }
    #define assertionEquals(lhs,rhs) if ((lhs)!=(rhs)) { \
      setAssertionOutputFormat; \
      std::cerr << "assertion in file " << __FILE__ << ", line " << __LINE__ << " failed: " << #lhs << "==" #rhs << std::endl << #lhs << "=" << lhs << "\t" << #rhs << "=" << rhs << std::endl; \
      std::cerr.flush(); \
      exit(ASSERTION_EXIT_CODE); \
    }
    #define assertionEquals2(lhs,rhs,larg,rarg) if ((lhs)!=(rhs)) { \
      setAssertionOutputFormat; \
      std::cerr << "assertion in file " << __FILE__ << ", line " << __LINE__ << " failed: " << #lhs << "==" #rhs << std::endl << #lhs << "=" << lhs << "\t" << #rhs << "=" << rhs << std::endl; \
      std::cerr << "left argument " << #larg << ": " << larg << std::endl; \
      std::cerr << "right argument " << #rarg << ": " << rarg << std::endl; \
      std::cerr.flush(); \
      exit(ASSERTION_EXIT_CODE); \
    }
    #define assertionEquals3(lhs,rhs,larg,rarg,three) if ((lhs)!=(rhs)) { \
      setAssertionOutputFormat; \
      std::cerr << "assertion in file " << __FILE__ << ", line " << __LINE__ << " failed: " << #lhs << "==" #rhs << std::endl << #lhs << "=" << lhs << "\t" << #rhs << "=" << rhs << std::endl; \
      std::cerr << "first argument " << #larg << ": " << larg << std::endl; \
      std::cerr << "second argument " << #rarg << ": " << rarg << std::endl; \
      std::cerr << "third argument " << #three << ": " << three << std::endl; \
      std::cerr.flush(); \
      exit(ASSERTION_EXIT_CODE); \
    }
    #define assertionNumericalEquals(lhs,rhs) if (!equals( (lhs),(rhs) )) { setAssertionOutputFormat; std::cerr << "assertion in file " << __FILE__ << ", line " << __LINE__ << " failed: " << #lhs << "==" #rhs << std::endl << #lhs << "=" << lhs << "\t == \t" << #rhs << "=" << rhs << std::endl; exit(ASSERTION_EXIT_CODE); }
    #define assertionNumericalEquals2(lhs,rhs,larg,rarg) if (!equals( (lhs),(rhs) )) { \
      setAssertionOutputFormat; \
      std::cerr << "assertion in file " << __FILE__ << ", line " << __LINE__ << " failed: " << #lhs << "==" #rhs << std::endl << #lhs << "=" << lhs << "\t == \t" << #rhs << "=" << rhs << std::endl; \
      std::cerr << "left argument " << #larg << ": " << larg << std::endl; \
      std::cerr << "right argument " << #rarg << ": " << rarg << std::endl; \
      std::cerr.flush(); \
      exit(ASSERTION_EXIT_CODE); \
    }
    #define assertionNumericalEquals3(lhs,rhs,a,b,c) if (!equals( (lhs),(rhs) )) { \
      setAssertionOutputFormat; \
      std::cerr << "assertion in file " << __FILE__ << ", line " << __LINE__ << " failed: " << #lhs << "==" #rhs << std::endl << #lhs << "=" << lhs << "\t == \t" << #rhs << "=" << rhs << std::endl; \
      std::cerr << "first argument " << #a << ": " << a << std::endl; \
      std::cerr << "second argument " << #b << ": " << b << std::endl; \
      std::cerr << "third argument " << #c << ": " << c << std::endl; \
      std::cerr.flush(); \
      exit(ASSERTION_EXIT_CODE); \
    }
    #define assertionMsg(expr,message) if (!(expr)) { \
      setAssertionOutputFormat; std::cerr << "assertion in file " << __FILE__ << ", line " << __LINE__ << " failed: " << #expr << std::endl << message << std::endl ; \
      std::cerr.flush(); \
      exit(ASSERTION_EXIT_CODE); \
    }
    #define assertionFail(message) { \
      setAssertionOutputFormat; std::cerr << "fail-assertion in file " << __FILE__ << ", line " << __LINE__ <<  std::endl << message << std::endl; \
      std::cerr.flush(); \
      exit(ASSERTION_EXIT_CODE); \
    }
  #else
    #define assertion(expr)
    #define assertion1(expr,param)
    #define assertionEquals(lhs,rhs)
    #define assertionEquals2(lhs,rhs,a,b)
    #define assertionEquals3(lhs,rhs,a,b,c)
    #define assertionNumericalEquals(lhs,rhs)
    #define assertionNumericalEquals2(lhs,rhs,a,b)
    #define assertionNumericalEquals3(lhs,rhs,a,b,c)
    #define assertionMsg(expr,message)
    #define assertionFail(message)
  #endif
  
  
  int aPowI(int i,int a);
  int threePowI(int i);
  int fourPowI(int i);
  int twoPowI(int i);

#endif
