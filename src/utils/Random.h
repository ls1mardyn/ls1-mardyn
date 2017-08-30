#ifndef RANDOM_H
#define RANDOM_H

// by Stefan Becker
#include <cmath>

/** Uniform (quasi)-random number generator class.
 *
 * Minimal random number generator of Park and Miller combinmed with a Marsaglia shift sequence.
 * (see Numerical Recipies in Fortran 90: The Art of Parallel Scientific Computing, Chapter B7)
 */
class Random
{
 public:
   Random(int seed = 8624);
   void init(int seed);

   float rnd();

   int getIX() { return this->ix; }

   // by Stefan Becker
   /** returns a gaussian distributed deviate with zero mean and a standard deviation of stdDeviation
    the returned value is in the range +/- infinity (better: smallest, largest double number)*/
   double gaussDeviate(double stdDeviation);

 private:
   int ix, iy;
   float am;
};

#endif
