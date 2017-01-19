#ifndef RANDOM_H
#define RANDOM_H


/** Uniform (quasi)-random number generator class.
 *
 * Minimal random number generator of Park and Miller combinmed with a Marsaglia shift sequence.
 * (see Numerical Recipies in Fortran 90: The Art of Parallel Scientific Computing, Chapter B7)
 */
class Random
{
 public:
   Random();
   void init(int seed);

   float rnd();

   int getIX() { return this->ix; }

 private:
   int ix, iy;
   float am;
};

#endif
