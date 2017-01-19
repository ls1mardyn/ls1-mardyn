#include "Random.h"

Random::Random() { this->init(8624); }

void Random::init(int seed)
{
   this->ix = (seed ^ (int)888889999) | (int)1;
   this->iy = seed ^ (int)777755555;
   // Calculate normalization factor
   this->am = 2.0 / (1.0 + (unsigned)((int)-1));
}

float Random::rnd()
{
   float rnd;
   const int IA = 16807;
   const int IM = 2147483647;
   const int IQ = 127773; 
   const int IR = 2836;

   /* Marsaglia shift sequence */
   ix ^= (ix >> 13);
   ix ^= (ix << 17);
   ix ^= (ix >> 5);

   int k = iy / IQ;
   iy = IA * (iy - k*IQ) - IR*k;
   if(iy < 0) iy += IM;
   rnd = am * ((IM & (ix ^ iy)) | (int)1);
   return rnd;
}
