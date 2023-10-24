#include "Random.h"

Random::Random() { this->init(8624); }

void Random::init(int seed)
{
   this->ix_muVT = (seed ^ (int)888889999) | (int)1;
   this->iy_muVT = seed ^ (int)777755555;
   // Calculate normalization factor
   this->am_muVT = 2.0 / (1.0 + (unsigned)((int)-1));
}

float Random::rnd()
{
   float rnd;
   const int IA = 16807;
   const int IM = 2147483647;
   const int IQ = 127773;
   const int IR = 2836;

   ix_muVT ^= (ix_muVT >> 13);
   ix_muVT ^= (ix_muVT << 17);
   ix_muVT ^= (ix_muVT >> 5);

   int k = iy_muVT / IQ;
   iy_muVT = IA * (iy_muVT - k*IQ) - IR*k;
   if(iy_muVT < 0) iy_muVT += IM;
   rnd = am_muVT * ((IM & (ix_muVT ^ iy_muVT)) | (int)1);
   return rnd;
}
