/*
 * FFTOrderReduction.cpp
 *
 *  Created on: Apr 18, 2016
 *  Author: gallardjm
 */

#include "FFTOrderReduction.h"



int FFTOrderReduction::_ReducedOrder[7];
int FFTOrderReduction::_order = 0;
  
void FFTOrderReduction::initialize(int order) {
  FFTOrderReduction::_order = order;
  
  int p = order+1;
  FFTOrderReduction::_ReducedOrder[0] = p;
  FFTOrderReduction::_ReducedOrder[1] = (p*3)/4;
  FFTOrderReduction::_ReducedOrder[2] = (p*5)/8;
  FFTOrderReduction::_ReducedOrder[3] = (p*9)/16;
  for(int i=4; i<7;i++)
        FFTOrderReduction::_ReducedOrder[i] = p/2;

/*  
  FFTOrderReduction::_ReducedOrder[0] = order+1;
  
  double r3 = sqrt(3);
  double r3_2 = 2*r3;
  
  //double coeff = 3.0 / pow(2, order+1);
  double coeff = pow(r3/2, order+1)/5;
  double d;
  
  for(int i = 1; i<7; ++i) {
    d = 4+i;
    //FFTOrderReduction::_ReducedOrder[i] = (int)ceil(log2((d-2)/(d+2)*coeff)/log2(2/d));
    FFTOrderReduction::_ReducedOrder[i] = (int)ceil(log2((d-r3_2)/(d+r3_2)*coeff)/log2(r3_2/d));
  }
*/  
  //for(int i=0;i<7;++i)
    //printf("ORED %i : %i \n",i+4,FFTOrderReduction::_ReducedOrder[i]);
}
  
int FFTOrderReduction::getM2LOrder(int x, int y, int z, int order) {
      
  if(FFTOrderReduction::_order != order)
    FFTOrderReduction::initialize(order);
   
  float dist2 = 4*(x*x + y*y + z*z);
  int dist_4 = floor(sqrt(dist2)) - 4;  
      
  return _ReducedOrder[dist_4];
      
}
 
