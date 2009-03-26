#include "utils/Globals.h"

int aPowI(int i,int a) {
  int result = 1;
  for (int d=0; d<i; d++) {
    result *= a;
  }
  return result;
}

int threePowI(int i) {
  return aPowI(i,3);
}

int fourPowI(int i) {
  return aPowI(i,4);
}

int twoPowI(int i) {
  return aPowI(i,2);
}
