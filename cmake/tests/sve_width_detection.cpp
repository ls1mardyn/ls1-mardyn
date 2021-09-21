// This test is here to determine the vectorization width of ARM SVE.

#include <arm_sve.h>
#include <iostream>

int main() { std::cout << svcntb() * 8; }