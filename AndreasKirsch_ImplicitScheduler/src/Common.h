// Common.h

#ifndef COMMON_H_
#define COMMON_H_

#include <string>

class Common;

//                                        "%Y-%m-%d_%H-%M-%S"
std::string gettimestring(const char* fmt = "%y%m%dT%H%M%S");

//! @brief Helper returning an aligned number.
//!
std::string aligned_number(int number, int num_digits, char c = ' ');

#endif /*COMMON_H_*/
