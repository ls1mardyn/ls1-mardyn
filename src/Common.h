// Common.h

#ifndef COMMON_H_
#define COMMON_H_

#include <string>

class Common;

//                                        "%Y-%m-%d_%H-%M-%S"
std::string gettimestring(const char* fmt="%y%m%dT%H%M%S");

#endif /*COMMON_H_*/
