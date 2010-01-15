// Common.cpp

#include "Common.h"
#include <sstream>
#include <ctime>
#include <cmath>
using namespace std;

string gettimestring(const char* fmt) {
   time_t rawtime;
   struct tm * timeinfo;
   char buffer [80];
   time ( &rawtime );
   timeinfo = localtime ( &rawtime );

   strftime (buffer,80,fmt,timeinfo);
   string date(buffer);

   return date;
}

std::string aligned_number( int number, int num_digits, char c ) {
  stringstream numstream;
  numstream.fill( c );
  numstream.width( num_digits );
  numstream << number;
  string numstr(numstream.str());
  return numstr;
}
