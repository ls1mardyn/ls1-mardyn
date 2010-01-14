// Common.cpp

#include "Common.h"
#include <ctime>
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
