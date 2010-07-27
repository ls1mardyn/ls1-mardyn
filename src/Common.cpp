// Common.cpp

#include "Common.h"

#include <sstream>
#include <ctime>
#include <cmath>

using namespace std;

string gettimestring(const char* fmt) {
	time_t rawtime;
	struct tm* timeinfo;
	char buffer[80];
	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, 80, fmt, timeinfo);
	string date(buffer);

	return date;
}

/**
 * Align a number to the right by padding with a filling character.
 *
 * @param number the number value to be aligned
 * @param num_digits the number of digits
 * @param c character used for filling up
 */
string aligned_number(int number, int num_digits, char c) {
	stringstream numstream;
	numstream.fill(c);
	numstream.width(num_digits);
	numstream << number;
	string numstr(numstream.str());
	return numstr;
}
