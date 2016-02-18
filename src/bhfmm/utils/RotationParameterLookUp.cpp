/*
 * RotationParameterLookUp.cpp
 *
 *  Created on: Dec 15, 2014
 *      Author: uwe
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "RotationParameterLookUp.h"

namespace bhfmm {

RotationParameterLookUp* RotationParameterLookUp::tab = NULL;

inline double factorial(double n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

void RotationParameterLookUp::initFromFile(char* filename) {
	using namespace std;
	fstream file;

	file.open(filename, ios::in);
	unsigned count = 0;
	if (file.is_open()) {
		while (file >> table[count]) {
			++count;
			if (count >= totalNumEntries) {break;}
		}
	} else {
		cerr << "Failed to open file " << filename << endl;
		exit(EXIT_FAILURE);
	}

	file.close();
}

void RotationParameterLookUp::initFromRecursion()
{

	for(unsigned l = 0; l < numSlices; ++l) {
		/* set diagonal to 1.0 */
		for (unsigned m = 0; m <=l; ++m) {
			acc(l,m,m) = 1.0;
		}
		/* set initial value for recursion */
		if (l>0) {
			acc(l,l-1,l) = sqrt(2*l);
			acc(l,l,l-1) = 1.0/acc_c(l,l-1,l);

		}
		for(int m = l-2; m >= 0; --m) {
			acc(l,m,l) = sqrt(static_cast<double>(l+m+1)/static_cast<double>(l-m)) * acc_c(l, m+1, l);
			acc(l,l,m) = 1.0/acc_c(l,m,l);
			for (int k = l-1;  k >= m + 1 ; --k) {
				acc(l,m,k) = acc_c(l,m,l) / acc_c(l,k,l);
				acc(l,k,m) = 1.0/acc_c(l,m,k);
			}
		}
	}
}

void RotationParameterLookUp::initFromDirEval()
{
	for(unsigned l = 0; l < numSlices; ++l) {
		for(unsigned m = 0; m <= l; ++m) {
			for (unsigned k = 0; k <= l; ++k) {
//				this->acc(l,m,k) = sqrt(static_cast<double> (factorial(l-k)*factorial(l+k))/ static_cast<double>(factorial(l-m)*factorial(l+m)));
				this->acc(l,m,k) = sqrt((factorial(l-k)*factorial(l+k))/ (factorial(l-m)*factorial(l+m)));

			}
		}
	}
}

void RotationParameterLookUp::print(int maxl) {
	using namespace std;
	int precisionSetting = cout.precision( );
	ios::fmtflags flagSettings = cout.flags();

	cout.setf(ios::fixed | ios::showpos | ios::showpoint);
	cout.precision(5);
	for(int l = 0; l <= maxl; ++l) {
		for(int m = 0; m <= l; ++m) {
			for (int k = 0; k <= l; ++k) {
				cout << left << setw(10) << this->acc_c(l,m,k);
			}
			cout << endl;
		}
		cout << endl;
	}

	cout.precision(precisionSetting);
	cout.flags(flagSettings);
}


}  // namespace bhfmm
