/*
 * WignerMatrix.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: uwe
 */

#include "utils/mardyn_assert.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "WignerMatrix.h"
#include "utils/Logger.h"


namespace bhfmm {

bool WignerMatrix::_sqrtFactorialLookUpInitialized = false;
std::vector<double> WignerMatrix::_sqrtFactorialLookUp;

inline int minus_one_pow(int m)
{
	/* -1^m */
	return 1 - ((m & 1) << 1);
}

void WignerMatrix::evaluate_0_Pip2(double theta)
{
	mardyn_assert(theta > 0 && theta <= M_PI/2);
	const int ord = order;

	/* g_lm (Dachsel Paper eq. (29))
	 *
	 * @Remark: it is possible to just use one double for g_l0 instead of an array
	 * for g_lm. In principle g_lm for m>0 could be temporarily stored in the matrix
	 * itself. However, for now,  using an array is easier to read and implement. */
	std::vector<double> g(ord+1); // ord+1 to have enough space for the last iteration l==order

	const double sint = sin(theta);
	const double cost = cos(theta);
	const double sint_per_1PLcost = sint/(1 + cost);

	/* calculate l=0 explicit to avoid branching */
	/* precompute g_lm (Dachsel Paper eq. (29)) */
	g[0] = 1.0; //
	acc(0,0,0) = 1.0; // g[0];

	for(int l = 1; l <= ord; ++l) {
		/* precompute g_lm (Dachsel Paper eq. (29)) */
		g[0] *= sqrt(double(2*l-1) / (2*l));
		for (int m = 1; m <= l; ++m) {
			g[m] = sqrt(double(l - m + 1) / (l + m))*g[m-1];
		}

		/* calculate last row (m=l)*/
		/* Dachsel Paper eq. (28) */
		acc(l,l,l) = minus_one_pow(2*l) * g[l]*pow(1 + cost, l); // l-m==0 !
		/* Dachsel Paper eq. (26) */

		for (int k = l; k > -l; --k) {
			const double factor = (l + k) / sqrt(double(l*(l+1)-k*(k-1))) * sint_per_1PLcost;
			acc(l,l,k-1) = factor*acc_c(l,l,k);
		}

		/* calculate remaining rows (m<l)*/
		for(int m = l-1; m >= 0; --m) {
			/* Dachsel Paper eq. (28) (k=l) */
			acc(l,m,l) = minus_one_pow(l+m) * g[m]*pow(1 + cost, m)*pow(sint,l-m);
			/* Dachsel Paper eq. (25) */
			for(int k = l; k > -l; --k) {
				const double factor1 = sqrt(double(l*(l+1) - m*(m+1))/(l*(l+1)-k*(k-1)));
				const double factor2 = (m+k)/sqrt(double(l*(l+1) - k*(k-1))) * sint_per_1PLcost;
				acc(l,m,k-1) = factor1*acc_c(l,m+1,k) + factor2*acc_c(l,m,k);
			} // for k
		} // for m
	} // for l
}

void WignerMatrix::evaluate(double theta)
{
	/* 2*PI periodicity, also for negative values */
	theta = std::fmod(theta, 2*M_PI);
	if (theta < 0) {
		theta = 2*M_PI + theta;
	}

	// ToDo: move cases for theta = 0,Pi, to methods that apply rotation to expansions
	// because rotation about these angles are special cases and can be performed more efficient.

	/* Dachsel Paper eq. (30) */
	if (theta == 0) {
		for(int l = 0; l <= int(order); ++l) {
			for(int m = 0; m <= l; ++m) {
				acc(l,m,m) = 1.0;
			}
		}

	} else if (theta > M_PI/2 && theta < M_PI) {
		evaluate_0_Pip2(M_PI-theta);
		mirror_k();
		apply_minus_one_pow<l_m>();

	} else if (theta == M_PI) {
		for(int l = 0; l <= int(order); ++l) {
			for(int m = 0; m <= l; ++m) {
				acc(l,m,-m) = 1.0;
			}
		}
		apply_minus_one_pow<l_k>();

	} else if (theta > M_PI && theta < 3*M_PI/2) {
		evaluate_0_Pip2(theta - M_PI);
		mirror_k();
		apply_minus_one_pow<l_k>();

	} else if (theta >= 3*M_PI/2 && theta < 2*M_PI){
		evaluate_0_Pip2(2*M_PI - theta);
		apply_minus_one_pow<m_k>();
	} else {
		evaluate_0_Pip2(theta);
	}

	scale();
}

void WignerMatrix::mirror_k() {
	for(int l = 0; l <= int(order); ++l) {
		for(int m = 0; m <= l; ++m) {
			for (int k = -l; k < 0; ++k) {
				std::swap(acc(l,m,k), acc(l,m,-k));
			}
		}
	}
}

template <class PowPair>
void WignerMatrix::apply_minus_one_pow() {
	for(int l = 0; l <= int(order); ++l) {
		for(int m = 0; m <= l; ++m) {
			for (int k = -l; k <= l; ++k) {
				acc(l,m,k) *= minus_one_pow(PowPair::sum(l,m,k));
			}

		}
	}
}


void WignerMatrix::print(int maxl) {
	int precisionSetting = std::cout.precision( );
	std::ios::fmtflags flagSettings = std::cout.flags();

	std::cout.setf(std::ios::fixed | std::ios::showpos | std::ios::showpoint);
	std::cout.precision(5);
	for(int l = 0; l <= maxl; ++l) {
		for(int m = 0; m <= l; ++m) {
			for (int k = -l; k <= l; ++k) {
				std::cout << std::left << std::setw(10) << this->acc_c(l,m,k);
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	std::cout.precision(precisionSetting);
	std::cout.flags(flagSettings);
}

void WignerMatrix::initializeSqrtLookUp() {

	unsigned numEntries = 2 * order + 1;
	_sqrtFactorialLookUp.resize(numEntries);

	// stores sqrt(i!)
	_sqrtFactorialLookUp[0] = 1.0;
	for(unsigned i = 1; i < numEntries; ++i) {
		_sqrtFactorialLookUp[i] = sqrt(static_cast<double>(i)) * _sqrtFactorialLookUp[i-1];
	}
}

void WignerMatrix::scale() {
	if(not _sqrtFactorialLookUpInitialized) {
		initializeSqrtLookUp();
	}

	if(_type == ROT_TYPE_L) {
		for (unsigned l = 0; l <= order; ++l) {
			for (unsigned m = 0; m <= l; ++m) {

				unsigned k = 0;
				{
					const double factor = lookUpFactor(l, m, k);
					acc(l, m, k) *= factor;
				}
				for (k = 1; k<=l; ++k) {
					const double factor = lookUpFactor(l, m, k);
					acc(l,m,-k) *= factor;
					acc(l,m,k) *= factor;
				}
			}
		}
	} else if (_type == ROT_TYPE_M) {
		//lookUpFactor(l, k, m) instead of lookUpFactor(l, m, k)
		for (unsigned l = 0; l <= order; ++l) {
			for (unsigned m = 0; m <= l; ++m) {

				unsigned k = 0; {
					const double factor = lookUpFactor(l, k, m);
					acc(l,m,k) *= factor;
				}
				for (k = 1; k<=l; ++k) {
					const double factor = lookUpFactor(l, k, m);
					acc(l,m,-k) *= factor;
					acc(l,m,k) *= factor;
				}
			}
		}
	} else {
		Log::global_log->error() << "type must be either L or M" << std::endl;
	}
}

}  // namespace bhfmm



