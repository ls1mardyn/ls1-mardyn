/*
 * DummyExpansion.cpp
 *
 *  Created on: Mar 31, 2016
 *      Author: gallardjm
 */

#include "DummyExpansion.h"

void DummyExpansion::evaluate_M_at_r(double X, double Y, double Z) {
	const int ord = order;

	const double r2 = X * X + Y * Y + Z * Z;
	const double inv_r2 = 1.0 / r2;
	const double inv_r = sqrt(inv_r2);

	// Equation 13.4.16 (a) and (b) Rappaport
	acc_C(0, 0) = inv_r;
	acc_S(0, 0) = 0.0;

	for (int l = 1; l <= ord; ++l) {
		int m;

		// 13.4.21
		for (m = 0; m <= l - 2; ++m) {
			const double factor1 = (2 * l - 1) * Z * inv_r2;
			const double factor2 = -(l - 1 + m) * (l - 1 - m) * inv_r2;

			acc_C(l, m) = factor1 * acc_c_C(l - 1, m)
					+ factor2 * acc_c_C(l - 2, m);
			acc_S(l, m) = factor1 * acc_c_S(l - 1, m)
					+ factor2 * acc_c_S(l - 2, m);
		}

		m = l - 1;
		// 13.4.21 when M(l-2,m) does not exist/is 0
		const double factor1 = (2 * l - 1) * Z * inv_r2;
		acc_C(l, m) = factor1 * acc_c_C(l - 1, m);
		acc_S(l, m) = factor1 * acc_c_S(l - 1, m);

		m = l;
		// 13.4.17, 13.4.18
		const double factor2 = -(2 * m - 1) * inv_r2;
		acc_C(m, m) = factor2
				* (X * acc_c_C(m - 1, m - 1) - Y * acc_c_S(m - 1, m - 1));
		acc_S(m, m) = factor2
				* (Y * acc_c_C(m - 1, m - 1) + X * acc_c_S(m - 1, m - 1));
	}
}
