/*
 * SolidHarmonicsExpansion.cpp
 *
 *  Created on: Nov 21, 2014
 *      Author: tchipevn, rajats
 */

#include "SolidHarmonicsExpansion.h"
#include <iomanip>
#include <algorithm>

namespace bhfmm {

// CONSTRUCTORS //
// note that order k requires k+1 rows!
SolidHarmonicsExpansion::SolidHarmonicsExpansion(int order, bool initializeToZero) :
		_order(order), _c(order + 1, initializeToZero), _s(order + 1, initializeToZero) {
}

SolidHarmonicsExpansion::SolidHarmonicsExpansion(const SolidHarmonicsExpansion & rhs) :
		_order(rhs._order), _c(rhs._c), _s(rhs._s) {
#ifdef FMM_FFT
	//copy the FFTData from FFTAccelerableExpansion inheritance if needed
	if (issetFFTData())
		_FFTData = rhs._FFTData->copyContainer();
#endif  /* FMM_FFT */
}

// DESTRUCTOR //
SolidHarmonicsExpansion::~SolidHarmonicsExpansion() {
}

void swap(SolidHarmonicsExpansion & s1, SolidHarmonicsExpansion & s2) {
	mardyn_assert(s1.getOrder() == s2.getOrder());
	swap(s1._c, s2._c);
	swap(s1._s, s2._s);
}

// OPERATORS //

SolidHarmonicsExpansion & SolidHarmonicsExpansion::operator=(SolidHarmonicsExpansion RHS) {
	mardyn_assert(this->_order == RHS._order);
	swap(*this, RHS);
	return *this;
}

SolidHarmonicsExpansion & SolidHarmonicsExpansion::operator+=(const SolidHarmonicsExpansion & RHS) {
	mardyn_assert(this->_order == RHS._order);
	this->_c += RHS._c;
	this->_s += RHS._s;
	return *this;
}

SolidHarmonicsExpansion operator+(SolidHarmonicsExpansion LHS, const SolidHarmonicsExpansion & RHS) {
	LHS += RHS;
	return LHS;
}

SolidHarmonicsExpansion & SolidHarmonicsExpansion::operator*=(double scalar) {
	this->_c *= scalar;
	this->_s *= scalar;
	return *this;
}

SolidHarmonicsExpansion operator*(double scalar, SolidHarmonicsExpansion RHS) {
	RHS *= scalar;
	return RHS;
}

// METHODS //
void SolidHarmonicsExpansion::clear() {
	_c.setToZero();
	_s.setToZero();
}

// math operators
void SolidHarmonicsExpansion::evaluateLOfR(Vector3<double> r) {
	const int ord = _order;

	const double r2 = r.L2NormSquare();
	const double X = r[0], Y = r[1], Z = r[2];

	// Equation 13.4.16 (c) and (d) Rappaport
	//NOTE: there is a typo in Rapaport: L^c_{00} should be 1
	acc_C(0, 0) = 1.0;
	acc_S(0, 0) = 0.0;

	for (int l = 1; l <= ord; ++l) {
		int m;

		// 13.4.22
		for (m = 0; m <= l - 2; ++m) {
			const int denom = (l + m) * (l - m);
			const double factor1 = (2 * l - 1) * Z / denom;
			const double factor2 = -r2 / denom;

			acc_C(l, m) = factor1 * acc_c_C(l - 1, m) + factor2 * acc_c_C(l - 2, m);
			acc_S(l, m) = factor1 * acc_c_S(l - 1, m) + factor2 * acc_c_S(l - 2, m);
		}

		m = l - 1;
		// 13.4.22 when L(l-2,m) does not is 0
		const double factor1 = (2 * l - 1) * Z / ((l + m) * (l - m)); //note that denom variable here would change

		acc_C(l, m) = factor1 * acc_c_C(l - 1, m);
		acc_S(l, m) = factor1 * acc_c_S(l - 1, m);

		m = l;
		// 13.4.19, 13.4.20
		const double factor2 = -0.5 / m;
		acc_C(m, m) = factor2 * (X * acc_c_C(m - 1, m - 1) - Y * acc_c_S(m - 1, m - 1));
		acc_S(m, m) = factor2 * (Y * acc_c_C(m - 1, m - 1) + X * acc_c_S(m - 1, m - 1));
	}
}

SolidHarmonicsExpansion evaluateLOfR(int order, Vector3<double> r) {
	const bool initializeToZero = false;
	SolidHarmonicsExpansion result(order, initializeToZero);
	result.evaluateLOfR(r);
	return result;
}

void SolidHarmonicsExpansion::evaluateMOfR(Vector3<double> r) {
	const int ord = _order;

	const double r2 = r.L2NormSquare();
	const double inv_r2 = 1.0 / r2;
	const double X = r[0], Y = r[1], Z = r[2];

	// Equation 13.4.16 (a) and (b) Rappaport
	acc_C(0, 0) = sqrt(inv_r2);
	acc_S(0, 0) = 0.0;

	for (int l = 1; l <= ord; ++l) {
		int m;

		// 13.4.21
		for (m = 0; m <= l - 2; ++m) {
			const double factor1 = (2 * l - 1) * Z * inv_r2;
			const double factor2 = -(l - 1 + m) * (l - 1 - m) * inv_r2;

			acc_C(l, m) = factor1 * acc_c_C(l - 1, m) + factor2 * acc_c_C(l - 2, m);
			acc_S(l, m) = factor1 * acc_c_S(l - 1, m) + factor2 * acc_c_S(l - 2, m);
		}

		m = l - 1;
		// 13.4.21 when M(l-2,m) does not exist/is 0
		const double factor1 = (2 * l - 1) * Z * inv_r2;
		acc_C(l, m) = factor1 * acc_c_C(l - 1, m);
		acc_S(l, m) = factor1 * acc_c_S(l - 1, m);

		m = l;
		// 13.4.17, 13.4.18
		const double factor2 = -(2 * m - 1) * inv_r2;
		acc_C(m, m) = factor2 * (X * acc_c_C(m - 1, m - 1) - Y * acc_c_S(m - 1, m - 1));
		acc_S(m, m) = factor2 * (Y * acc_c_C(m - 1, m - 1) + X * acc_c_S(m - 1, m - 1));
	}
}

SolidHarmonicsExpansion evaluateMOfR(int order, Vector3<double> r) {
	const bool initializeToZero = false;
	SolidHarmonicsExpansion result(order, initializeToZero);
	result.evaluateMOfR(r);
	return result;
}

void SolidHarmonicsExpansion::convoluteLL(const SolidHarmonicsExpansion & LE1, const SolidHarmonicsExpansion & LE2) {
	const int ord = (int) _order;

	for (int l = 0; l <= ord; ++l) {                           // l  = 0 : order

		for (int m = 0; m <= l; ++m) {                             // m  = 0 : l
			double temp_result_C = 0.0;
			double temp_result_S = 0.0;

			for (int l_prime = 0; l_prime <= l; ++l_prime) {       // l' = 0 : l
				const int l_diff = l - l_prime;

				const int m_prime_start = std::max(-l_prime, m - l_diff);
				const int m_prime_end = std::min(l_prime, m + l_diff);

				for (int m_prime = m_prime_start; m_prime <= m_prime_end; ++m_prime) { // m' = max (-l', m-(l-l')) : min(l', m+(l-l'))

					const int m_diff = m - m_prime;

					double val_1_c, val_1_s;
					LE1.signed_acc_const_CS(l_prime, m_prime, val_1_c, val_1_s);

					double val_2_c, val_2_s;
					LE2.signed_acc_const_CS(l_diff, m_diff, val_2_c, val_2_s);

					temp_result_C += val_1_c * val_2_c - val_1_s * val_2_s;
					temp_result_S += val_1_s * val_2_c + val_1_c * val_2_s;
				}
			}

			this->acc_C(l, m) = temp_result_C;
			this->acc_S(l, m) = temp_result_S;
		}
	}
}

SolidHarmonicsExpansion convoluteLL(const SolidHarmonicsExpansion & LE1, const SolidHarmonicsExpansion & LE2) {
	const bool initializeToZero = false;
	const int order = LE1.getOrder();
	SolidHarmonicsExpansion L_result(order, initializeToZero);
	L_result.convoluteLL(LE1, LE2);
	return L_result;
}

void SolidHarmonicsExpansion::convoluteLL_Z(const SolidHarmonicsExpansion & LE1, const SolidHarmonicsExpansion & LE2) {
	// assume that LE2(l,m) = 0 everywhere except m = 0;
	// assume that the imaginary part of LE2 = 0;
	const int ord = (int) _order;

	for (int l = 0; l <= ord; ++l) {

		for (int m = 0; m <= l; ++m) {

			double temp_result_C = 0.0;
			double temp_result_S = 0.0;

			const int length = l - m;
			for (int n = 0; n <= length; ++n) {
				const int l_LE1 = m + n;

				const int l_LE2 = length - n;

				double val_1_c, val_1_s;
				LE1.signed_acc_const_CS(l_LE1, m, val_1_c, val_1_s);

				double val_2_c, val_2_s;
				LE2.signed_acc_const_CS(l_LE2, 0, val_2_c, val_2_s);
				temp_result_C += val_1_c*val_2_c;
				temp_result_S += val_1_s*val_2_c;
			}
			this->acc_C(l, m) = temp_result_C;
			this->acc_S(l, m) = temp_result_S;
		}
	}
}

SolidHarmonicsExpansion convoluteLL_Z(const SolidHarmonicsExpansion & LE1, const SolidHarmonicsExpansion & LE2) {
	const bool initializeToZero = true;
	const int order = LE1.getOrder();
	SolidHarmonicsExpansion L_result(order, initializeToZero);
	L_result.convoluteLL_Z(LE1, LE2);
	return L_result;
}

void SolidHarmonicsExpansion::convoluteLM(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME) {
	const int ord = (int) _order;

	for (int l = 0; l <= ord; ++l) {                             // l = 0: order

		for (int m = 0; m <= l; ++m) {                             // m  = 0 : l
			double temp_result_C = 0.0;
			double temp_result_S = 0.0;

			const int l_prime_end = ord - l;
			for (int l_prime = 0; l_prime <= l_prime_end; ++l_prime) { // l' = 0 : (order - l)
				const int l_sum = l + l_prime;

				const int m_prime_start = -l_prime; // max(-l', -m-(l+l')) is always achieved at -l'
				const int m_prime_end = l_prime; // min( l', -m+(l+l')) is always achieved at  l'

				for (int m_prime = m_prime_start; m_prime <= m_prime_end;
						++m_prime) { // m' = -l' : l'

					const int m_sum = m + m_prime;

					double val_L_c, val_L_s;
					LE.signed_acc_const_CS(l_prime, m_prime, val_L_c, val_L_s);

					double val_M_c, val_M_s;
					ME.signed_acc_const_CS(l_sum, m_sum, val_M_c, val_M_s);

					temp_result_C += val_L_c * val_M_c + val_L_s * val_M_s;
					temp_result_S += val_L_c * val_M_s - val_L_s * val_M_c;
				}
			}
			this->acc_C(l, m) = temp_result_C;
			this->acc_S(l, m) = temp_result_S;
		}
	}
}

SolidHarmonicsExpansion convoluteLM(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME) {
	const bool initializeToZero = false;
	const int order = LE.getOrder();
	SolidHarmonicsExpansion M_result(order, initializeToZero);
	M_result.convoluteLM(LE, ME);
	return M_result;
}

void SolidHarmonicsExpansion::convoluteLM_Z(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME) {
	// assume that ME(l,m) = 0 everywhere except m = 0;
	// assume that the imaginary part of ME = 0;
	const int ord = (int) _order;

	for (int l = 0; l <= ord; ++l) {

		for (int m = 0; m <= l; ++m) {
			if (m + l > ord)
				continue; // nothing has to be added
			double temp_result_C = 0.0;
			double temp_result_S = 0.0;

			const int length = ord - l - m;
			for (int n = 0; n <= length; ++n) {
				const int l_ME = l + m + n; // always l+m zeros until the first non-zero element
				const int m_LE = -m;
				const int l_LE = m + n;

				double val_L_c, val_L_s;
				LE.signed_acc_const_CS(l_LE, m_LE, val_L_c, val_L_s);

				double val_M_c, val_M_s;
				ME.signed_acc_const_CS(l_ME, 0, val_M_c, val_M_s);

				temp_result_C += val_L_c * val_M_c;
				temp_result_S -= val_L_s * val_M_c;
			}
			this->acc_C(l, m) = temp_result_C;
			this->acc_S(l, m) = temp_result_S;
		}
	}
}

SolidHarmonicsExpansion convoluteLM_Z(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME) {
	const bool initializeToZero = true;
	const int order = LE.getOrder();
	SolidHarmonicsExpansion M_result(order, initializeToZero);
	M_result.convoluteLM_Z(LE, ME);
	return M_result;
}

void SolidHarmonicsExpansion::convoluteL_ZM(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME) {
	const int ord = (int) _order;

	for (int l = 0; l <= ord; ++l) {

		for (int m = 0; m <= l; ++m) {
			double temp_result_C = 0.0;
			double temp_result_S = 0.0;

			const int length = ord - l;
			for (int n = 0; n <= length; ++n) {
				const int l_LE = n;
				const int l_ME = n + l;
				const int m_ME = m;


				double val_L_c, val_L_s;
				LE.signed_acc_const_CS(l_LE, 0, val_L_c, val_L_s);

				double val_M_c, val_M_s;
				ME.signed_acc_const_CS(l_ME, m_ME, val_M_c, val_M_s);

				temp_result_C += val_L_c * val_M_c;
				temp_result_S += val_L_c * val_M_s;
			}
			this->acc_C(l, m) = temp_result_C;
			this->acc_S(l, m) = temp_result_S;
		}
	}
}

SolidHarmonicsExpansion convoluteL_ZM(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME) {
	const bool initializeToZero = true;
	const int order = LE.getOrder();
	SolidHarmonicsExpansion M_result(order, initializeToZero);
	M_result.convoluteL_ZM(LE, ME);
	return M_result;
}

void SolidHarmonicsExpansion::scaleL(double factor) {
	double f = 1.0;
	for (int l = 0; l <= _order; ++l) {
		for (int m = 0; m <= l; ++m) {
			acc_C(l, m) *= f;
			acc_S(l, m) *= f;
		}
		f *= factor;
	}
}

void SolidHarmonicsExpansion::scaleM(double factor) {
	const double factor_inverse = 1.0 / factor;

	double f = factor_inverse;
	for (int l = 0; l <= _order; ++l) {
		for (int m = 0; m <= l; ++m) {
			acc_C(l, m) *= f;
			acc_S(l, m) *= f;
		}
		f *= factor_inverse;
	}
}

void SolidHarmonicsExpansion::rotatePhi(const double* CosSinPhi, int negate) {
	using std::isnan; // C++11 required
	/* Dachsel Paper eq. (5) */

	for (int l = 0; l <= _order; ++l) {

		for (int m = 0; m <= l; ++m) {
			const double val_M_c = acc_C(l, m);
			const double val_M_s = acc_S(l, m);

			mardyn_assert(!std::isnan(val_M_c));
			mardyn_assert(!std::isnan(val_M_s));
			mardyn_assert(!std::isnan(this->acc_C(l, m)));
			mardyn_assert(!std::isnan(this->acc_S(l, m)));

			acc_C(l, m) = (CosSinPhi[2*m]*val_M_c + negate*CosSinPhi[2*m+1]*val_M_s);
			acc_S(l, m) = (CosSinPhi[2*m]*val_M_s - negate*CosSinPhi[2*m+1]*val_M_c);

			mardyn_assert(!std::isnan(this->acc_C(l, m)));
			mardyn_assert(!std::isnan(this->acc_S(l, m)));
		}
	}
}

void SolidHarmonicsExpansion::convoluteWL(const WignerMatrix& W, const SolidHarmonicsExpansion & LE) {
	/* Dachsel Paper eq. (4) */

	for (int l = 0; l <= _order; ++l) {

		for (int m = 0; m <= l; ++m) {
			double temp_result_C = 0.0;
			double temp_result_S = 0.0;
			for (int k = -l; k<=l; ++k) {
				const double d = W.acc_c(l,m,k); // pre-factor already included

				double val_L_c, val_L_s;
				LE.signed_acc_const_CS(l, k, val_L_c, val_L_s);

				temp_result_C += val_L_c * d;
				temp_result_S += val_L_s * d;
			}

			acc_C(l, m) = temp_result_C;
			acc_S(l, m) = temp_result_S;
		}
	}
}

void SolidHarmonicsExpansion::convoluteWM(const WignerMatrix& W, const SolidHarmonicsExpansion & ME) {
	/* Dachsel Paper eq. (5) */

	for (int l = 0; l <= _order; ++l) {

		for (int m = 0; m <= l; ++m) {
			double temp_result_C = 0.0;
			double temp_result_S = 0.0;
			for (int k = -l; k<=l; ++k) {
				const double d = W.acc_c(l,m,k); // pre-factor already included

				double val_M_c, val_M_s;
				ME.signed_acc_const_CS(l, k, val_M_c, val_M_s);

				temp_result_C += val_M_c * d;
				temp_result_S += val_M_s * d;
			}

			acc_C(l, m) = temp_result_C;
			acc_S(l, m) = temp_result_S;
		}
	}
}

SolidHarmonicsExpansion rotatePhi(const SolidHarmonicsExpansion & E, const double* CosSinPhi, int negate) {

	SolidHarmonicsExpansion M_result(E);
	M_result.rotatePhi(CosSinPhi, negate);
	return M_result;
}

SolidHarmonicsExpansion rotateThetaL(const SolidHarmonicsExpansion & LE, const WignerMatrix& W) {
	mardyn_assert(W.getType() == bhfmm::ROT_TYPE_L);
	const bool initializeToZero = false;
	const int order = LE.getOrder();
	SolidHarmonicsExpansion L_result(order, initializeToZero);
	L_result.convoluteWL(W, LE);
	return L_result;

}

SolidHarmonicsExpansion rotateThetaM(const SolidHarmonicsExpansion & ME, const WignerMatrix& W) {
	mardyn_assert(W.getType() == bhfmm::ROT_TYPE_M);
	const bool initializeToZero = false;
	const int order = ME.getOrder();
	SolidHarmonicsExpansion M_result(order, initializeToZero);
	M_result.convoluteWM(W, ME);
	return M_result;
}


void SolidHarmonicsExpansion::clearMonopole() {
	acc_C_seq(0) = 0.0;
	acc_S_seq(0) = 0.0;
}

void SolidHarmonicsExpansion::setAtMinusR() {
	// Property 13.4.6
	double minus_one_power_l = 1.0;
	for (int l = 0; l <= _order; ++l) {
		for (int m = 0; m <= l; ++m) {
			acc_C(l, m) *= minus_one_power_l;
			acc_S(l, m) *= minus_one_power_l;
		}
		minus_one_power_l *= -1.0;
	}
}

SolidHarmonicsExpansion setAtMinusR(SolidHarmonicsExpansion E) {
	E.setAtMinusR();
	return E;
}

double potentialML(const SolidHarmonicsExpansion & ME, const SolidHarmonicsExpansion & LE) {
	const int ord = std::min(ME.getOrder(),LE.getOrder());

	double u = 0.0;

	// 13.4.10
	for (int l = 0; l <= ord; ++l) {
		int m = 0;
		u += ME.acc_c_C(l, m) * LE.acc_c_C(l, m) + ME.acc_c_S(l, m) * LE.acc_c_S(l, m);

		for (m = 1; m <= l; ++m) {
			u += 2.0 * (ME.acc_c_C(l, m) * LE.acc_c_C(l, m) + ME.acc_c_S(l, m) * LE.acc_c_S(l, m));
		}
	}

	return u;
}

// The following gradient stencil:
//     X 0 X
//       X
// applied on a lower triangular matrix becomes complicated to unroll without if statements.
// I'm leaving it with them.
Vector3<double> forceGradLAndM(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME) {
	// works because of:
	// -potential identity
	// -gradient identity
	// -property about entries with m<0
	// -the fact that S(l,0)=0

	const int ord = (int) LE._order;

	Vector3<double> fc;
	Vector3<double> fs;
	Vector3<double> f = 0.0;

	for (int l = 1; l <= ord; l++) {
		for (int m = 0; m <= l; m++) {
			if (m < l - 1) {
				fc[0] = LE.acc_c_C(l - 1, m + 1);
				fc[1] = LE.acc_c_S(l - 1, m + 1);

				fs[0] = LE.acc_c_S(l - 1, m + 1);
				fs[1] = -LE.acc_c_C(l - 1, m + 1);

			} else {
				fc[0] = 0.0;
				fc[1] = 0.0;

				fs[0] = 0.0;
				fs[1] = 0.0;
			}

			if (m < l) {
				fc[2] = LE.acc_c_C(l - 1, m);
				fs[2] = LE.acc_c_S(l - 1, m);

			} else {
				fc[2] = 0.0;
				fs[2] = 0.0;

			}

			if (m > 0) {
				fc[0] -= LE.acc_c_C(l - 1, m - 1);
				fc[1] += LE.acc_c_S(l - 1, m - 1);
				fc[2] *= 2.0;

				fs[0] -= LE.acc_c_S(l - 1, m - 1);
				fs[1] -= LE.acc_c_C(l - 1, m - 1);
				fs[2] *= 2.0;
			}

			f += fc * ME.acc_c_C(l, m) + fs * ME.acc_c_S(l, m);
		}
	}
	// the force is the negative gradient!
	return f * -1.0;
}

Vector3<double> forceLAndGradM(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME) {
	mardyn_assert(ME._order == LE._order + 1);

	const int ord = (int) ME._order;


	Vector3<double> fc;
	Vector3<double> fs;
	Vector3<double> f = 0.0;

	for (int l = 0; l < ord; l++) {
		for (int m = 0; m <= l; m++) {
			fc[0] = ME.acc_c_C(l + 1, m + 1);
			fc[1] = ME.acc_c_S(l + 1, m + 1);

			fs[0] = ME.acc_c_S(l + 1, m + 1);
			fs[1] = -ME.acc_c_C(l + 1, m + 1);

			fc[2] = -ME.acc_c_C(l + 1, m);
			fs[2] = -ME.acc_c_S(l + 1, m);

			if (m > 0) {
				fc[0] -= ME.acc_c_C(l + 1, m - 1);
				fc[1] += ME.acc_c_S(l + 1, m - 1);
				fc[2] *= 2.0;

				fs[0] -= ME.acc_c_S(l + 1, m - 1);
				fs[1] -= ME.acc_c_C(l + 1, m - 1);
				fs[2] *= 2.0;
			}

			f += LE.acc_c_C(l, m) * fc + LE.acc_c_S(l, m) * fs;
		}
	}
	// the force is the negative gradient!
	return f * -1.0;
}

void SolidHarmonicsExpansion::print() const
{
	int precisionSetting = std::cout.precision( );
	std::ios::fmtflags flagSettings = std::cout.flags();

	std::cout.setf(std::ios::dec | std::ios::showpos | std::ios::showpoint);
	std::cout.precision(5);

	for (int l = 0; l <= _order; ++l) {
		for (int m = 0; m <= l; ++m) {
			std::cout << "(" << std::setw(10) << acc_c_C(l,m) << ", " << std::setw(10) << acc_c_S(l,m) << ")   ";
		}
		std::cout << std::endl;
	}

	std::cout.precision(precisionSetting);
	std::cout.flags(flagSettings);
}

} /* namespace bhfmm */
