/*
 * SolidHarmonicsExpansion.h
 *
 *  Created on: Nov 21, 2014
 *      Author: tchipevn, rajats
 */

#ifndef SOLIDHARMONICSEXPANSION_H_
#define SOLIDHARMONICSEXPANSION_H_

#include "bhfmm/utils/Vector3.h"
#include "bhfmm/expansions/SolidHarmonicsStorage.h"
#include "bhfmm/utils/WignerMatrix.h"
#ifdef FMM_FFT
#include "bhfmm/fft/FFTAccelerableExpansion.h"
#endif  /* FMM_FFT */
#include <cstdlib>
#include <cmath>

namespace bhfmm {
class SolidHarmonicsExpansion;

/**
 * swap function for the copy and swap idiom
 * @param s1
 * @param s2
 */
void swap(SolidHarmonicsExpansion & s1, SolidHarmonicsExpansion & s2);

/**
 * the potential due to an M and an L expansion
 * implements 13.4.10 @see SolidHarmonicsExpansion
 *
 * @param ME the M expansion
 * @param LE the L expansion
 * @return potential
 */
double potentialML(const SolidHarmonicsExpansion & ME, const SolidHarmonicsExpansion & LE);

/**
 * the force due to a local expansion: needed for L2P
 * implements 13.4.33-35
 *
 * @param LE gradient of L expansion is calculated
 * @param ME the M expansion
 * @return force vector
 */
Vector3<double> forceGradLAndM(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME);

/**
 * the force due to a multipole expansion: needed for M2P
 * The formula is not given in Rapaport. Took it from
 * Perez-Jorda and Yang.
 * "A concise redefinition of the solid spherical harmonics and its use in fast multipole methods".
 * In: The Journal of chemical physics 104.20 (1996)
 *
 * @param LE L expansion
 * @param ME gradient of M expansion is calculated
 * @return
 */
Vector3<double> forceLAndGradM(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME);

/**
 * Implements a multipole or local expansion in Solid Harmonics.
 * All formulas follow the notation and code of D.C. Rapaport,
 * "The Art of Molecular Dynamics Simulation", Chapter 13.4.
 *
 * Makes use of two triangular matrices (SolidHarmonicsStorage objects)
 * for the real and imaginary part (or (-imaginary), see formula 13.4.9 in Rapaport).
 */
class SolidHarmonicsExpansion
#ifdef FMM_FFT
		: public FFTAccelerableExpansion
#endif  /* FMM_FFT */
{
public:
	/**
	 * constructor
	 * @param order order of the expansion
	 * @param initializeToZero if true, values are set to zero
	 */
	SolidHarmonicsExpansion(int order, bool initializeToZero = true);

	/**
	 * copy constructor
	 * @param rhs
	 */
	SolidHarmonicsExpansion(const SolidHarmonicsExpansion & rhs);

	/**
	 * destructor
	 */
	virtual ~SolidHarmonicsExpansion();

	/**
	 * swap function for the copy and swap idiom
	 * @param s1
	 * @param s2
	 */
	friend void swap(SolidHarmonicsExpansion & s1, SolidHarmonicsExpansion & s2);

	// OPERATORS //
	/**
	 * operator= entrywise copy
	 * @param RHS
	 * @return
	 */
	SolidHarmonicsExpansion & operator =(SolidHarmonicsExpansion RHS);

	/**
	 * operator+= entrywise addition
	 * @param RHS
	 * @return
	 */
	SolidHarmonicsExpansion & operator+=(const SolidHarmonicsExpansion & RHS);

	/**
	 * operator *= scalar multiplication
	 * @param scalar
	 * @return
	 */
	SolidHarmonicsExpansion & operator*=(double scalar);

	// METHODS //
	/**
	 * clear expansion - set all entries to zero
	 */
	void clear();

	/**
	 * evaluate L expansion at a given position
	 * implements 13.4.16-20
	 *
	 * @param r the position vector
	 */
	void evaluateLOfR(Vector3<double> r);

	/**
	 * evaluate M expansion at a given position
	 * implements 13.4.16-18
	 * @param r the position vector
	 */
	void evaluateMOfR(Vector3<double> r);

	/**
	 * convolution of two L-type expansions.
	 * result is again an L expansion.
	 * needed for M2M
	 * implements 13.4.38
	 *
	 * @param LE1
	 * @param LE2
	 */
	void convoluteLL(const SolidHarmonicsExpansion & LE1, const SolidHarmonicsExpansion & LE2);

	/**
	 * convolution of two L-type expansions, where LE2 in considered
	 * to be evaluated at a direction vector that is parallel to
	 * the z-axis.
	 * result is again an L expansion.
	 * needed for M2M
	 * implements 13.4.38
	 *
	 * @param LE1
	 * @param LE2
	 */
	void convoluteLL_Z(const SolidHarmonicsExpansion & LE1, const SolidHarmonicsExpansion & LE2);

	/**
	 * convolution of an L and an M expansion
	 * result is an M expansion.
	 * needed for M2M, L2L
	 * implements 13.4.42
	 *
	 * @param LE
	 * @param ME
	 */
	void convoluteLM(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME);

	/**
	 * convolution of an L and an M expansion, where M in considered
	 * to be evaluated at a direction vector that is parallel to
	 * the z-axis.
	 * result is again an L expansion.
	 * needed for M2L
	 *
	 * @param LE1
	 * @param LE2
	 */
	void convoluteLM_Z(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME);

	/**
	 * convolution of an L and an M expansion, where L in considered
	 * to be evaluated at a direction vector that is parallel to
	 * the z-axis.
	 * result is again an L expansion.
	 * needed for M2M
	 *
	 * @param LE1
	 * @param LE2
	 */
	void convoluteL_ZM(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME);

	/**
	 * scale an L expansion to a different pseudoparticle size
	 * needed for periodic boundary conditions (coming soon)
	 *
	 * @param factor the factor by which the expansion should be scaled
	 */
	void scaleL(double factor);

	/**
	 * scale an L expansion to a different pseudoparticle size
	 * needed for periodic boundary conditions (coming soon)
	 *
	 * @param factor
	 */
	void scaleM(double factor);

	void rotatePhi(const double* CosSinPhi, int negate);

	void convoluteWL(const WignerMatrix & W, const SolidHarmonicsExpansion & LE);

	void convoluteWM(const WignerMatrix & W, const SolidHarmonicsExpansion & ME);

	/**
	 * set monopole entry to zero
	 * needed for periodic boundary conditions
	 */
	void clearMonopole();

	/**
	 * given an L or M expansion, evaluated at position (R),
	 * set it to the value it would have, if evaluated at (-R)
	 * implements 13.4.6 from Rapaport.
	 */
	void setAtMinusR();

	/**
	 * the potential due to an M and an L expansion
	 * implements 13.4.10 @see SolidHarmonicsExpansion
	 *
	 * @param ME the M expansion
	 * @param LE the L expansion
	 * @return potential
	 */
	friend double potentialML(const SolidHarmonicsExpansion & ME, const SolidHarmonicsExpansion & LE);

	/**
	 * the force due to a local expansion: needed for L2P
	 * implements 13.4.33-35
	 *
	 * @param LE gradient of L expansion is calculated
	 * @param ME the M expansion
	 * @return force vector
	 */
	friend Vector3<double> forceGradLAndM(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME);

	/**
	 * the force due to a multipole expansion: needed for M2P
	 * The formula is not given in Rapaport. Took it from
	 * Perez-Jorda and Yang.
	 * "A concise redefinition of the solid spherical harmonics and its use in fast multipole methods".
	 * In: The Journal of chemical physics 104.20 (1996)
	 *
	 * @param LE L expansion
	 * @param ME gradient of M expansion is calculated
	 * @return
	 */
	friend Vector3<double> forceLAndGradM(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME);

	/**
	 * @return the order of the expansion
	 */
	int getOrder() const {
		return _order;
	}

	/**
	 * Access to C values directly. Intended to be used only for MPI communication.
	 * @param l
	 * @param m
	 * @return value at (l,m)
	 */
	double & getC(int l, int m) {
		return acc_C(l,m);
	}
  //for FFTAccelerableExpansion
  double & get_C(unsigned l, unsigned m) {
		return acc_C(l,m);
	}

	/**
	 * Access to S values directly. Intended to be used only for MPI communication.
	 * @param l
	 * @param m
	 * @return value at (l,m)
	 */
	double & getS(int l, int m) {
		return acc_S(l,m);
	}
  //for FFTAccelerableExpansion
  double & get_S(unsigned l, unsigned m) {
		return acc_S(l,m);
	}

	int getNumEntries() const {
		return _c.getTotalNumValues() + _s.getTotalNumValues();
	}

	/**
	 * For debugging
	 */
	void print() const;

	/**
	 * write values of expansion to a double buffer for sending through MPI
	 * @param buf buffer to write to
	 * @param position index at which next entry should be written, note that it's value is updated!
	 */
	void writeValuesToMPIBuffer(std::vector<double> buf, int& position) const {
		const int end = _c.getTotalNumValues();
		for (int i = 0; i < end; ++i) {
			buf[position++] = acc_c_C_seq(i);
		}
		for (int i = 0; i < end; ++i) {
			buf[position++] = acc_c_S_seq(i);
		}
	}

	void readValuesFromMPIBuffer(std::vector<double> buf, int& position) {
		const int end = _c.getTotalNumValues();
		for (int i = 0; i < end; ++i) {
			acc_C_seq(i) = buf[position++];
		}
		for (int i = 0; i < end; ++i) {
			acc_S_seq(i) = buf[position++];
		}
	}
	void addValuesFromMPIBuffer(std::vector<double> buf, int& position) {
		const int end = _c.getTotalNumValues();
		for (int i = 0; i < end; ++i) {
			acc_C_seq(i) += buf[position++];
		}
		for (int i = 0; i < end; ++i) {
			acc_S_seq(i) += buf[position++];
		}
	}

private:
	//private accessors to C and S terms

	/**
	 * @param l order
	 * @param m phase
	 * @return C term at entry (l,m) by reference
	 */
	double & acc_C(int l, int m) {
		return _c.getValue(l, m);
	}

	/**
	 * @param l order
	 * @param m phase
	 * @return S term at entry (l,m) by reference
	 */
	double & acc_S(int l, int m) {
		return _s.getValue(l, m);
	}

	/**
	 * @param r
	 * @return C term at sequential entry (r) by reference
	 */
	double & acc_C_seq(int r) {
		return _c.getValueSequential(r);
	}

	/**
	 * @param r
	 * @return S term at sequential entry (r) by reference
	 */
	double & acc_S_seq(int r) {
		return _s.getValueSequential(r);
	}

	/**
	 * @param l order
	 * @param m phase
	 * @return C term at entry (l,m) by value
	 */
	double acc_c_C(int l, int m) const {
		return _c.getValueConst(l, m);
	}

	/**
	 * @param l order
	 * @param m phase
	 * @return S term at entry (l,m) by value
	 */
	double acc_c_S(int l, int m) const {
		return _s.getValueConst(l, m);
	}

	/**
	 * @param r
	 * @return C term at sequential entry (r) by value
	 */
	double acc_c_C_seq(int r) const {
		return _c.getValueConstSequential(r);
	}

	/**
	 * @param r
	 * @return S term at sequential entry (r) by value
	 */
	double acc_c_S_seq(int r) const {
		return _s.getValueConstSequential(r);
	}

	/**
	 * "signed" accessor: gives access to (l,m) with negative phase
	 * by property 13.4.11-12
	 *
	 * @param l order
	 * @param m phase (can also be negative here !)
	 * @param c returned C part at entry (l,m)
	 * @param s returned S part at entry (l,m)
	 */
	void signed_acc_const_CS(int l, int m, double & c, double & s) const {
		const int ind = _c.index(l, abs(m));
		c = _c.getValueConstSequential(ind);
		s = _s.getValueConstSequential(ind);
		if (m >= 0) {
			return;
		} else {
			const int minus_one_pow_m = (m & 1) ? -1.0 : 1.0;
			c *= minus_one_pow_m;
			s *= -minus_one_pow_m;
		}
	}

protected:
	int _order;

private:
	/**
	 * C terms of an expansion:
	 * real part
	 * see 13.4.8-9
	 */
	SolidHarmonicsStorage _c;

	/**
	 * S terms of an expansion:
	 * imaginary part of an M expansion or (-imaginary) part of an L expansion
	 * see 13.4.8-9
	 */
	SolidHarmonicsStorage _s; // S-terms: imaginary or (-imaginary) part, see 13.4.8

};

/**
 * entrywise addition
 * @param LHS
 * @param RHS
 * @return
 */
SolidHarmonicsExpansion operator+(SolidHarmonicsExpansion LHS, const SolidHarmonicsExpansion & RHS);

/**
 * scalar multiplication
 * @param scalar
 * @param RHS
 * @return
 */
SolidHarmonicsExpansion operator*(double scalar, SolidHarmonicsExpansion RHS);

/**
 * @see SolidHarmonicsExpansion::evaluateLOfR
 * @param order
 * @param r
 * @return
 */
SolidHarmonicsExpansion evaluateLOfR(int order, Vector3<double> r);

/**
 * @see SolidHarmonicsExpansion::evaluateMOfR
 * @param order
 * @param r
 * @return
 */
SolidHarmonicsExpansion evaluateMOfR(int order, Vector3<double> r);

/**
 * @see SolidHarmonicsExpansion::convoluteLL
 * @param LE1
 * @param LE2
 * @return
 */
SolidHarmonicsExpansion convoluteLL(const SolidHarmonicsExpansion & LE1, const SolidHarmonicsExpansion & LE2);

/**
 * @see SolidHarmonicsExpansion::convoluteLL_Z
 * @param LE1
 * @param LE2
 * @return
 */
SolidHarmonicsExpansion convoluteLL_Z(const SolidHarmonicsExpansion & LE1, const SolidHarmonicsExpansion & LE2);

/**
 * @see SolidHarmonicsExpansion::convoluteLM
 * @param LE
 * @param ME
 * @return
 */
SolidHarmonicsExpansion convoluteLM(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME);

/**
 * @see SolidHarmonicsExpansion::convoluteLM_Z
 * @param LE
 * @param ME
 * @return
 */
SolidHarmonicsExpansion convoluteLM_Z(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME);

/**
 * @see SolidHarmonicsExpansion::convoluteL_ZM
 * @param LE
 * @param ME
 * @return
 */
SolidHarmonicsExpansion convoluteL_ZM(const SolidHarmonicsExpansion & LE, const SolidHarmonicsExpansion & ME);

/**
 * @see SolidHarmonicsExpansion::setAtMinusR
 * @param E
 * @return
 */

SolidHarmonicsExpansion rotatePhi(const SolidHarmonicsExpansion & E, const double* CosSinPhi, int negate);

SolidHarmonicsExpansion rotateThetaL(const SolidHarmonicsExpansion & LE, const WignerMatrix& W);

SolidHarmonicsExpansion rotateThetaM(const SolidHarmonicsExpansion & ME, const WignerMatrix& W);

SolidHarmonicsExpansion setAtMinusR(SolidHarmonicsExpansion E);

} /* namespace bhfmm */

#endif /* SOLIDHARMONICSEXPANSION_H_ */
