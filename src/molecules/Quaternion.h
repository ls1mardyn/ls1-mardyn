#ifndef QUATERNION_H_
#define QUATERNION_H_

#include <cmath>
#include <array>
#include "utils/mardyn_assert.h"
/**
 @author Martin Bernreuther
 */
class Quaternion {
public:
	Quaternion(double qw = 1., double qx = 0., double qy = 0., double qz = 0.)
			: m_qw(qw), m_qx(qx), m_qy(qy), m_qz(qz) {
	}

	/** @brief create a quaternion, representing a rotation around a rotation axis (vector n) by an angel alpha
	 * @param alpha_rad the rotation angle in radiant
	 * @param n rotation axis represented as a vector (will be normalized)
	 */
	Quaternion(const double& alpha_rad, const std::array<double,3>& n);

	double qw() const { return m_qw; }
	double qx() const { return m_qx; }
	double qy() const { return m_qy; }
	double qz() const { return m_qz; }
	double magnitude2() const {
		return m_qw * m_qw + m_qx * m_qx + m_qy * m_qy + m_qz * m_qz;
	}
	void scale(double s) {
		m_qw *= s;
		m_qx *= s;
		m_qy *= s;
		m_qz *= s;

	}
	void scaleinv(double s) {
		mardyn_assert(s != 0.);
		scale( 1./ s);
	}
	void normalize() {
		scaleinv(sqrt(magnitude2()));
	}
	void conjugate() {
		m_qx = -m_qx;
		m_qy = -m_qy;
		m_qz = -m_qz;
	}
	void inverse() {
		conjugate();
		scaleinv(magnitude2());
	}
	void add(const Quaternion& q) {
		m_qw += q.m_qw;
		m_qx += q.m_qx;
		m_qy += q.m_qy;
		m_qz += q.m_qz;
	}


	void multiply_left(const Quaternion& q);

	void operator *=(const Quaternion& q);
	/**
	 * apply the rotation represented by tis quaternion to d
	 * @param d the vector to be rotated
	 * @return result vector of the rotation
	 */
	std::array<double, 3> rotate(const std::array<double, 3>& d) const;

	/**
	 * apply the rotation represented by tis quaternion to d. The result vector
	 * is stored to d.
	 * @param d the vector to be rotated
	 */
	void rotateInPlace(std::array<double, 3>& d) const {
		std::array<double, 3> dcopy = d;
		d = rotate(dcopy);
	}
	std::array<double, 3> rotateinv(const std::array<double, 3>& d) const;
	void rotateinvInPlace(std::array<double, 3>& d) const {
		std::array<double, 3> dcopy = d;
		d = rotateinv(dcopy);
	}
	//void differentiate_dbl(const double w[3], Quaternion& dqdt) const;
	void differentiate(const std::array<double, 3>& w, Quaternion& dqdt) const;
	//  { differentiate_dbl(w,dqdt); dqdt.scale(.5); }
	void getRotMatrix(double R[3][3]) const;
	void getRotinvMatrix(double R[3][3]) const;

	bool isNormalized() const {
		return fabs(magnitude2() - 1.0) <= 1e-15;
	}
	void check() const{
		using std::isfinite;
		mardyn_assert(std::isfinite(m_qw));
		mardyn_assert(std::isfinite(m_qx));
		mardyn_assert(std::isfinite(m_qy));
		mardyn_assert(std::isfinite(m_qz));
	}

private:
	double m_qw, m_qx, m_qy, m_qz; // components

};

#endif /*QUATERNION_H_*/
