#include "Quaternion.h"
#include <cmath>

#define VECTOR3_ZERO ((std::array<double, 3>){ 0, 0, 0 })

Quaternion::Quaternion(const double& alpha_rad, const std::array<double,3>& n)
{
	// normalize rotation axis vector n
	double length = 0.0;
	for(uint16_t d=0; d<3; ++d)
		length += n.at(d)*n.at(d);
	double inv_length = 1./sqrt(length);
	std::array<double,3> n_norm;
	for(uint16_t d=0; d<3; ++d)
		n_norm.at(d) = n.at(d) * inv_length;

	double alpha_rad_half = 0.5*alpha_rad;
	double sin_alpha_rad_half = sin(alpha_rad_half);
	m_qw = cos(alpha_rad_half);
	m_qx = sin_alpha_rad_half * n_norm.at(0);
	m_qy = sin_alpha_rad_half * n_norm.at(1);
	m_qz = sin_alpha_rad_half * n_norm.at(2);
}

void Quaternion::operator *=(const Quaternion& q) {
	double qw = m_qw*q.m_qw - m_qx*q.m_qx - m_qy*q.m_qy - m_qz*q.m_qz;
	double qx = m_qw*q.m_qx + m_qx*q.m_qw + m_qy*q.m_qz - m_qz*q.m_qy;
	double qy = m_qw*q.m_qy + m_qy*q.m_qw + m_qz*q.m_qx - m_qx*q.m_qz;
	m_qz = m_qw*q.m_qz + m_qz*q.m_qw + m_qx*q.m_qy - m_qy*q.m_qx;
	m_qy = qy;
	m_qx = qx;
	m_qw = qw;
}

void Quaternion::multiply_left(const Quaternion& q) {
	double qw = q.m_qw*m_qw - q.m_qx*m_qx - q.m_qy*m_qy - q.m_qz*m_qz;
	double qx = q.m_qw*m_qx + q.m_qx*m_qw + q.m_qy*m_qz - q.m_qz*m_qy;
	double qy = q.m_qw*m_qy + q.m_qy*m_qw + q.m_qz*m_qx - q.m_qx*m_qz;
	m_qz = q.m_qw*m_qz + q.m_qz*m_qw + q.m_qx*m_qy - q.m_qy*m_qx;
	m_qy = qy;
	m_qx = qx;
	m_qw = qw;
}

std::array<double, 3> Quaternion::rotate(const std::array<double, 3>& d) const {
	std::array<double, 3> drot VECTOR3_ZERO;
	if(VECTOR3_ZERO != d) { // GH#267
		double ww = m_qw*m_qw;
		double xx = m_qx*m_qx;
		double yy = m_qy*m_qy;
		double zz = m_qz*m_qz;
		double wx = m_qw*m_qx;
		double wy = m_qw*m_qy;
		double wz = m_qw*m_qz;
		double xy = m_qx*m_qy;
		double xz = m_qx*m_qz;
		double yz = m_qy*m_qz;
		//          1-2*(yy+zz)
		drot[0] = (ww+xx-yy-zz)*d[0] + 2.*(xy-wz)*d[1] + 2.*(wy+xz)*d[2];
		//                            1-2*(xx+zz)
		drot[1] = 2.*(wz+xy)*d[0] + (ww-xx+yy-zz)*d[1] + 2.*(yz-wx)*d[2];
		//                                              1-2*(xx+yy)
		drot[2] = 2.*(xz-wy)*d[0] + 2.*(wx+yz)*d[1] + (ww-xx-yy+zz)*d[2];
	}
	return drot;
}

std::array<double, 3> Quaternion::rotateinv(const std::array<double, 3>& d) const {
	std::array<double, 3> drot VECTOR3_ZERO;
	if (VECTOR3_ZERO != d) { // GH#267
		double ww = m_qw*m_qw;
		double xx = m_qx*m_qx;
		double yy = m_qy*m_qy;
		double zz = m_qz*m_qz;
		double wx = m_qw*m_qx;
		double wy = m_qw*m_qy;
		double wz = m_qw*m_qz;
		double xy = m_qx*m_qy;
		double xz = m_qx*m_qz;
		double yz = m_qy*m_qz;
		//          1-2*(yy+zz)
		drot[0] = (ww+xx-yy-zz)*d[0] + 2.*(xy+wz)*d[1] + 2.*(xz-wy)*d[2];
		//                            1-2*(xx+zz)
		drot[1] = 2.*(xy-wz)*d[0] + (ww-xx+yy-zz)*d[1] + 2.*(yz+wx)*d[2];
		//                                              1-2*(xx+yy)
		drot[2] = 2.*(xz+wy)*d[0] + 2.*(yz-wx)*d[1] + (ww-xx-yy+zz)*d[2];
	}
	return drot;
}

/*
 void Quaternion::differentiate_dbl(const double w[3], Quaternion& dqdt) const {
 dqdt.m_qw=-m_qx*w[0]-m_qy*w[1]-m_qz*w[2];
 dqdt.m_qx= m_qw*w[0]-m_qz*w[1]+m_qy*w[2];
 dqdt.m_qy= m_qz*w[0]+m_qw*w[1]-m_qx*w[2];
 dqdt.m_qz=-m_qy*w[0]+m_qx*w[1]+m_qw*w[2];
 }
 */

void Quaternion::differentiate(const std::array<double, 3>& w, Quaternion& dqdt) const {
	if(VECTOR3_ZERO == w) { // GH#267
		dqdt.m_qw = 0; dqdt.m_qx = 0; dqdt.m_qy = 0; dqdt.m_qz = 0;
		return;
	}

	dqdt.m_qw = .5 * ( -m_qx*w[0] - m_qy*w[1] - m_qz*w[2] );
	dqdt.m_qx = .5 * (  m_qw*w[0] - m_qz*w[1] + m_qy*w[2] );
	dqdt.m_qy = .5 * (  m_qz*w[0] + m_qw*w[1] - m_qx*w[2] );
	dqdt.m_qz = .5 * ( -m_qy*w[0] + m_qx*w[1] + m_qw*w[2] );
}

void Quaternion::getRotMatrix(double R[3][3]) const {
	double ww = m_qw*m_qw;
	double xx = m_qx*m_qx;
	double yy = m_qy*m_qy;
	double zz = m_qz*m_qz;
	double wx = m_qw*m_qx;
	double wy = m_qw*m_qy;
	double wz = m_qw*m_qz;
	double xy = m_qx*m_qy;
	double xz = m_qx*m_qz;
	double yz = m_qy*m_qz;
	R[0][0] = ww+xx-yy-zz;
	R[0][1] = 2.*(xy-wz);
	R[0][2] = 2.*(xz+wy);
	R[1][0] = 2.*(xy+wz);
	R[1][1] = ww-xx+yy-zz;
	R[1][2] = 2.*(yz-wx);
	R[2][0] = 2.*(xz-wy);
	R[2][1] = 2.*(yz+wx);
	R[2][2] = ww-xx-yy+zz;
}

void Quaternion::getRotinvMatrix(double R[3][3]) const {
	double ww = m_qw*m_qw;
	double xx = m_qx*m_qx;
	double yy = m_qy*m_qy;
	double zz = m_qz*m_qz;
	double wx = m_qw*m_qx;
	double wy = m_qw*m_qy;
	double wz = m_qw*m_qz;
	double xy = m_qx*m_qy;
	double xz = m_qx*m_qz;
	double yz = m_qy*m_qz;
	R[0][0] = ww+xx-yy-zz;
	R[0][1] = 2.*(xy+wz);
	R[0][2] = 2.*(xz-wy);
	R[1][0] = 2.*(xy-wz);
	R[1][1] = ww-xx+yy-zz;
	R[1][2] = 2.*(yz+wx);
	R[2][0] = 2.*(xz+wy);
	R[2][1] = 2.*(yz-wx);
	R[2][2] = ww-xx-yy+zz;
}
