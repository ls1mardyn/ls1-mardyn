#include "Quaternion.h"

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

void Quaternion::rotate(const double d[3], double drot[3]) const {
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

void Quaternion::rotateinv(const double d[3], double drot[3]) const {
	//cout<<"quat after "<<m_qw<<" "<<m_qx<<" "<<m_qy<<" "<<m_qz<<endl;
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

/*
 void Quaternion::differentiate_dbl(const double w[3], Quaternion& dqdt) const {
 dqdt.m_qw=-m_qx*w[0]-m_qy*w[1]-m_qz*w[2];
 dqdt.m_qx= m_qw*w[0]-m_qz*w[1]+m_qy*w[2];
 dqdt.m_qy= m_qz*w[0]+m_qw*w[1]-m_qx*w[2];
 dqdt.m_qz=-m_qy*w[0]+m_qx*w[1]+m_qw*w[2];
 }
 */

void Quaternion::differentiate(const double w[3], Quaternion& dqdt) const {
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
