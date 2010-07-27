/***************************************************************************
 *   Copyright (C) 2005 by Martin Bernreuther   *
 *   Martin.Bernreuther@informatik.uni-stuttgart.de   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "molecules/Quaternion.h"

void Quaternion::operator *=(const Quaternion& q) {
	double qw=m_qw*q.m_qw-m_qx*q.m_qx-m_qy*q.m_qy-m_qz*q.m_qz;
	double qx=m_qw*q.m_qx+m_qx*q.m_qw+m_qy*q.m_qz-m_qz*q.m_qy;
	double qy=m_qw*q.m_qy+m_qy*q.m_qw+m_qz*q.m_qx-m_qx*q.m_qz;
	m_qz=m_qw*q.m_qz+m_qz*q.m_qw+m_qx*q.m_qy-m_qy*q.m_qx;
	m_qy=qy;
	m_qx=qx;
	m_qw=qw;
}

void Quaternion::rotate(const double d[3], double drot[3]) const {
	double ww=m_qw*m_qw;
	double xx=m_qx*m_qx;
	double yy=m_qy*m_qy;
	double zz=m_qz*m_qz;
	double xy=m_qx*m_qy;
	double zw=m_qz*m_qw;
	double xz=m_qx*m_qz;
	double yw=m_qy*m_qw;
	//       1-2*(yy+zz)
	drot[0]=(ww+xx-yy-zz)*d[0]+2.*(xy+zw)*d[1]+2.*(xz-yw)*d[2];
	double yz=m_qy*m_qz;
	double xw=m_qx*m_qw;
	//                       1-2*(xx+zz)
	drot[1]=2.*(xy-zw)*d[0]+(ww-xx+yy-zz)*d[1]+2.*(yz+xw)*d[2];
	//                                       1-2*(xx+yy)
	drot[2]=2.*(xz+yw)*d[0]+2.*(yz-xw)*d[1]+(ww-xx-yy+zz)*d[2];
}

void Quaternion::rotateinv(const double d[3], double drot[3]) const {
	double ww=m_qw*m_qw;
	double xx=m_qx*m_qx;
	double yy=m_qy*m_qy;
	double zz=m_qz*m_qz;
	double xy=m_qx*m_qy;
	double zw=m_qz*m_qw;
	double xz=m_qx*m_qz;
	double yw=m_qy*m_qw;
	//       1-2*(yy+zz)
	drot[0]=(ww+xx-yy-zz)*d[0]+2.*(xy-zw)*d[1]+2.*(xz+yw)*d[2];
	double yz=m_qy*m_qz;
	double xw=m_qx*m_qw;
	//                       1-2*(xx+zz)
	drot[1]=2.*(xy+zw)*d[0]+(ww-xx+yy-zz)*d[1]+2.*(yz-xw)*d[2];
	//                                       1-2*(xx+yy)
	drot[2]=2.*(xz-yw)*d[0]+2.*(yz+xw)*d[1]+(ww-xx-yy+zz)*d[2];
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
	dqdt.m_qw=.5*(-m_qx*w[0]-m_qy*w[1]-m_qz*w[2]);
	dqdt.m_qx=.5*( m_qw*w[0]-m_qz*w[1]+m_qy*w[2]);
	dqdt.m_qy=.5*( m_qz*w[0]+m_qw*w[1]-m_qx*w[2]);
	dqdt.m_qz=.5*(-m_qy*w[0]+m_qx*w[1]+m_qw*w[2]);
}

void Quaternion::getRotMatrix(double R[3][3]) const {
	double ww=m_qw*m_qw;
	double xx=m_qx*m_qx;
	double yy=m_qy*m_qy;
	double zz=m_qz*m_qz;
	double xy=m_qx*m_qy;
	double zw=m_qz*m_qw;
	double xz=m_qx*m_qz;
	double yw=m_qy*m_qw;
	R[0][0]=ww+xx-yy-zz;
	R[0][1]=2.*(xy+zw);
	R[0][2]=2.*(xz-yw);
	double yz=m_qy*m_qz;
	double xw=m_qx*m_qw;
	R[1][0]=2.*(xy-zw);
	R[1][1]=ww-xx+yy-zz;
	R[1][2]=2.*(yz+xw);
	R[2][0]=2.*(xz+yw);
	R[2][1]=2.*(yz-xw);
	R[2][2]=ww-xx-yy+zz;
}

void Quaternion::getRotinvMatrix(double R[3][3]) const {
	double ww=m_qw*m_qw;
	double xx=m_qx*m_qx;
	double yy=m_qy*m_qy;
	double zz=m_qz*m_qz;
	double xy=m_qx*m_qy;
	double zw=m_qz*m_qw;
	double xz=m_qx*m_qz;
	double yw=m_qy*m_qw;
	R[0][0]=ww+xx-yy-zz;
	R[0][1]=2.*(xy-zw);
	R[0][2]=2.*(xz+yw);
	double yz=m_qy*m_qz;
	double xw=m_qx*m_qw;
	R[1][0]=2.*(xy+zw);
	R[1][1]=ww-xx+yy-zz;
	R[1][2]=2.*(yz-xw);
	R[2][0]=2.*(xz-yw);
	R[2][1]=2.*(yz+xw);
	R[2][2]=ww-xx-yy+zz;
}
