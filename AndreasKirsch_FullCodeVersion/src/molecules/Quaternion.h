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
#ifndef QUATERNION_H_
#define QUATERNION_H_

#include "cmath"

/**
 @author Martin Bernreuther
 */
class Quaternion {
public:
	Quaternion(double qw = 0., double qx = 0., double qy = 0., double qz = 0.)
			: m_qw(qw), m_qx(qx), m_qy(qy), m_qz(qz) {
	}

	double qw() const { return m_qw; }
	double qx() const { return m_qx; }
	double qy() const { return m_qy; }
	double qz() const { return m_qz; }
	double magnitude2() {
		return m_qw * m_qw + m_qx * m_qx + m_qy * m_qy + m_qz * m_qz;
	}
	void scale(double s) {
		m_qw *= s;
		m_qx *= s;
		m_qy *= s;
		m_qz *= s;
	}
	void scaleinv(double s) {
		m_qw /= s;
		m_qx /= s;
		m_qy /= s;
		m_qz /= s;
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

	void operator *=(const Quaternion& q);
	void rotate(const double d[3], double drot[3]) const;
	void rotate(double d[3]) const {
		double drot[3];
		rotate(d, drot);
		d[0] = drot[0];
		d[1] = drot[1];
		d[2] = drot[2];
	}
	void rotateinv(const double d[3], double drot[3]) const;
	void rotateinv(double d[3]) const {
		double drot[3];
		rotateinv(d, drot);
		d[0] = drot[0];
		d[1] = drot[1];
		d[2] = drot[2];
	}
	//void differentiate_dbl(const double w[3], Quaternion& dqdt) const;
	void differentiate(const double w[3], Quaternion& dqdt) const;
	//  { differentiate_dbl(w,dqdt); dqdt.scale(.5); }
	void getRotMatrix(double R[3][3]) const;
	void getRotinvMatrix(double R[3][3]) const;

private:
	double m_qw, m_qx, m_qy, m_qz; // components

};

#endif /*QUATERNION_H_*/
