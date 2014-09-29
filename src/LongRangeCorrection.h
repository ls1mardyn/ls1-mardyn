/*************************************************************************
 * Copyright (C) 2012 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
 *                                                                       *
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or (at *
 * your option) any later version.                                       *
 *                                                                       *
 * This program is distributed in the hope that it will be useful, but   *
 * WITHOUT ANY WARRANTY; without even the implied warranty of            * 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      *
 * General Public License for more details.                              *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program; if not, write to the Free Software           *
 * Foundation, 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.   *
 *************************************************************************/

#ifndef __LONGRANGECORRECTION_H__
#define __LONGRANGECORRECTION_H__

#include <cmath>

//class Domain;
//class Planar;
//class Homogeneous;

class LongRangeCorrection{

public:
	LongRangeCorrection() {}
//	~LongRangeCorrection() {}
//	void initializeLongRange();
	virtual void calculateLongRange() = 0;
/*
private:
	unsigned _type;
	Planar* _planar;
	Homogeneous* _homogen;
*/
  
};

#endif /* __HOMOGENEOUS_H__ */
