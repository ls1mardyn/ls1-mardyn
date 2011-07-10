/** \file array2d.h
  * \brief 2 dimensional array
  * \author Martin Bernreuther <martin@ipvs.uni-stuttgart.de>
  * \version 0.1
  * \date 03.01.2006
  * license: GPL
*/

/***************************************************************************
 *   Copyright (C) 2006 by Martin Bernreuther   *
 *   Martin.Bernreuther@ipvs.uni-stuttgart.de   *
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

#ifndef ARRAY2D_H_
#define ARRAY2D_H_

#include <vector>
#include <cassert>

template<class T>
class Array2D : private std::vector<T> {
public:
	/// Constructor
	/**
	    parameter dim0  std::size_t dimension, first index
	    parameter dim1  std::size_t dimension, second index
	 **/
	Array2D(std::size_t dim0=0, std::size_t dim1=0)
			: std::vector<T>(dim0*dim1) {
		assert(this->size()==dim0*dim1); m_dim[0]=dim0; m_dim[1]=dim1;
	}
	/// get dimension for first index
	/** retval  std::size_t dimension
	 **/
	std::size_t dim0() const { return m_dim[0]; }
	/// get dimension for second index
	/** retval  std::size_t dimension
	 **/
	std::size_t dim1() const { return m_dim[1]; }
	/// get dimension
	/** retval  std::size_t dimension for dth index
	    parameter d unsigned char dimension
	 **/
	std::size_t dim(unsigned char d) const { assert(d<3); return m_dim[d]; }
	/// get dimension
	/** retval  std::size_t dimension for dth index
	 **/
	template <unsigned char d> std::size_t dim() const { assert(d<3); return m_dim[d]; }

	/// change dimensions
	/**
	    parameter dim0  std::size_t dimension, first index
	    parameter dim1  std::size_t dimension, second index
	 **/
	void redim(std::size_t dim0, std::size_t dim1) {
		this->resize(dim0*dim1);
		assert(this->size()==dim0*dim1);
		m_dim[0]=dim0;
		m_dim[1]=dim1;
	}

	/// access element
	/**
	    parameter i0  std::size_t first index
	    parameter i1  std::size_t second index
	 **/
	T& operator()(std::size_t i0, std::size_t i1) { return (*this)[indices2hash(i0,i1)]; }


	/// Convert cell coordinates to index/hash value
	/** retval  Tindex  index/hash
	    parameter i0  std::size_t first index
	    parameter i1  std::size_t second index
	 **/
	std::size_t indices2hash(std::size_t i0, std::size_t i1) const {
		assert(i0<dim0());
		assert(i1<dim1());

		//return i0*m_dim[1]+i1;  // "column-order"
		return i1*m_dim[0]+i0;  // "row-order"
	}

	/// Convert hash value to indices
	/**
	    parameter h std::size_t hash value
	    parameter i0  std::size_t&  first index
	    parameter i1  std::size_t&  second index
	 **/
	void hash2indices(std::size_t h, std::size_t& i0, std::size_t& i1) const {
		assert(h<this->size());

		// "column-order"
		//i1=h%m_dim[1];
		//i0=h/m_dim[1];
		// "row-order"
		i0=h%m_dim[0];
		i1=h/m_dim[0];
	}


protected:
	std::size_t m_dim[2];
};
#endif /*ARRAY2D_H_*/
