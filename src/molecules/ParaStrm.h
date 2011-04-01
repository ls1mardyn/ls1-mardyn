/** \file parastrm.h
  * \brief Parameter stream
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
#ifndef PARASTRM_H_
#define PARASTRM_H_

#include <cstdlib>
#include <cassert>
#include <string.h>
#include <cstddef>

/** Parameter stream
 *
 * @author Martin Bernreuther
 */
class ParaStrm {
public:
	//enum Eparatype {UNDEF=0,DOUBLE,INT};

    /// Constructor generating a new, empty parameter stream.
	ParaStrm() : m_size(0)/*, m_pos(0)*/, m_pstrm(NULL), m_readpos(NULL) { }

    /// Constructor using an existing parameter stream (resets the reading "pointer" to the beginning of the stream).
	ParaStrm( const ParaStrm& param_stream ) : m_size( param_stream.m_size ), m_pstrm(NULL), m_readpos(NULL) {
        m_pstrm = (char *) malloc(m_size);
        assert(m_pstrm);
        memcpy( m_pstrm, param_stream.m_pstrm, m_size );
        m_readpos = m_pstrm;
	}

    /// Destructor
	~ParaStrm() {
		if (m_pstrm)
			free(m_pstrm);
	}

	/// end of stream reached?
	bool eos() const {
		return m_readpos >= m_pstrm + m_size;
	}
	/// stream empty?
	bool empty() const {
		return m_size == 0;
	}
	/// reset reading "pointer" to the beginning of the stream
	void reset_read() {
		m_readpos = m_pstrm; /*m_pos=0;*/
	}

    /// Copy parameter stream (resets the reading "pointer" to the beginning of the stream).
    ParaStrm& operator=( const ParaStrm& param_stream ) {
        m_size  = param_stream.m_size;
        m_pstrm = (char *) malloc(m_size);
        assert(m_pstrm);
        memcpy( m_pstrm, param_stream.m_pstrm, m_size );
        m_readpos = m_pstrm;
        return *this;
    }

	/// add value to the (end of) stream
	template<class T> ParaStrm& operator <<(T p) {
		ptrdiff_t oldsize = m_size;
		// determine relative reading position
		ptrdiff_t readpos = m_readpos - m_pstrm;
		// enlarge vector
		m_size += sizeof(p);
		if (!m_pstrm)
			m_pstrm = (char*) malloc(m_size);
		else
			m_pstrm = (char*) realloc(m_pstrm, m_size);
		// realloc might change m_pstrm!
		assert(m_pstrm);
		// save value
		*((T*) (m_pstrm + oldsize)) = p;
		//m_ptypes.push_back(UNDEF);
		// adjust m_readpos to maybe changed m_pstrm
		m_readpos = m_pstrm + readpos;
		return *this;
	}

	/// read value at actual stream position and move position
	template<class T> ParaStrm& operator >>(T& p) {
		assert(m_readpos+sizeof(T)<=m_pstrm+m_size);
		// get value
		p = *((T*) m_readpos);
		// move reading position
		m_readpos = (char*) ((T*) m_readpos + 1);
		//++m_pos;
		return *this;
	}

private:
	size_t m_size;  /**< Size of stored stream data in bytes. */
	char *m_pstrm;  /**< Base pointer of data storage. */
	char *m_readpos;  /**< Pointer to the current data item. */
	// we also could check for valid types...
	//std::vector<Eparatype> m_ptypes;
	//std::size_t m_pos;
};

#endif /*PARASTRM_H_*/
