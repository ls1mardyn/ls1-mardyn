/***************************************************************************
 *   Copyright (C) 2010 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
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

#ifndef PARTICLEPAIRSHANDLER_H_
#define PARTICLEPAIRSHANDLER_H_

class RDF;

typedef enum {
    MOLECULE_MOLECULE = 0,      /**< molecule molecule */
    MOLECULE_HALOMOLECULE = 1,  /**< molecule - halo molecule */
    MOLECULE_MOLECULE_FLUID = 2 /**< molecule - molecule (fluid) */
} PairType;


//! @brief interface for defining the action performed when processing a pair
//! @author Martin Buchholz
//! 
//! The idea of a ParticleContainer is, that the container itself only knows
//! about how to efficiently store and access Particles (or neighbouring pairs
//! of Particles) It doesn't know anything (exception follows) about the 
//! internal structure of a Particle, or what action should be performed for 
//! a pair of particles (an action could e.g. be to calculate the force).
//! The only thing that has to be known from a particle is it's position.
//! 
//! An application e.g. wants to find all neighbouring particle pairs and
//! calculate the forces between them. The retrieval of the pairs has to be
//! done by the ParticleContainer, but the force calculation has to be
//! performed somewhere else. That's where this interface comes into play.
//! There are typically three things to be done:
//! - Do some initial stuff before the pair processing
//! - Do something for each pair
//! - Do something after all pairs have been processed
//!
//! The ParticleContainer has an instance of a class implementing this interface
//! as a member variable. When the method of the ParticleContainer that traverses
//! the pairs is called, it first has to call the method "init()" of it's
//! ParticlePairsHandler, then for each pair "processPair(...)" and at the end
//! "finish()".
//! A class implementing this interface now serves as an adapter between the
//! particleContainer an some other part of the programm (e.g. force calculation)
class ParticlePairsHandler {
public:
	//! Constructor
	ParticlePairsHandler() : _rdf( 0 ) {
	}

	//! Destructor
	virtual ~ParticlePairsHandler() {
	}

	//! @brief things to be done before particle pairs are processed
	virtual void init() = 0;

	//! @brief things to be done after particle pairs are processed
	virtual void finish() = 0;

	//! @brief things to be done for each particle pair
	//!
	//! @param particle1 first particle
	//! @param particle2 second particle
	//! @param distanceVector[3] distance between the two particles
	//! @param pairType describes whether the pair is a original pair(0) or a duplicated pair(1)
	//!                 for details about pair types see comments on traversePairs() in ParticleContainer
	virtual double processPair(Molecule& particle1, Molecule& particle2, double distanceVector[3], PairType pairType, double dd, bool calculateLJ) = 0;
	virtual void preprocessTersoffPair(Molecule& particle1, Molecule& particle2, bool pairType) = 0;
	virtual void processTersoffAtom(Molecule& particle1, double params[15], double delta_r) = 0;

	/**
	 * @todo it is not clean to have particleHandlers need to know about the rdf.
	 *       however, this more or less reflects the previous design, so I do it just in the old way
	 *       for the moment.
	 *
	 *       Once we generalized the treatment of the potentials in the ParticleContainer somehow,
	 *       we should make the RDF a particlePairsHandler of its own.
	 */
	//virtual void recordRDF() = 0;
	void setRDF(RDF* rdf) {
		this->_rdf = rdf;
	}

protected:
	RDF* _rdf;
};

#endif /*PARTICLEPAIRSHANDLER_H_*/

