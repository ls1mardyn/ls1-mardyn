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

#ifndef PARTICLECONTAINER_H_
#define PARTICLECONTAINER_H_

#include <list>

class CellProcessor;
class ChemicalPotential;
class Domain;
class DomainDecompBase;
class Molecule;
class ParticleContainer;
class ParticlePairsHandler;
class XMLfileUnits;

//! @brief This Interface is used to get access to particles and pairs of particles
//! @author Martin Buchholz
//!
//! A particleContainer is used to store Particles with a short-range potential
//! in a way that the access to pairs of neighbouring particles is efficient.
//! Neighbouring particles are particles which have a distance shorter than
//! a given cutoff radius.
//! A common task when using a PariticleContainer is to do something for all
//! particles. This can be done using the methods begin, next and end, e.g.:
//! \code
//! ParticleContainer* partContPtr;
//! Molecule* partPtr;
//! for(partPtr = partContPtr->begin(); partPtr != partContPtr->end(); partPtr = partContPtr->next()){
//!   partPtr->doSomething();
//! }
//! \endcode
//!
//!
//! Particles stored in this container can either belong to the local process
//! or they can be duplicates (from neighbouring processes or periodic boundary).
//! As the simulated regions are often cuboids, a bounding box is defined through
//! two opposing corners (_boundingBoxMin[3] and _boundingBoxMax[3]).
//! Particles inside this bounding box belong to this process, those outside don't.
//! An exception to this is when particles are moved in a time step. It has to be
//! ensured that particles which leave the bounding box are properly handled.
//!
//! For non-cuboid regions, the bounding box still has to be defined as it gives
//! an approximation for the region that is covered by the ParticleContainer.
//!
//! This interface doesn't implement the datastructure, it just tells which
//! methods a class implementing this kind of datastructure has to provide to
//! be used by the framework. Such a class should
//! be implemented as a subclass of this class.
class ParticleContainer {
private:
	/* Copy operator private */
	ParticleContainer& operator=(const ParticleContainer&);

public:
	//! @brief The constructor
	//! @param bBoxMin coordinates of the lowest (in all coordinates) corner of the bounding box
	//! @param bBoxMax coordinates of the highest (in all coordinates) corner of the bounding box
	ParticleContainer(double bBoxMin[3], double bBoxMax[3]);

	//! @brief Default constructor
	ParticleContainer(){}
	//! @brief The destructor
	virtual ~ParticleContainer();

	virtual void readXML(XMLfileUnits& xmlconfig) = 0;

	//! @brief rebuild the datastructure
	//!
	//! Load-balancing decompositions change the position and size of the local region
	//! during runtime. Therefore, the datastructure needs to be rebuild completely.
	//! This method basically does what the constructor does as well, with the difference,
	//! that there are already particles stored, and particles which don't belong to the
	//! new region have to be deleted after rebuild
	virtual void rebuild(double bBoxMin[3], double bBoxMax[3]);

	//! @brief do necessary updates resulting from changed particle positions
	//!
	//! For some implementations of the interface ParticleContainer, the place
	//! where particles are stored might e.g. depend on the spacial position of the
	//! particle. So when some externel method (e.g. Leap-Frog) changes the spacial
	//! position of a particle, the representation within the particleContainer becomes
	//! invalid. This method restores a valid representation.
	virtual void update() = 0;

	//! @brief add a single Molecules to the ParticleContainer.
	//!
	//! It is important, that the Particle is entered in "the front" of the container.
	//! This is important when the container is traversed with the "next" method.
	//! E.g. a method traversing the container which adds copies of particles
	//! (periodic boundary) must not run over the added copies.
	//! This method has to be implemented in derived classes
	//! @param particle reference to the particle which has to be added
	virtual void addParticle(Molecule& particle) = 0;

	//! @brief traverse pairs which are close to each other
	//!
	//! Only interactions between particles which have a distance which is not
	//! larger than the cutoff radius are to be considered. \n
	//! This method has to be implemented in derived classes
	//! Precondition: All Particles of the process + halo molecules are stored
	//! Task: Run over all pairs (Each pair exactely once!) of particles (within cutoffradius)
	//! Important: Some pairs might be "duplicated": All pairs which cross the boundary occur twice
	//! (second time at the periodic image). Those pairs are from the point of view of the datastructure
	//! two different pairs, but they both times connect the same particles.
	//! For a pair which occurs twice, it has to be made sure, that one gets the status "original pair"
	//! and the other one "duplicated pair".
	//! For each pair found, there is an action executed, but it is a different action for
	//! original and duplicated pairs. Details about how to handle pairs can be found
	//! in the documentation for the class ParticlePairsHandler
	//! @param particlePairsHandler specified concrete action to be done for each pair
//	virtual void traversePairs(ParticlePairsHandler* particlePairsHandler) = 0;

	virtual void traverseCells(CellProcessor& cellProcessor) = 0;

	//! @return the number of particles stored in this container
	//!
	//! This number may includes particles which are outside of
	//! the bounding box
	virtual unsigned long getNumberOfParticles() = 0;

	//! @brief returns one coordinate of the lower corner of the bounding box
	//!
	//! @param dimension the coordinate which should be returned
	double getBoundingBoxMin(int dimension) const;


	double getHaloWidthNumCells();
	//! @brief returns one coordinate of the higher corner of the bounding box
	//!
	//! @param dimension the coordinate which should be returned
	double getBoundingBoxMax(int dimension) const;

	//! @brief Returns a pointer to the first particle in the Container
	virtual Molecule* begin() = 0;

	//! @brief Returns a pointer to the next particle in the Container
	//!
	//! The class internally has to store the Particle to which is currently pointed
	//! With the call of next, this internal pointer is advanced to the next particle
	//! and this new pointer is returned
	virtual Molecule* next() = 0;

	//! @brief Has to return the same as next() after it already pointed to the last particle
	virtual Molecule* end() = 0;

	virtual Molecule* deleteCurrent() = 0;

    /* TODO can we combine this with the update method? */
	//! @brief delete all Particles which are not within the bounding box
	virtual void deleteOuterParticles() = 0;

	//! @brief returns the width of the halo strip (for the given dimension index)
	//! @todo remove this method, because a halo_L shouldn't be necessary for every ParticleContainer
	//!       e.g. replace it by the cutoff-radius
	virtual double get_halo_L(int index) const = 0;

	//! @brief appends pointers to all particles in the halo region to the list
	virtual void getHaloParticles(std::list<Molecule*> &haloParticlePtrs) = 0;

	//! @brief fills the given list with pointers to all particles in the given region
	//! @param lowCorner minimum x-, y- and z-coordinate of the region
	//! @param highwCorner maximum x-, y- and z-coordinate of the region
	virtual void getRegion(double lowCorner[3], double highCorner[3], std::list<Molecule*> &particlePtrs) = 0;

	virtual double getCutoff() = 0;

    /* TODO: This goes into the component class */
	//! @brief counts all particles inside the bounding box
	virtual unsigned countParticles(unsigned int cid) = 0;

    /* TODO: This goes into the component class */
	//! @brief counts particles in the intersection of bounding box and control volume
	virtual unsigned countParticles(unsigned int cid, double* cbottom, double* ctop) = 0;

    /* TODO: Have a look on this */
	virtual void deleteMolecule(unsigned long molid, double x, double y, double z) = 0;

    /* TODO goes into grand canonical ensemble */
	virtual double getEnergy(ParticlePairsHandler* particlePairsHandler, Molecule* m1, CellProcessor& cellProcessor) = 0;
	virtual int localGrandcanonicalBalance() = 0;
	virtual int grandcanonicalBalance(DomainDecompBase* comm) = 0;
	virtual void grandcanonicalStep(ChemicalPotential* mu, double T, Domain* domain, CellProcessor& cellProcessor) = 0;

	//! @brief Update the caches of the molecules.
	void updateMoleculeCaches();

protected:

	//!  coordinates of the left, lower, front corner of the bounding box
	double _boundingBoxMin[3];
	//! coordinates of the right, upper, back corner of the bounding box
	double _boundingBoxMax[3];

};

#endif /* PARTICLECONTAINER_H_ */
