/*
 * LinkedCellsCUDAAdapter.h
 *
 *  Created on: Jun 8, 2011
 *      Author: andreas
 */

#ifndef LINKEDCELLSCUDAPROXY_H_
#define LINKEDCELLSCUDAPROXY_H_

#include "particleContainer/ParticleContainer.h"
#include "particleContainer/LinkedCells.h"
#include "particleContainer/Cell.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "Domain.h"

#include "behaviorProbe.h"
#include "benchmark.h"

#include "moleculeInteraction.h"

extern "C" {
	extern char _binary_cuda_kernel_ptx_start, _binary_cuda_kernel_ptx_end;
}

class LinkedCellsCUDADecorator : public ParticleContainer {
	// null handler
	class Handler : public ParticlePairsHandler {
	public:
		//! @brief things to be done before particle pairs are processed
		virtual void init() {}

		//! @brief things to be done after particle pairs are processed
		virtual void finish() {}

		//! @brief things to be done for each particle pair
		//!
		//! @param particle1 first particle
		//! @param particle2 second particle
		//! @param distanceVector[3] distance between the two particles
		//! @param pairType describes whether the pair is a original pair(0) or a duplicated pair(1)
		//!                 for details about pair types see comments on traversePairs() in ParticleContainer
		virtual double processPair(Molecule& particle1, Molecule& particle2, double distanceVector[3], int pairType, double dd, bool calculateLJ) {
			assert( false );
			return 0.0f;
		}
		virtual void preprocessTersoffPair(Molecule& particle1, Molecule& particle2, bool pairType) {
			assert( false );
		}
		virtual void processTersoffAtom(Molecule& particle1, double params[15], double delta_r) {
			assert( false );
		}

		virtual void recordRDF() {
			//assert( false );
		}
	};

public:
	//! @brief The constructor
	//! @param partPairsHandler specified concrete action to be done for each pair
	//! @param bBoxMin coordinates of the lowest (in all coordinates) corner of the bounding box
	//! @param bBoxMax coordinates of the highest (in all coordinates) corner of the bounding box
	LinkedCellsCUDADecorator(Domain &domain, double bBoxMin[3], double bBoxMax[3], double cutoffRadius,
			 ParticlePairsHandler* partPairsHandler)
	: ParticleContainer(partPairsHandler, bBoxMin, bBoxMax),
	  _domain( domain ),
	  _linkedCells(bBoxMin, bBoxMax, cutoffRadius, cutoffRadius, cutoffRadius, 1.0, partPairsHandler),
	  _moleculeInteraction( cuda().loadModuleData( &_binary_cuda_kernel_ptx_start, &_binary_cuda_kernel_ptx_end, MAX_REGISTER_COUNT ), _linkedCells, _domain )
	{}

	//! @brief The destructor
	virtual ~LinkedCellsCUDADecorator() {
	}

	//! @brief rebuild the datastructure
	//!
	//! Load-balancing decompositions change the position and size of the local region
	//! during runtime. Therefore, the datastructure needs to be rebuild completely.
	//! This method basically does what the constructor does as well, with the difference,
	//! that there are already particles stored, and particles which don't belong to the
	//! new region have to be deleted after rebuild
	virtual void rebuild(double bBoxMin[3], double bBoxMax[3]) {
		ParticleContainer::rebuild( bBoxMin, bBoxMax );
		_linkedCells.rebuild(bBoxMin, bBoxMax);
	}

	//! @brief do necessary updates resulting from changed particle positions
	//!
	//! For some implementations of the interface ParticleContainer, the place
	//! where particles are stored might e.g. depend on the spacial position of the
	//! particle. So when some externel method (e.g. Leap-Frog) changes the spacial
	//! position of a particle, the representation within the particleContainer becomes
	//! invalid. This method restores a valid representation.
	virtual void update() {
		_linkedCells.update();
	}

	//! @brief add a single Molecules to the ParticleContainer.
	//!
	//! It is important, that the Particle is entered in "the front" of the container.
	//! This is important when the container is traversed with the "next" method.
	//! E.g. a method traversing the container which adds copies of particles
	//! (periodic boundary) must not run over the added copies.
	//! This method has to be implemented in derived classes
	//! @param particle reference to the particle which has to be added
	virtual void addParticle(Molecule& particle) {
		_linkedCells.addParticle(particle);
	}

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
	virtual void traversePairs() {
#ifdef COMPARE_TO_CPU
		_linkedCells.traversePairs();

		double cpuPotential = _domain.getLocalUpot();
		double cpuVirial = _domain.getLocalVirial();
#endif
		double potential;
		double virial;
		_moleculeInteraction.calculate( potential, virial );

		BehaviorProbe::FMset();

		// update the domain values
		_domain.setLocalUpot( potential );
		_domain.setLocalVirial( virial );

#ifdef COMPARE_TO_CPU
		printf( "CPU Potential: %f CPU Virial: %f\n", cpuPotential, cpuVirial );
#endif
		printf( "CUDA Potential: %f CUDA Virial: %f\n", potential, virial );
	}

	//! @return the number of particles stored in this container
	//!
	//! This number may includes particles which are outside of
	//! the bounding box
	virtual unsigned long getNumberOfParticles() {
		return _linkedCells.getNumberOfParticles();
	}

	//! @brief Returns a pointer to the first particle in the Container
	virtual Molecule* begin() {
		return _linkedCells.begin();
	}

	//! @brief Returns a pointer to the next particle in the Container
	//!
	//! The class internally has to store the Particle to which is currently pointed
	//! With the call of next, this internal pointer is advanced to the next particle
	//! and this new pointer is returned
	virtual Molecule* next() {
		return _linkedCells.next();
	}

	//! @brief Has to return the same as next() after it already pointed to the last particle
	virtual Molecule* end() {
		return _linkedCells.end();
	}

	virtual Molecule* deleteCurrent() {
		return _linkedCells.deleteCurrent();
	}

	//! @brief delete all Particles which are not within the bounding box
	virtual void deleteOuterParticles() {
		_linkedCells.deleteOuterParticles();
	}

	//! @brief returns the width of the halo strip (for the given dimension index)
	//! @todo remove this method, because a halo_L shouldn't be necessary for every ParticleContainer
	//!       e.g. replace it by the cutoff-radius
	virtual double get_halo_L(int index) const {
		return _linkedCells.get_halo_L(index);
	}

	//! @brief appends pointers to all particles in the boundary region to the list
	virtual void getBoundaryParticles(std::list<Molecule*> &boundaryParticlePtrs) {
		_linkedCells.getBoundaryParticles(boundaryParticlePtrs);
	}

	//! @brief appends pointers to all particles in the halo region to the list
	virtual void getHaloParticles(std::list<Molecule*> &haloParticlePtrs) {
		_linkedCells.getHaloParticles(haloParticlePtrs);
	}

	//! @brief fills the given list with pointers to all particles in the given region
	//! @param lowCorner minimum x-, y- and z-coordinate of the region
	//! @param highwCorner maximum x-, y- and z-coordinate of the region
	virtual void getRegion(double lowCorner[3], double highCorner[3], std::list<Molecule*> &particlePtrs) {
		_linkedCells.getRegion(lowCorner, highCorner, particlePtrs);
	}

	virtual double getCutoff() {
		return _linkedCells.getCutoff();
	}
	virtual double getLJCutoff() {
		return _linkedCells.getLJCutoff();
	}

	virtual double getTersoffCutoff() {
		return _linkedCells.getTersoffCutoff();
	}

	virtual void countParticles(Domain* d) {
		_linkedCells.countParticles(d);
	}

	//! @brief counts all particles inside the bounding box
	virtual unsigned countParticles(int cid) {
		return _linkedCells.countParticles(cid);
	}

	//! @brief counts particles in the intersection of bounding box and control volume
	virtual unsigned countParticles(int cid, double* cbottom, double* ctop) {
		return _linkedCells.countParticles(cid, cbottom, ctop);
	}

	virtual void deleteMolecule(unsigned long molid, double x, double y, double z) {
		return _linkedCells.deleteMolecule(molid, x, y, z);
	}

	virtual double getEnergy(Molecule* m1) {
		return _linkedCells.getEnergy(m1);
	}

	virtual int localGrandcanonicalBalance() {
		return _linkedCells.localGrandcanonicalBalance();
	}

	virtual int grandcanonicalBalance(DomainDecompBase* comm) {
		return _linkedCells.grandcanonicalBalance(comm);
	}

	virtual void grandcanonicalStep(ChemicalPotential* mu, double T) {
		return _linkedCells.grandcanonicalStep(mu, T);
	}

	//! @brief find out whether m1 is before m2 (in some global ordering)
	//!
	//! At the boundary between two processes (if used in parallel mode), the forces
	//! for pairs which cross the boundary are calculated twice (once by each proc who
	//! owns one of the particles). But the contribution to macroscopic value must be
	//! counted only once, which is done by the process who owns the "first" particle.
	//! As order criterion, the spacial position is used in this method. The particles
	//! with lower x-coordinate is first (if equal, then y- or z-coordinate).
	//! For pairs which are completely on one process, the first particle can be
	//! determined from the cell structure. But for pairs on different procs, the
	//! corresponding cell discretisations might be different as well, and therefore
	//! the cell structure must not be used to determine the order.
	virtual bool isFirstParticle(Molecule& m1, Molecule& m2) {
		return _linkedCells.isFirstParticle(m1, m2);
	}

	//! @brief Update the caches of the molecules.
	void updateMoleculeCaches() {
		_linkedCells.updateMoleculeCaches();
	}

private:
	Domain &_domain;

	// internal particle container
	LinkedCells _linkedCells;
	MoleculeInteraction _moleculeInteraction;
};

#endif /* LINKEDCELLSCUDAADAPTER_H_ */
