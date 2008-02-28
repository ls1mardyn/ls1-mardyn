#ifndef PARTICLECONTAINER_H_
#define PARTICLECONTAINER_H_

#include <vector>
#include <list>

#include "utils/Log.h"

namespace datastructures {

  template<class ParticleType>
  class ParticlePairsHandler;
  
  template<class ParticleType>
  class ParticleContainer; 

}

namespace parallel {
  class DomainDecompBase;
}

using namespace std;

//! @brief This Interface is used to get access to particles and pairs of particles
//! @author Martin Buchholz
//! 
//! A particleContainer is used to store Particles with a short-range potential
//! in a way that the access to pairs of neighbouring particles is efficient.
//! Neighbouring particles are particles which have a distance shorter than
//! a given cutoff radius. 
//! A common task when using a PariticleContainer is to do do something for
//! all particles. This can be done using the methods begin, next and end, e.g.:
//! \code
//! ParticleContainer* partContPtr; 
//! ParticleType* partPtr; 
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
//! An exception to this is when particles are moved in a time step. Is has to be
//! ensured, that particles which leave the bounding box are properly handled.
//! 
//! For non-cuboid regions, the bounding box still has to be defined as it gives
//! an approximation for the region that is covered by the ParticleContainer.
//! 
//! This interface doesn't implement the datastructure, it just tells which
//! methods a class implementing this kind of datastructure has to provide to 
//! be used by the framework. Such a class should
//! be implemented as a subclass of this class.
template<class ParticleType>
class datastructures::ParticleContainer {
  public:
    //! @brief The constructor
    //! @param partPairsHandler specified concrete action to be done for each pair
    //! @param bBoxMin coordinates of the lowest (in all coordinates) corner of the bounding box
    //! @param bBoxMax coordinates of the highest (in all coordinates) corner of the bounding box
    ParticleContainer(datastructures::ParticlePairsHandler<ParticleType>& partPairsHandler, double bBoxMin[3], double bBoxMax[3]);
    
    //! The destructor
    virtual ~ParticleContainer();
    
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
    virtual void addParticle(ParticleType& particle) = 0;
  
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
    virtual void traversePairs() = 0;
    
    //! @return the number of particles stored in this container
    //!
    //! This number may includes particles which are outside of
    //! the bounding box
    virtual unsigned long getNumberOfParticles() = 0;
    
    //! @brief returns one coordinate of the lower corner of the bounding box
    //!
    //! @param dimension the coordinate which should be returned
    double getBoundingBoxMin(int dimension);

    //! @brief returns one coordinate of the higher corner of the bounding box
    //!
    //! @param dimension the coordinate which should be returned
    double getBoundingBoxMax(int dimension);
    
    //! @brief Returns a pointer to the first particle in the Container
    virtual ParticleType* begin() = 0;

    //! @brief Returns a pointer to the next particle in the Container
    //!
    //! The class internally has to store the Particle to which is currently pointed
    //! With the call of next, this internal pointer is advanced to the next particle
    //! and this new pointer is returned
    virtual ParticleType* next() = 0;

    //! @brief Has to return the same as next() after it already pointed to the last particle
    virtual ParticleType* end() = 0;
    
    //! @brief delete all Particles which are not within the bounding box
    virtual void deleteOuterParticles() = 0;
    
    //! @brief returns the width of the halo strip (for the given dimension index)
    //! @todo remove this method, because a halo_L shouldn't be necessary for every ParticleContainer
    virtual double get_halo_L(int index);
    
  protected:

    //! A ParticlePairsHandler is used to process pairs of particles
    datastructures::ParticlePairsHandler<ParticleType>& _particlePairsHandler;
    
    //!  coordinates of the left, lower, front corner of the bounding box
    double _boundingBoxMin[3];
    //! coordinates of the right, upper, back corner of the bounding box
    double _boundingBoxMax[3]; 
           
  private:
    //! Logging interface
    static utils::Log _log;
    

    
};

#include "datastructures/ParticleContainer.cpph"

#endif /*PARTICLECONTAINER_H_*/
