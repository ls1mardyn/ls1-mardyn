/*
 * ParticleContainerFactory.h
 *
 * @Date: 21.09.2010
 * @Author: eckhardw
 */

#ifndef PARTICLECONTAINERFACTORY_H_
#define PARTICLECONTAINERFACTORY_H_

class ParticleContainer;

/**
 * This is a factory class to create different types of particleContainers,
 * initialized and ready for testing.
 */
class ParticleContainerFactory {

public:

	enum type {LinkedCell};

	/**
	 * creates a small, empty particleContainer. For the moment, only LinkedCell is supported.
	 * The caller is responsible for deleting the pointer!
	 */
	static ParticleContainer* createEmptyParticleContainer(type type);

};

#endif /* PARTICLECONTAINERFACTORY_H_ */
