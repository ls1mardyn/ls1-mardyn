/*
 * ParticleContainerFactory.h
 *
 * @Date: 21.09.2010
 * @Author: eckhardw
 */

#ifndef PARTICLECONTAINERFACTORY_H_
#define PARTICLECONTAINERFACTORY_H_

#include <string>

class ParticleContainer;
class Domain;
class DomainDecompBase;

/**
 * This is a factory class to create different types of particleContainers,
 * initialized and ready for testing.
 */
class ParticleContainerFactory {

public:

	enum Type {LinkedCell};

	/**
	 * creates a small, empty particleContainer. For the moment, only LinkedCell is supported.
	 * The caller is responsible for deleting the pointer!
	 */
	static ParticleContainer* createEmptyParticleContainer(Type type);

	/**
	 * return a (more or less - rather less ;) initialized Particlecontainer. Actually, only the
	 * given file is read, and the domain is set up correspondingly to the file, and after that the
	 * container is initialized from the file. For the domain, the parameter streams are initialized,
	 * but no further initialization methods are called on any class!
	 *
	 * The caller is responsible for deleting the pointer!
	 */
	static ParticleContainer* createInitializedParticleContainer(
			Type type, Domain* domain, DomainDecompBase* domainDecomposition, double cutoff, const std::string& fileName, bool binary);

};

#endif /* PARTICLECONTAINERFACTORY_H_ */
