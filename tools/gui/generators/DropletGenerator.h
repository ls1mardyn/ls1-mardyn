/*
 * DropletGenerator.h
 *
 *  Created on: May 14, 2011
 *      Author: kovacevt
 */

class Domain;

#include "MDGenerator.h"
#include "common/ComponentParameters.h"
#include "parallel/DomainDecompBase.h"

#ifndef DROPLETGENERATOR_H_
#define DROPLETGENERATOR_H_

class DropletGenerator: public MDGenerator {

private:
	double _temperature;
	long int numOfMolecules;
	std::vector<Component>& _components;

	double fluidDensity;
	double rho;
	double gasDensity;
	double fluidVolume;
	double maxSphereVolume;
	double numSphereSizes;
	double simBoxLength[3];

	//! @brief each element is a sphere (vector containing x,y,z and r)
	std::vector<std::vector<double> > localClusters;

public:
	DropletGenerator();
	virtual ~DropletGenerator();

	/**
	 *
	 * Apart from storing gasDensity, fluidDensity and clusterFileName
	 * to the member variables, the values are used to calculate the
	 * new size of the simulation box
	 * @param gasDensity the density of the gas
	 * @param fluidDensity the density of the fluid
	 * @param fluidVolumePercent determines what percentage of the volume is covered
	 *        by fluid. Given in percent(so int would be more appropriate!).
	 * @param numSphereSizes parameter for DropletGenerator.
	 *
	 * @TODO correct data types
	 */
	void setClusterParameters(double gasDensity, double fluidDensity,
			double fluidVolumePercent, double maxSphereVolume,
			int numSphereSizes);

	/**
	 * @brief creates an inhomogeneous distribution (gas with fluid drops)
	 *
	 * Only particles inside of the specified bounding box (bBoxMin,
	 * bBoxMax) and belonging to this process (see domainDecomp) are
	 * created
	 *
	 * Generates dropplets, called from GenerateDrawableMolecules()
	 *
	 * @return the highest id of a molecule generated.
	 */
	unsigned long generateMoleculesCluster(ParticleContainer* particleContainer,
			std::vector<double> &bBoxMin, std::vector<double> &bBoxMax, Domain* domain,
			DomainDecompBase* domainDecomp);

	//Domain* domain,DomainDecompBase* domainDecomp);

	//! @brief creates those spheres that are close to this domain
	//!
	//! All spheres from the cluster file are read. Then all 27 periodic
	//! copies are created and for each of them it is checked, whether the
	//! sphere intersects this process' domain. If that's the case, the
	//! sphere (the corresponding periodic copy) is added to _localClusters

	void readLocalClusters(Domain* domain, DomainDecompBase* domainDecomp);

	//! @brief returns true if the given position is part of a previous sphere
	//!
	//! When adding the particles for all spheres, some spheres might intersect
	//! At a given position (and close to it), only one particle is allowed to
	//! exist (otherwise, infinitive forces will occur). So when adding the particles
	//! for one sphere, all previous spheres (lower clusterid) must be checked whether
	//! they contain the position given by (x,y,z).
	//! This method checks all spheres with id < clusterid whether they possess (x,y,z),
	//! and returns true if (at least) one of them does, otherwise it returns false
	bool belongsToPreviousCluster(double x, double y, double z, int clusterid);

	//! @brief returns true if the given position is close to any sphere
	//!
	//! belongsToPreviousCluster is used when adding a new spheres to check the
	//! previous spheres. This method is useful when some other particle is added,
	//! which is not part of a cluster (so not part of the fluid phase). Then, not
	//! only has to be checked, whether any of the spheres contains this position,
	//! but also wheter it is close (closer than offset) to one of the spheres.
	//! This is because the particles from  the fluid phase and the gas phase are
	//! created on two different grids. In one of the phases, if two points are not
	//! equal, the are guaranteed to have a certain distance given by the grid. But
	//! two points from different grids (different phases) can be arbitrarily close
	//! together. So a new added particle, which is not part of a sphere but very
	//! close to it, could also cause almost infinitive forces. This method returns
	//! true if the distance of the position (x,y,z) to any of the spheres is
	//! less (or equal) than offset. Otherwise it returns false.
	bool closeToAnyCluster(double x, double y, double z, double offset);

	//! @brief counts the number of molecules of each component type.
	//!
	//! This method is usually only needed once in the beginning of the simulation and only
	//! if the particles were not read in from a single file but read in from one file per proc or
	//! if the particles were created by each proc seperately.
	//! @param moleculeContainer container for the molecules
	//! @param compCount vector which has to have the size which equals the number of components
	//!                  this method will will the vector with the number of molecules for each
	//!                  of the components (in the global domain)
	//! @return the number of molecules in the global domain is returned
	//!
	//! @TODO move this method to the DropletGenerator!
	unsigned long countMolecules(DomainDecompBase* domainDecomp, ParticleContainer* moleculeContainer, std::vector<unsigned long> &compCount);


	//! @brief returns a guaranteed distance of (x,y,z) to the local domain
	//!
	//! This method is e.g. used by a particle generator which creates clusters (nuclei, drops).
	//! Only if the cluster is close (cluster radius larger then guaranteedDistance) to the
	//! domain the particles have to be created.
	//! @param x x-coordinate of the position for which the guaranteed distance is returned
	//! @param y y-coordinate of the position for which the guaranteed distance is returned
	//! @param z z-coordinate of the position for which the guaranteed distance is returned
	//! @param domain might be needed to get the bounding box
	//!
	double guaranteedDistance(double x, double y, double z, DomainDecompBase* domainDecomp, Domain* domain);


	virtual void readPhaseSpaceHeader(Domain* domain, double timestep);

	//! @brief read the phase space components and header information
	unsigned long readPhaseSpace(ParticleContainer* particleContainer,
			Domain* domain,
			DomainDecompBase* domainDecomp);

	std::vector<ParameterCollection*> getParameters();

	//void generatePreview();

	void setParameter(Parameter* p);

	virtual bool validateParameters();

	/**
	 * set the orientation of a molecule according to their position on
	 * an alpha-fcc-lattice (code copied from Martin Buchholz's branch - PartGen.cpp)
	 * as described in "Molecular Simulation of Liquids".
	 */
	void getFCCOrientation(int q_type, double q[4]);

};

#endif /* DROPLETGENERATOR_H_ */
