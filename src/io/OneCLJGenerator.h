#ifndef GEN1CLJ_H_
#define GEN1CLJ_H_

#include "io/InputBase.h"
#include <fstream>
#include <vector>

/**
 * @author Martin Buchholz, Wolfgang Eckhardt
 */

using namespace std;


//! @brief generates homogeneous and inhomogeneous particle distributions
//! 
//! This class is intended to be used in parallel environments (even though
//! it does not contain any parallelisation itself). With this class,
//! each process is able to create particles in the own domain. Either
//! homogeneous or inhomogeneous distribution can be created.
//! For inhomogeneous distributions, a cluster file has to specified which
//! contains a list of spheres (given by x, y, z and radius). In those
//! spheres, fluid (high density) is created, in the rest of the domain,
//! gas (low density) is created.
class OneCLJGenerator : public InputBase{

public:
	//! epsilon parameter of the 1CLJ center
	static const double eps;
	//! sigma parameter of the 1CLJ center
	static const double sigma;
	//! mass of the 1CLJ center
	static const double mass;

	/**! Initialize the corresponding member variables and calculates the size of
	 *   the simulation box from rho and N
	 *   @param mode must be "Homogeneous" or "Cluster"
	 *   @param N the number of molecules
	 *   @param rho the density
	 *   @param T the temperature
	 */
	OneCLJGenerator(string mode, unsigned long N, double T);


	//! @brief stores the cluster file and related parameters
	//!
	//! Apart from storing gasDensity, fluidDensity and clusterFileName
	//! to the member variables, the values are used to calculate the
	//! new size of the simulation box
	//!
	//! @param gasDensity the density of the gas
	//! @param fluidDensity the density of the fluid
	//! @param fluidVolumePercent determines what percentage of the volume is covered
	//!        by fluid. Given in percent(so int would be more appropriate!).
	//! @param numSphereSizes parameter for DropletGenerator.
	//!
	//! @TODO correct data types
	//!
	void setClusterParameters(double gasDensity, double fluidDensity, double fluidVolumePercent,
	        double maxSphereVolume, double numSphereSizes);

	/**! @brief store the density
	 *   Moreover, calculate the new size of the simulation box.
	 */
	void setHomogeneuosParameter(double rho);

	//! @brief NOP for this class
	void setPhaseSpaceFile(string filename) {
	}

	//! @brief NOP for this class
	void setPhaseSpaceHeaderFile(string filename) {
	}

	//! @brief read the phase space components and header information
	void readPhaseSpaceHeader(Domain* domain, double timestepLength);

	//! @brief read the actual phase space information
	virtual unsigned long readPhaseSpace(ParticleContainer* particleContainer, std::list<ChemicalPotential>* lmu,
	        Domain* domain, DomainDecompBase* domainDecomp);



	//! @brief creates a homogeneous distribution of particles
	//!
	//! Only particles inside of the specified bounding box (bBoxMin,
	//! bBoxMax) and belonging to this process (see domainDecomp) are
	//! created
	void createHomogeneousDist(ParticleContainer* particleContainer, vector<double> &bBoxMin, vector<double> &bBoxMax,
	        Domain* domain, DomainDecompBase* domainDecomp);

	//! @brief creates an inhomogeneous distribution (gas with fluid drops)
	//!
	//! Only particles inside of the specified bounding box (bBoxMin,
	//! bBoxMax) and belonging to this process (see domainDecomp) are
	//! created
	void createClusters(ParticleContainer* particleContainer, vector<double> &bBoxMin, vector<double> &bBoxMax,
	        Domain* domain, DomainDecompBase* domainDecomp);


 private:

	//! @brief creates a Molecule and adds it to the particleContainer
	//!
	//! Initial velocity, orientation and angular velocity are created
	//! randomly. The Component id is also created randomly but considers
	//! the ratio given in the config file for the probability
	void addParticle(int id, double x, double y, double z, ParticleContainer* particleContainer, Domain* domain,
	        DomainDecompBase* domainDecomp);

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


//###########################################################
//### methods which create random numbers                 ###
//###########################################################


	//! @brief returns a random double value between a and b
	double randdouble(double a, double b);


//###########################################################
//### member variables                                    ###
//###########################################################

	//! @brief density of the gas phase
	double _gasDensity;
	//! @brief density of the fluid phase (drops)
	double _fluidDensity;
	//! @brief file in which the fluid drops are defined (lines with x,y,z,r)
	//string _clusterFile;
	double _fluidVolume;
	double _maxSphereVolume;
	double _numSphereSizes;

	//! @brief each element is a sphere (vector containing x,y,z and r)
	vector<vector<double> > _localClusters;

	//! @brief total number of molecules
	int _numberOfMolecules;

	string _mode;
	//! @brief density of a homogeneous distribution (reduced unit)
	double _rho;
	//! @brief temperature of the simulated material (reduced unit)
	double _temperature;
	//! @brief length of the global simulation box
	vector<double> _simBoxLength;

	unsigned long int _moleculeCountOffset;

};

#endif /*Gen1CLJ_H_*/
