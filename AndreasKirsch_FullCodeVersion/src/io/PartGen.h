#ifndef PARTGEN_H_
#define PARTGEN_H_

#include "io/InputBase.h"
#include <fstream>
#include <vector>
#include <list>

class ChemicalPotential;

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
class PartGen : public InputBase {

public:
	//! @brief reads in the configuration given in infile
	PartGen();

	//! @brief set the phase space file name
	void setPhaseSpaceFile(std::string filename);

	//! @brief set the phase space header file name (can be identical to the
	//         phase space file
	void setPhaseSpaceHeaderFile(std::string filename);

	//! @brief read the phase space components and header information
	void readPhaseSpaceHeader(Domain* domain, double timestep);

	//! @brief read the actual phase space information
	unsigned long readPhaseSpace(ParticleContainer* particleContainer, std::list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp);


	//! @brief stores the cluster file and related parameters
	//!
	//! Apart from storing gasDensity, fluidDensity and clusterFileName
	//! to the member variables, the values are used to calculate the
	//! new size of the simulation box
	void setClusterFile(double gasDensity, double fluidDensity,
			double volPercOfFluid, std::string clusterFileName);

	//! @brief creates a homogeneous distribution of particles
	//!
	//! Only particles inside of the specified bounding box (bBoxMin,
	//! bBoxMax) and belonging to this process (see domainDecomp) are
	//! created
	void createHomogeneousDist(ParticleContainer* particleContainer,
			std::vector<double> &bBoxMin, std::vector<double> &bBoxMax,
			Domain* domain, DomainDecompBase* domainDecomp);

	//! @brief creates an inhomogeneous distribution (gas with fluid drops)
	//!
	//! Only particles inside of the specified bounding box (bBoxMin,
	//! bBoxMax) and belonging to this process (see domainDecomp) are
	//! created
	void createClusters(ParticleContainer* particleContainer,
			std::vector<double> &bBoxMin, std::vector<double> &bBoxMax, Domain* domain,
			DomainDecompBase* domainDecomp);

	//! @brief prints some values of the config to stdout
	void printConfig();

	//###########################################################
	//### various get methods                                 ###
	//###########################################################

	double getTemperature();
	double getLSimBox(int dim);
	double getEpsilonRF();

	int getNumComps();
	int getNumSites(int comp);
	int getNumDipoles(int comp);
	int getNumQuadrupoles(int comp);

	int getNumCharges(int comp);
	int getNumTersoff(int comp);

	double getSitePos(int comp, int site, int dim);
	double getDipolePos(int comp, int site, int dim);
	double getQuadrupolePos(int comp, int site, int dim);

	double getChargePos(int comp, int site, int dim);
	double getTersoffPos(int comp, int site, int dim);

	double getEpsilonSite(int comp, int site);
	bool getShiftSite(int comp, int site);
	double getMSite(int comp, int site);
	double getSigmaSite(int comp, int site);
	double getEMyBody(int comp, int site, int dim);
	double getAbsMy(int comp, int site);
	double getEQBody(int comp, int site, int dim);
	double getAbsQ(int comp, int site);

	double getChargeMass(int comp, int site);
	double getCharge(int comp, int site);
	double getTersoffMass(int comp, int site);
	double getTersoffA(int comp, int site);
	double getTersoffB(int comp, int site);
	double getTersoff_lambda(int comp, int site);
	double getTersoff_mu(int comp, int site);
	double getTersoffR(int comp, int site);
	double getTersoffS(int comp, int site);
	double getTersoff_c(int comp, int site);
	double getTersoff_d(int comp, int site);
	double getTersoff_h(int comp, int site);
	double getTersoff_n(int comp, int site);
	double getTersoff_beta(int comp, int site);

	double getIDummy(int comp, int dim);
	double getEta(int comp1, int comp2);
	double getXi(int comp1, int comp2);

private:

	//! @brief creates a Molecule and adds it to the particleContainer
	//!
	//! Initial velocity, orientation and angular velocity are created
	//! randomly. The Component id is also created randomly but considers
	//! the ratio given in the config file for the probability
	void addParticle(int id, double x, double y, double z,
			ParticleContainer* particleContainer,
			Domain* domain, DomainDecompBase* domainDecomp);

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
	//### methods for processing lines or parts of lines      ###
	//### from the ITT-style input file                       ###
	//###########################################################

	//! @brief ignore numLines Lines in the given stream
	void ignoreLines(std::ifstream &inpfStream, int numLines);

	//! @brief get the next parameter value (int) from the input stream
	//!
	//! The next tokens in the steam must have the following format:
	//! <something without ws><some whitespace>=<some whitespace><a integer value>
	//! <whatever>
	//! The method returns the integer value, the stream is updated to point to the
	//! beginning of the next line
	int getIntParamValue(std::ifstream &inpfStream);

	//! @brief get the next parameter value (double) from the input stream
	//!
	//! For an explanation, see the documentation of getIntParamValue(...)
	double getDoubleParamValue(std::ifstream &inpfStream);

	//! @brief finds the next " = " in the stream and jumps at the place after it
	//!
	//! The "=" has to be surrounded by whitespace, the input stream will point
	//! to the next token after the "="
	void removePrefix(std::ifstream &inpfStream);

	//! @brief get 2 double values from the stream
	//!
	//! The stream has to point to the first value and the next value has
	//! to be separated from the first value (and everything after) by whitespace
	void getDoubleParamValues(std::ifstream &inpfStream, double &val1, double &val2);

	//! @brief get 3 double values from the stream (see 2-value-version)
	void getDoubleParamValues(std::ifstream &inpfStream, double &val1, double &val2,
			double &val3);

	//! @brief get 4 double values from the stream (see 2-value-version)
	void getDoubleParamValues(std::ifstream &inpfStream, double &val1, double &val2,
			double &val3, double &val4);

	//#############################################################
	//### method to read sections from the ITT-style input file ###
	//#############################################################

	//! @brief read in eps/kB_ref, mRef and sigmaRef
	//!
	//! The first two lines of the input stream are ignored, the next
	//! three lines must contain the three values
	void readRefValues(std::ifstream &inpfStream);

	//! @brief read in state point
	//!
	//! The first three lines of the input stream are ignored, the next
	//! line must contain the number of components, the next numcomp lines
	//! must contain the number of molecules per comp,
	//! the following two lines must contain rho and T
	void readStatePoint(std::ifstream &inpfStream);


	//! @brief read in information about time steps and cutoff radius
	void readAlgorithm(std::ifstream &inpfStream, std::vector<double> &simBoxRatio);

	//###########################################################
	//### methods for operations on matrices and vectors      ###
	//###########################################################

	//! @brief init the 3D matrix to the given size and with the given value
	void set3DMatrixValues(std::vector<std::vector<std::vector<double> > > &matrix, int numx,
			int numy, int numz, double value);

	//! @brief init the matrix to the given size and with the given value
	void setMatrixValues(std::vector<std::vector<double> > &matrix, int numrows,
			int numcols, double value);

	//! @brief init the int-vector to the given size and with the given value
	void setVectorValues(std::vector<int> &vect, int num, int value);

	//! @brief init the double-vector to the given size and with the given value
	void setVectorValues(std::vector<double> &vect, int num, double value);

	//! @brief return the scalar produkt of two vectors
	double dotprod(std::vector<double> &v1, std::vector<double> &v2);

	//! @brief multiply the given vector with the given value
	void vecmult(std::vector<double> &v, double value);

	//! @brief return the maximum of the three values
	double max(double a, double b, double c);

	//! @brief add the second vector to the first (only the first is modified)
	void vecadd(std::vector<double> &v1, std::vector<double> &v2);

	//! @brief multiply the given matrix with the given value
	void matmult(std::vector<std::vector<double> > &m, double value);

	//! @brief return the product of the two given matrices
	std::vector<std::vector<double> > matmult(std::vector<std::vector<double> > &m1,
			std::vector<std::vector<double> > &m2);

	//! @brief return the product of the row vector v with the matrix m
	void matmult(std::vector<double> &v, std::vector<std::vector<double> > &m);

	//! @brief multiply the given 3D-matrix with the given value
	void mat3dmult(std::vector<std::vector<std::vector<double> > > &m, double value);

	//! @brief lvst ein 3x3 LGS mit rechter Seite Null
	//! @todo enhance comments
	void solveLGS(std::vector<std::vector<double> > &A, std::vector<double> &x);

	//! @brief calcuate the eigenvectors of the given Matrix (3x3 symmetric matrix)
	//! @todo enhance comments
	void getEigenvecs(std::vector<std::vector<double> > &m,
			std::vector<std::vector<double> > &eigenvecs);

	//! @brief transforms the molecule to a principle axis system
	void principleAxisTransform();

	//! @brief erstellt die Massentraegheitsmatrix
	//! @todo enhance comments
	void createTraegheitsMatrix(std::vector<std::vector<double> > &matrix,
			std::vector<std::vector<double> > &sitesPos,
			std::vector<double> &masses, int numsites);

	//###########################################################
	//### methods which create random numbers                 ###
	//###########################################################

	//! @brief returns a random int number between a and b (incl. both)
	int randint(int a, int b);

	//! @brief returns a or b, each with 50% probability
	double randchoice(double a, double b);

	//! @brief returns a random double value between a and b
	double randdouble(double a, double b);

	//! @brief returns a random component id
	//!
	//! The probability for each component is proportional
	//! to the number of molecules with that id
	int randCompID();


	//###########################################################
	//### member variables                                    ###
	//###########################################################

	// physical constants
	//! @brief atomic mass, not reduced, [atomic mass] = kg
	static const double _atomicMassDim;
	//! @brief Avogadro constant, not reduced, [Avogadro] = 1/mol
	static const double _avogadroDim;
	//! @brief Boltzmann constant, not reduced,  [Boltzmann] = J/K
	//! @note note that for simulation purposes, k equals 1 instead of 1.3e-23.
	static const double _boltzmannDim;
	//! @brief Coulomb constant in units of [Coco] = Jm/e^2
	static const double _cocoDim;
	//! @brief Debye in units of [my] = me
	static const double _debyeDim;

	//! @brief density of the gas phase
	double _gasDensity;
	//! @brief density of the fluid phase (drops)
	double _fluidDensity;
	//! @brief file in which the fluid drops are defined (lines with x,y,z,r)
	std::string _clusterFile;
	//! @brief each element is a sphere (vector containing x,y,z and r)
	std::vector<std::vector<double> > _localClusters;

	// Reference Values (units not reduced!)
	double _epsilonRefDim;
	double _mRefDim;
	double _sigmaRefDim;

	//
	//! @brief number of different components
	int _numberOfComponents;
	//! @brief total number of molecules
	int _numberOfMolecules;
	//! @brief number of molecules from each molecule type
	std::vector<int> _numMolsPerComp;

	//! @brief density of a homogeneous distribution (reduced unit)
	double _rho;
	//! @brief temperature of the simulated material (reduced unit)
	double _temperature;
	//! @brief describes the relative length of the simulation box
	std::vector<double> _simBoxRatio;
	//! @brief length of the global simulation box
	std::vector<double> _simBoxLength;

	//! @brief epsilon reaction field
	double _epsilonRF;

	//! @brief each element of the vector is the mass of a component
	std::vector<double> _massOfComps;
	//! @brief each element is the number of DOF of a component
	std::vector<int> _degreesOfFreedom;

	// Mixing parameters for interactions between different components
	std::vector<std::vector<double> > _etaLB;
	std::vector<std::vector<double> > _xiLB;

	// Variables used to represent the components
	std::vector<int> _numSites;
	std::vector<int> _numDipoles;
	std::vector<int> _numQuadrupoles;

	std::vector<int> _numCharges;
	std::vector<int> _numTersoff;

	std::vector< std::vector<double> > _iBodyDummy;
	std::vector< std::vector< std::vector<double> > > _rSiteBody;
	std::vector< std::vector< std::vector<double> > > _rDipoleBody;
	std::vector< std::vector< std::vector<double> > > _rQuadrupoleBody;
	std::vector< std::vector<double> > _epsilonSite;
	std::vector< std::vector<double> > _mSite;
	std::vector< std::vector<double> > _sigmaSite;
	std::vector< std::vector<bool> > _shiftSite;
	std::vector< std::vector< std::vector<double> > > _eMyBody;
	std::vector< std::vector<double> > _absMy;
	std::vector< std::vector< std::vector<double> > > _eQBody;
	std::vector< std::vector<double> > _absQ;

	std::vector< std::vector< std::vector<double> > > _rChargeBody;
	std::vector< std::vector< std::vector<double> > > _rTersoffBody;
	std::vector< std::vector<double> > _mCharge;
	std::vector< std::vector<double> > _qCharge;
	std::vector< std::vector<double> > _mTersoff;
	std::vector< std::vector<double> > _ATersoff;
	std::vector< std::vector<double> > _BTersoff;
	std::vector< std::vector<double> > _lambdaTersoff;
	std::vector< std::vector<double> > _muTersoff;
	std::vector< std::vector<double> > _RTersoff;
	std::vector< std::vector<double> > _STersoff;
	std::vector< std::vector<double> > _cTersoff;
	std::vector< std::vector<double> > _dTersoff;
	std::vector< std::vector<double> > _hTersoff;
	std::vector< std::vector<double> > _nTersoff;
	std::vector< std::vector<double> > _betaTersoff;

	std::vector< std::vector<double> > _iBody;
	std::vector< std::vector<double> > _invIBody;

	// double _cutoffRadius;
	double _LJCutoffRadius;
	// double _tersoffCutoffRadius;
};

#endif /*PARTGEN_H_*/
