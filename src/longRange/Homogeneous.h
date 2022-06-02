
#ifndef HOMOGENEOUS_H__
#define HOMOGENEOUS_H__

#include "LongRangeCorrection.h"
#include "molecules/Component.h"
#include "molecules/Comp2Param.h"
#include <cmath>

class Simulation;
class Domain;
class LongRangeCorrection;

class Homogeneous: public LongRangeCorrection{

public:
	Homogeneous(double cutoffRadius, double cutoffRadiusLJ,  Domain* domain, ParticleContainer* particleContainer, Simulation* simulation);
	~Homogeneous() override = default;

	void init() override;
	void readXML(XMLfileUnits& xmlconfig) override {}
	void calculateLongRange() override;
	void writeProfiles(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) override {}

	// Get potential energy correction per molecule
	double getUpotCorr(Molecule* /* mol */) override;

private:
	/* TODO: Comments on all the functions */
	// Long range correction for the Lennard-Jones interactions based on Lustig (1988)
	double _TICCu(int n,double rc,double sigma2);
	double _TICSu(int n,double rc,double sigma2,double tau);
	double _TISSu(int n,double rc,double sigma2,double tau1,double tau2);
	double _TICCv(int n,double rc,double sigma2);
	double _TICSv(int n,double rc,double sigma2,double tau);
	double _TISSv(int n,double rc,double sigma2,double tau1,double tau2);
	
	//! Components resp. molecule types
	std::vector<Component>* _components{nullptr};
	//! parameter streams for each possible pair of molecule-types
	Comp2Param _comp2params;
	
	Domain* _domain{nullptr};
	ParticleContainer* _particleContainer{nullptr};

	double _uPotCorrPerMol; // Correction of potential energy per particle

	double _cutoff{std::numeric_limits<double>::quiet_NaN()};
	double _cutoffLJ{std::numeric_limits<double>::quiet_NaN()};

	double _upotCorrLJ_no_num_molecules{std::numeric_limits<double>::quiet_NaN()};
	double _virialCorrLJ_no_num_molecules{std::numeric_limits<double>::quiet_NaN()};
	double _mySelbstTerm_no_num_molecules{std::numeric_limits<double>::quiet_NaN()};
};

#endif /* __HOMOGENEOUS_H__ */
