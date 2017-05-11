
//Calculation of the surface tension in a system with planar interfaces needs a Long Range Correction.
//
//The correction terms are based on Janecek (2006) and Lustig (1988).

#ifndef PLANAR_H_
#define PLANAR_H_

#include "LongRangeCorrection.h"

#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <cstdint>

#include "molecules/MoleculeForwardDeclaration.h"
class Domain;
class ParticleContainer;

class Planar:public LongRangeCorrection {
public:
	Planar(double cutoffT,double cutoffLJ,Domain* domain,  DomainDecompBase* domainDecomposition, ParticleContainer* particleContainer, unsigned slabs, Simulation _simulation);
	virtual ~Planar() {}

	virtual void init();
	virtual void readXML(XMLfileUnits& xmlconfig);
	virtual void calculateLongRange();
	double lrcLJ(Molecule* mol);
	// For non-equilibrium simulations the density profile must not be smoothed, therefore the density profile from the actual time step is used.
	void directDensityProfile();
	void SetSmoothDensityProfileOption(bool bVal) {_smooth = bVal;}
	virtual void writeProfiles(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep);

private:

	void centerCenter(double sig,double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj); 
	void centerSite(double sig,double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj);
	void siteSite(double sig,double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj);
	void dipoleDipole(unsigned ci,unsigned cj,unsigned si,unsigned sj);

	unsigned _slabs;
	unsigned numComp;
	unsigned *numLJ;
	unsigned *numDipole;
	unsigned numLJSum;
	unsigned numDipoleSum;
	unsigned *numLJSum2;
	unsigned *numDipoleSum2;
	bool _smooth;
//	bool _dipole;
	double *uLJ;
	double *vNLJ;
	double *vTLJ;
	double *fLJ;
	double *rho_g;
	double *rho_l;
	double *fDipole;
	double *uDipole;
	double *vNDipole;
	double *vTDipole;
	double *rhoDipole;
	double *rhoDipoleL;
	double *muSquare;
	double *eLong;
	double cutoff;
	double delta;
	unsigned cutoff_slabs;
	int frequency;
	double ymax;
	double boxlength[3];
	double V;
	int sint;
	double temp;
	unsigned simstep;
	
	ParticleContainer* _particleContainer;
	Domain* _domain;
	DomainDecompBase* _domainDecomposition;
	
	// write control
	uint64_t _nStartWritingProfiles;
	uint64_t _nWriteFreqProfiles;
	uint64_t _nStopWritingProfiles;
};


#endif /*Planar_H_*/
