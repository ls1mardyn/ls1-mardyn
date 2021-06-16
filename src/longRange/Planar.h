
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
	Planar(double cutoffT, double cutoffLJ, Domain* domain, DomainDecompBase* domainDecomposition, ParticleContainer* particleContainer, unsigned slabs, Simulation* simulation);
	virtual ~Planar();

	virtual void init();
	virtual void readXML(XMLfileUnits& xmlconfig);
	virtual void calculateLongRange();
	double lrcLJ(Molecule* mol);
	// For non-equilibrium simulations the density profile must not be smoothed, therefore the density profile from the actual time step is used.
	void directDensityProfile();
	void SetSmoothDensityProfileOption(bool bVal) {_smooth = bVal;}
	virtual void writeProfiles(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep);

private:
	template<typename T>
	void resizeExactly(std::vector<T>& v, unsigned int numElements) const {
		v.reserve(numElements);
		v.resize(numElements);
	}

	void centerCenter(double sig,double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj); 
	void centerSite(double sig,double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj);
	void siteSite(double sig,double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj);
	void dipoleDipole(unsigned ci,unsigned cj,unsigned si,unsigned sj);

	unsigned _slabs;
	unsigned numComp;
	std::vector<unsigned> numLJ;
	std::vector<unsigned> numDipole;
	unsigned numLJSum;
	unsigned numDipoleSum;
	std::vector<unsigned> numLJSum2;
	std::vector<unsigned> numDipoleSum2;
	bool _smooth;
//	bool _dipole;
	std::vector<double> uLJ;
	std::vector<double> vNLJ;
	std::vector<double> vTLJ;
	std::vector<double> vNDLJ;
	std::vector<double> fLJ;
	std::vector<double> rho_g;
	std::vector<double> rho_l;
	std::vector<double> fDipole;
	std::vector<double> uDipole;
	std::vector<double> vNDipole;
	std::vector<double> vTDipole;
	std::vector<double> rhoDipole;
	std::vector<double> rhoDipoleL;
	std::vector<double> muSquare;
	std::vector<double> eLong;
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
