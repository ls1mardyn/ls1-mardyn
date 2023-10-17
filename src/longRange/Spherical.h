
//Calculation of the surface tension in a system with spherical interfaces needs a Long Range Correction.
//
//The correction terms are based on Janecek (2006) and Lustig (1988).

#ifndef SPHERICAL_H_
#define SPHERICAL_H_

#include "LongRangeCorrection.h"

#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <cstdint>

#include "molecules/MoleculeForwardDeclaration.h"
class Domain;
class ParticleContainer;

class Spherical:public LongRangeCorrection {
public:
	Spherical(double cutoffT,double cutoffLJ,Domain* domain,  DomainDecompBase* domainDecomposition, ParticleContainer* particleContainer, Simulation* simulation);
	virtual ~Spherical();

	virtual void init();
	virtual void readXML(XMLfileUnits& xmlconfig);
	virtual void calculateLongRange();
	virtual void writeProfiles(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep);

private:
	template<typename T>
	void resizeExactly(std::vector<T>& v, unsigned int numElements) const {
		v.reserve(numElements);
		v.resize(numElements);
	}

	// void centerCenter(double sig,double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj); 
	// void centerSite(double sig,double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj);
	// void siteSite(double sig,double eps,unsigned ci,unsigned cj,unsigned si, unsigned sj);

	void calcDensityProfiles();

	double RhoP(double r, double rhov, double rhol, double D0, double R0);
	double SICCu(int n, double r);
	double SICSu(int n, double r, double tau);
	double CS(int n, double r, double tau);
	double SISSu(int n, double r, double tau1, double tau2);
	double SS(int n, double r, double tau1, double tau2);
	double SSLN(double r, double tau1, double tau2);
	double TICCu(int n, double rc, double sigma2);
	double TICSu(int n, double rc, double sigma2, double tau);
	double TISSu(int n, double rc, double sigma2, double tau1, double tau2);
	double TICCp(int n, double rc, double sigma2);
	double TICSp(int n, double rc, double sigma2, double tau);
	double TISSp(int n, double rc, double sigma2, double tau1, double tau2);

	unsigned int numComp;
	unsigned long globalNumMols;
	double rc;
	double rcmax;
	double UpotKorrLJ;
	double VirialKorrLJ;
	double TempRho;
	bool droplet;
	std::string _outputPrefix;
	double _T;

	unsigned int NShells;
	unsigned int NSMean;
	double _drShells;
	double _deltaShells;
	std::vector<double> RShells;
	std::vector<double> RShells2;
	std::vector<double> RShells3;
	std::vector<double> VShells;
	std::vector<double> rhoShellsTemp;
	std::vector<double> rhoShellsTemp_global;
	std::vector<double> rhoShellsMean;
	std::vector<double> rhoShellsAvg;
	std::vector<double> rhoShellsAvg_global;
	std::vector<double> TShellsTemp;
	std::vector<double> TShellsAvg;
	std::vector<double> TShellsAvg_global;
	std::vector<double> rhoShells;
	std::vector<double> rhoShellsT;
	std::vector<double> rhoShells_global;
	std::vector<double> PartShells;
	std::vector<double> UShells_Mean;
	std::vector<double> UShells_Mean_global;
	std::vector<double> FShells_Mean;
	std::vector<double> FShells_Mean_global;
	std::vector<double> PNShells_Mean;
	std::vector<double> PNShells_Mean_global;
	std::vector<double> PTShells_Mean;
	std::vector<double> PTShells_Mean_global;
	std::vector<double> VirShells_Mean;
	std::vector<double> VirShells_Mean_global;
	std::vector<double> VirShells_Corr;
	std::vector<double> VirShells_Corr_global;
	std::vector<double> VirShells_N;
	std::vector<double> VirShells_N_global;
	std::vector<double> VirShells_T;
	std::vector<double> VirShells_T_global;
	std::vector<unsigned> numLJ;
	unsigned int numLJSum;
	std::vector<unsigned> numLJSum2;
	std::vector<double> ksi;
	std::vector<double> FcorrX;
	std::vector<double> FcorrY;
	std::vector<double> FcorrZ;
	std::vector<double> FcorrX_global;
	std::vector<double> FcorrY_global;
	std::vector<double> FcorrZ_global;
	std::vector<double> lowerS;
	std::vector<double> interS;
	std::vector<double> upperS;
	std::vector<double> eLong;
	double boxlength[3];
	double systemcenter[3];

	std::string filenameTanhParams;
	std::string filenameGlobalCorrs;
	std::string filenameThermData;
	
	ParticleContainer* _particleContainer;
	Domain* _domain;
	DomainDecompBase* _domainDecomposition;
	
};


#endif /*SPHERICAL_H_*/
