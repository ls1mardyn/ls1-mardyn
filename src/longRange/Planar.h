
//Calculation of the surface tension in a system with planar interfaces needs a Long Range Correction.
//
//The correction terms are based on Janecek (2006) and Lustig (1988).

#ifndef PLANAR_H_
#define PLANAR_H_

#include "LongRangeCorrection.h"

#include "utils/ObserverBase.h"
#include "utils/Region.h"

#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <cstdint>

#include "molecules/MoleculeForwardDeclaration.h"
class Domain;
class ParticleContainer;

class Planar : public LongRangeCorrection, public ObserverBase, public ControlInstance {
public:
	Planar(double cutoffT, double cutoffLJ, Domain* domain, DomainDecompBase* domainDecomposition, ParticleContainer* particleContainer, unsigned slabs, Simulation* simulation);
	~Planar() override = default;

	/** @brief Read in XML configuration for Planar and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
		<longrange type="planar">
			<region> <!-- y coordinates of left and right boundaries within correction of force is applied to particles; pot. energy and virial correction is always applied since it is necessary for the correct calculation of the state values -->
				<left refcoordsID="INT">FLOAT</left> <!-- Reference of coordinate can be set (see DistControl); 0: origin (default) | 1:left interface | 2:right interface -->
				<right refcoordsID="INT">FLOAT</right>
			</region>
			<slabs>INT</slabs> <!-- Domain is divided into INT slabs -->
			<smooth>0</smooth>
			<frequency>10</frequency> <!-- Frequency at which LRC is recalculated -->
			<writecontrol> <!-- Parameters to control output in file -->
				<start>900000</start>
				<frequency>100000</frequency>
				<stop>5000000</stop>
			</writecontrol>
		</longrange>
	   \endcode
	 */

	void init() override;
	void readXML(XMLfileUnits& xmlconfig) override;
	void calculateLongRange() override;
	// For non-equilibrium simulations the density profile must not be smoothed, therefore the density profile from the actual time step is used.
	void directDensityProfile();
	void SetSmoothDensityProfileOption(bool bVal) {_smooth = bVal;}
	void writeProfiles(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep) override;
	
	// Get potential energy correction per molecule
	double getUpotCorr(Molecule* mol) override;

	// Observer, ControlInstance
	SubjectBase* getSubject();
	void update(SubjectBase* subject) override;
	std::string getShortName() override {return "Planar";}

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
	struct RegionPos {
		int refPosID[2];  // kind of reference position, see DistControl
		double refPos[2]; // left and right boundary (y coord) set in config.xml
		double actPos[2]; // left and right boundary (y coord) within the correction is applied
	} _region;
	SubjectBase* _subject;
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
