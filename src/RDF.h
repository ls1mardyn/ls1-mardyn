/*
 * RDF.h
 *
 * @Date: 11.02.2011
 * @Author: eckhardw
 */

#ifndef RDF_H
#define RDF_H

#include <cmath>
#include <string>
#include <vector>

#include "molecules/Molecule.h"

class Domain;
class DomainDecompBase;
class Component;

/**
 * This class calculates the Radial Distribution Function (RDF).
 *
 * The RDF "describes how the atomic density varies as a function of the distance
 * from one particular atom." (see http://en.wikipedia.org/wiki/Radial_distribution_function ).
 * For example, it should be possible to recognize the aggregate state of a system
 * (see http://matdl.org/matdlwiki/index.php/softmatter:Radial_Distribution_Function ).
 *
 * Here the RDF is only determined for atom pairs of which their distance is within
 * the cutoff-radius.
 *
 * Calculation:
 * - calculate the distance of the pair, discretize it to intervalls
 * with length dr (i.e. sort the pairs into bins to obtain a histogram).
 * - for each bin: calculate the number density (i.e. number of particles per volume)
 *   of the corresponding shell
 * - divide the number density by the number density of the system.
 *
 */
class RDF {

	friend class RDFTest;

public:

	/**
	 * @todo Wouldn't make sense to calculate the parameter intervalLength?
	 *       intervalLength = cutoffRadius / bins
	 */
	RDF(double intervalLength, unsigned int bins, const std::vector<Component>& components);

	virtual ~RDF();

	//! @todo put this in the constructor (when the transition to the xml file is done),
	//! or create a seperate output component.
	void setOutputTimestep(unsigned int timestep);

	//! @todo put this in the constructor (when the transition to the xml file is done),
	//! or create a seperate output component.
	void setOutputPrefix(std::string prefix);

	//! plot all the statistics calculated to one or several files
	void doOutput(DomainDecompBase* domainDecomposition, const Domain* domain, unsigned long simStep);

	//! increment the counter indicating for how many iterations
	//! the molecule pairs have been counted.
	void tickRDF() {
		_numberOfRDFTimesteps++;
	}

	//! count the number of molecules per component
	//! @todo: remove it and replace it by component.getNumMolecules()
	void accumulateNumberOfMolecules(std::vector<Component>& components) const;

	void observeRDF(Molecule& mi, Molecule& mj, double dd, double distanceVector[3]) const {
		observeRDF(dd, mi.componentid(), mj.componentid());

		if(_doCollectSiteRDF) {
			double drs[3];
			double dr2;
			unsigned si = mi.numSites();
			unsigned sj = mj.numSites();
			if(si+sj > 2) {
				for(unsigned m = 0; m < si; m++) {
					for(unsigned n = 0; n < sj; n++) {
						const double* dii = mi.site_d(m);
						const double* djj = mj.site_d(n);
						SiteSiteDistance(distanceVector, dii, djj, drs, dr2);
						observeRDF(dr2, mi.componentid(), mj.componentid(), m, n);
					}
				}
			}
		}
	}

	/**
	 * This method "really" counts the number of pairs within a certain distance.
	 */
	void observeRDF(double dd, unsigned i, unsigned j) const {
		if(_numberOfRDFTimesteps <= 0) return;
		if(dd > _maxDistanceSquare) return;
		if(i > j) { this->observeRDF(dd, j, i); return; }
		unsigned l = (unsigned)floor(sqrt(dd)/this->_intervalLength);
		this->_localDistribution[i][j-i][l] ++;
	}

	/**
	 * Count center pairing for particle pair for molecules i and j and centers
	 * m_i, n_j at distance dd.
	 */
	inline void observeRDF(double dd, unsigned i, unsigned j, unsigned m, unsigned n) const {
		if(_numberOfRDFTimesteps <= 0) return;
		if(dd > _maxDistanceSquare) return;
		if(i > j) {
			this->observeRDF(dd, j, i, n, m);
			return;
		}

		unsigned l = (unsigned)floor(sqrt(dd)/this->_intervalLength);
		this->_localSiteDistribution[i][j-i][m][n][l] ++;
                if((i == j) && (m != n)) this->_localSiteDistribution[i][j-i][n][m][l] ++;
//		std::cout << "Observed RDF i=" << i << " j=" << j << " m=" << m << " n=" << n << std::endl;
	}

	bool siteRDF() {
		return this->_doCollectSiteRDF;
	}

	//! set all values counted to 0, except the accumulated ones.
	void reset();

private:


	//! Performs a reduction of the local rdf data of all nodes
	//! to update the "global" fields
	void collectRDF(DomainDecompBase* domainDecomp);

	//! Update the "accumulatedXXX"-fields from the "global"-variables.
	//! @note consequently, collectRDF should be called just before.
	void accumulateRDF();

	void writeToFile(const Domain* domain, const char* prefix, unsigned int i, unsigned int j) const;

	//! The length of an interval
	//! Only used for the output to scale the "radius"-axis.
	const double _intervalLength;

	//! The number of bins, i.e. the number of intervals in which the cutoff
	//! radius will be subdivided.
	unsigned int _bins;

	//! number of different components (i.e. molecule types).
	const unsigned int _numberOfComponents;

	//! components vector
	const std::vector<Component>& _components;

	//! number of timesteps over which the counters are being accumulated
	//! since the last calculation of the RDF.
	int _numberOfRDFTimesteps;

	//! number of timesteps over which the "accumulated"-counters are being
	//! accumulated
	int _accumulatedNumberOfRDFTimesteps;

	//! the maximum distance uthis->_localSiteDistribution[i][j-i][m][n][l]p to which particle pairs are counted, squared
	double _maxDistanceSquare;

	/**
	 * holds the numberOfMolecules of component i at _globalCtr[i], globally.
	 *
	 * @todo remove it, as it can be retrieved via Component::getNumMolecules()
	 */
	unsigned long* _globalCtr;

	//! holds the numberOfMolecules of component i at _globalAccumulatedCtr[i],
	//! globally and accumulated over all timesteps for which particles where
	//! counted for the RDF.
	unsigned long* _globalAccumulatedCtr;

	//! holds the distribution of the neighbouring particles, locally for this process,
	//! i.e. the number of particles of components m and n in bin b: _localDistribution[m][n][b];
	unsigned long ***_localDistribution;

	//! holds the distribution of the neighbouring particles, globally.
	unsigned long ***_globalDistribution;

	//! holds the distribution of the neighbouring particles, globally accumulated.
	unsigned long ***_globalAccumulatedDistribution;

	bool _doCollectSiteRDF;

	unsigned long *****_localSiteDistribution;

	unsigned long *****_globalSiteDistribution;

	unsigned long *****_globalAccumulatedSiteDistribution;

	/**
	 * aggregation interval for the RDF data
	 *
	 * @todo I guess this is redundant with _numberOfRDFTimesteps?
	 */
	unsigned int _RDFOutputTimesteps;

	//! the timestep, the respective component IDs, and "rdf" are
	//! appended to this prefix
	std::string _RDFOutputPrefix;
};

#endif /* RDF_H */
