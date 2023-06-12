#ifndef RDF_H
#define RDF_H

#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "plugins/PluginBase.h"
#include "molecules/Molecule.h"
#include "utils/CommVar.h"

class Component;
class RDFCellProcessor;

/** @brief This class calculates the Radial Distribution Function (RDF).
 *
 * The RDF "describes how the atomic density varies as a function of the distance
 * from one particular atom." (see http://en.wikipedia.org/wiki/Radial_distribution_function ).
 * For example, it should be possible to recognize the aggregate state of a system
 * (see http://matdl.org/matdlwiki/index.php/softmatter:Radial_Distribution_Function ).
 *
 * \note The RDF is only determined for molecule pairs within the cutoff radius of the force and
 *       potential calculation. This means that bins outside the cut-off will always be computed
 *       to be zero.
 *
 * Calculation:
 * - calculate the distance of the pair, discretize it to intervalls
 * with length dr (i.e. sort the pairs into bins to obtain a histogram).
 * - for each bin: calculate the number density (i.e. number of particles per volume)
 *   of the corresponding shell
 * - divide the number density by the number density of the system.

 *Update: Optionally, the RDF can be additionally resolved in angular direction, by chosing angularbins > 1 in the input. Phi is the angle between the central molecules orientation vector and the vector connecting the pair of molecules.
 *The angular coordinate is given and discretized in terms of the cosine of the angle cos(phi), to ensure equal control volume sizes for equal distance R
 */
class RDF : public PluginBase {

	friend class RDFTest;

public:

	RDF();
	virtual ~RDF();

	/** @brief Read in XML configuration for RDFWriter.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <outputplugin name="RDF">
	     <writefrequency>INTEGER</writefrequency>
	     <outputprefix>STRING</outputprefix>
	     <bins>INTEGER</bins>
	     <intervallength>DOUBLE</intervallength>
	   </outputplugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);

	void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep);

	void init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain);

	void finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain);

	std::string getPluginName() {
		return std::string("RDF");
	}
	static PluginBase* createInstance() { return new RDF(); }

	//! @todo put this in the constructor (when the transition to the xml file is done),
	//! or create a seperate output component.
	void setOutputTimestep(unsigned int timestep) { _writeFrequency = timestep; }

	//! @todo put this in the constructor (when the transition to the xml file is done),
	//! or create a seperate output component.
	void setOutputPrefix(std::string prefix) { _outputPrefix = prefix; }

	//! plot all the statistics calculated to one or several files
	void endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomposition, Domain *domain,
                 unsigned long simStep);

	//! increment the counter indicating for how many iterations
	//! the molecule pairs have been counted.
	void tickRDF() {
		_numberOfRDFTimesteps++;
	}

	//! count the number of molecules per component
	//! @todo: remove it and replace it by component.getNumMolecules()
	void accumulateNumberOfMolecules(std::vector<Component>& components);

	void observeARDFMolecule(double dd, double cosPhi, double cosPhiReverse, unsigned cid1, unsigned cid2) {

		if(dd > _maxDistanceSquare) { return; }
		size_t distanceBinID = floor( sqrt(dd) / binwidth() );
		size_t angularBinID = floor( (-cosPhi + 1.)/ angularbinwidth() );
		size_t angularBinIDReverse = floor((-cosPhiReverse + 1.)/ angularbinwidth() );
		size_t binID = distanceBinID * _angularBins + angularBinID;
		size_t binIDReverse = distanceBinID * _angularBins + angularBinIDReverse;
		#if defined _OPENMP
		#pragma omp atomic
		#endif
		_ARDFdistribution.local[cid1][cid2][binID]++;
		#if defined _OPENMP
		#pragma omp atomic
		#endif
		_ARDFdistribution.local[cid2][cid1][binIDReverse]++;
	}

	void observeRDF(Molecule const& mi, Molecule const& mj, double dd) {
		observeRDFMolecule(dd, mi.componentid(), mj.componentid());

		if(isEnabledSiteRDF()) {
			double drs[3];
			double dr2;
			unsigned si = mi.numSites();
			unsigned sj = mj.numSites();
			if(si+sj > 2) {
				for(unsigned m = 0; m < si; m++) {
					// when interacting two molecules with the same component id, the interaction of site m with site n is calculated twice if m!=n.
					// this duplication is corrected by correctionFactor in the writeToFile(...) function.
					// fixing this here, by starting n at m will break the unit tests and some symmetry.
					for(unsigned n = 0; n < sj; n++) {
						const std::array<double,3> dii = mi.site_d_abs(m);
						const std::array<double,3> djj = mj.site_d_abs(n);
						SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
						observeRDFSite(dr2, mi.componentid(), mj.componentid(), m, n);
					}
				}
			}
		}
	}

	/**
	 * This method "really" counts the number of molecule pairs within a certain distance.
	 */
	void observeRDFMolecule(double dd, unsigned i, unsigned j) {
		if(dd > _maxDistanceSquare) { return; }
		if(i > j) { std::swap(j, i); }
		size_t binId = floor( sqrt(dd) / binwidth() );
		#if defined _OPENMP
		#pragma omp atomic
		#endif
		_distribution.local[i][j-i][binId]++;
	}

	/**
	 * Count center pairing for particle pair for molecules i and j and centers
	 * m_i, n_j at distance dd.
	 */
	inline void observeRDFSite(double dd, unsigned i, unsigned j, unsigned m, unsigned n) {
		if(dd > _maxDistanceSquare) { return; }
		if (i > j) {
			std::swap(j, i);
			std::swap(m, n);
		}

		unsigned int binId = floor( sqrt(dd) / binwidth() );
		#if defined _OPENMP
		#pragma omp atomic
		#endif
		_siteDistribution.local[i][j-i][m][n][binId] ++;
		if((i == j) && (m != n)){
			#if defined _OPENMP
			#pragma omp atomic
			#endif
			_siteDistribution.local[i][j-i][n][m][binId] ++;
		}
	}

	bool isEnabledSiteRDF() const { return _doCollectSiteRDF; }
	bool doARDF() const { return _doARDF; }

	void reset();  //!< reset all values to 0, except the accumulated ones.

private:

	template<typename T>
	void resizeExactly(std::vector<T>& v, unsigned int numElements) const {
		v.reserve(numElements);
		v.resize(numElements);
	}

	void init();
	unsigned int numBins() const { return _bins; }
	unsigned int numARDFBins() const { return _ARDFBins; }
	double binwidth() const { return _intervalLength; }
	double angularbinwidth() const { return _angularIntervalLength; }
	void collectRDF(DomainDecompBase* domainDecomp);  //!< update global values from local once

	//! Update the "accumulatedXXX"-fields from the "global"-variables.
	//! @note consequently, collectRDF should be called just before.
	void accumulateRDF();

	void writeToFile(const Domain* domain, const std::string& filename, unsigned int i, unsigned int j) const;
	void writeToFileARDF(const Domain* domain, const std::string& filename, unsigned int i, unsigned int j) const;
	//! The length of an interval
	//! Only used for the output to scale the "radius"-axis.
	double _intervalLength;

	//! The length of an angular interval
	//! Only used for the output to scale the "phi"-axis.
	double _angularIntervalLength;

	//! The number of bins, i.e. the number of intervals in which the cutoff
	//! radius will be subdivided.
	unsigned long _bins;

	//! The number of bins in angular direction in case the angular RDF is
	//! being calculated
	unsigned long _angularBins {1};

	//! The total number of bins for the ARDF is the product of the radial bins and
	//! the angular bins
	unsigned long _ARDFBins;

	//! number of different components (i.e. molecule types).
	unsigned int _numberOfComponents;

	//! components vector
	std::vector<Component>* _components;

	//! sample the RDF every this many iterations
	int _samplingFrequency;

	//! number of timesteps over which the counters are being accumulated
	//! since the last calculation of the RDF.
	int _numberOfRDFTimesteps;

	//! number of timesteps over which the "accumulated"-counters are being
	//! accumulated
	int _accumulatedNumberOfRDFTimesteps;

	//! the maximum distance up to which particle pairs are counted, squared
	double _maxDistanceSquare;

	/**
	 * holds the numberOfMolecules of component i at _globalCtr[i], globally.
	 * accumulates over time steps as the number of molecules may change
	 */
	std::vector<unsigned long> _globalCtr;

	//! holds the numberOfMolecules of component i at _globalAccumulatedCtr[i],
	//! globally and accumulated over all timesteps for which particles where
	//! counted for the RDF.
	std::vector<unsigned long> _globalAccumulatedCtr;

	//! holds the distribution of the neighbouring particles, locally for this process,
	//! i.e. the number of particles of components m and n in bin b: _distribution.local[m][n][b] or _distribution.global[m][n][b];
	CommVar<std::vector<std::vector<std::vector<unsigned long>>>> _distribution;
	CommVar<std::vector<std::vector<std::vector<unsigned long>>>> _ARDFdistribution;
	//! holds the distribution of the neighbouring particles, globally accumulated.
	std::vector<std::vector<std::vector<unsigned long>>> _globalAccumulatedDistribution;
	std::vector<std::vector<std::vector<unsigned long>>> _globalAccumulatedARDFDistribution;

	bool _doCollectSiteRDF;
	bool _doARDF;
	// vector indices:
	// first: component i
	// second: component j, but shifted by -i, so that we only have to save it once for each component interaction (i<->j is the same as j<->i)
	// third: site m of first component
	// fourth: site n of second component
	// fifth: bin to sort the particles into to get the actual RDF.
	CommVar<std::vector<std::vector<std::vector<std::vector<std::vector<unsigned long>>>>>> _siteDistribution;

	std::vector<std::vector<std::vector<std::vector<std::vector<unsigned long>>>>> _globalAccumulatedSiteDistribution;

	unsigned int _writeFrequency;  //!< aggregation and output writing interval for the RDF data
	std::string _outputPrefix;  //!< output prefix for rdf files

	bool _initialized;
	bool _readConfig;

	RDFCellProcessor * _cellProcessor;
};

#endif /* RDF_H */
