#ifndef DRIFTCTRL_H_
#define DRIFTCTRL_H_

#include "plugins/PluginBase.h"
#include "utils/CommVar.h"

#include <string>
#include <map>
#include <list>
#include <cstdint>
#include <array>
#include <vector>

#define FORMAT_SCI_MAX_DIGITS std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10)

struct BinVectors {
	CommVar<std::vector<uint64_t> > numParticles;
	std::array<CommVar<std::vector<double> >,3> velocity;
	std::array<std::vector<double>,3> velo_corr;
};

class ParticleContainer;
class DomainDecompBase;
class Domain;
class ChemicalPotential;
class CavityEnsemble;
class DriftCtrl : public PluginBase
{
public:
	// constructor and destructor
	DriftCtrl();
	~DriftCtrl();

	/** @brief Read in XML configuration for DriftCtrl and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
		<plugin name="DriftCtrl">
			<control>
				<freq>
					<sample>INT</sample>		<!-- sample frequency -->
					<control>INT</control>		<!-- control frequency -->
					<write>INT</write>			<!-- write frequency -->
				</freq>
			</control>
			<target>
				<cid>INT</cid>														<!-- target component id -->
				<drift> <vx>DOUBLE</vx> <vy>DOUBLE</vy> <vz>DOUBLE</vz> </drift>	<!-- target drift vector (vx,vy,vz) -->
			</target>
			<range>
				<yl>DOUBLE</yl> <yr>DOUBLE</yr>		<!-- range in which drift will be controled -->
				<subdivision>
					<binwidth>DOUBLE</binwidth>		<!-- binwidth of bins the range is subdivided -->
				</subdivision>
			</range>
		</plugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	void init(ParticleContainer *particleContainer,
			  DomainDecompBase *domainDecomp, Domain *domain) override;

	void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;
	void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override {}
	void endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain, unsigned long simstep) override {}
	void finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) override {}
	std::string getPluginName() override {return std::string("DriftCtrl");}
	static PluginBase* createInstance() {return new DriftCtrl();}

private:
	struct Control{
		struct Freq{
			uint32_t sample, control, write;
		} freq;
	} _control;
	struct Range{
		double yl, yr, width;
		struct Subdivision {
			uint32_t numBins;
			struct BinWidth {
				double init, actual;
			} binWidth;
		} subdivision;
	} _range;
	
	struct Target {
		std::array<double,3> drift;
		uint32_t cid;
	} _target;
	
	std::vector<BinVectors> _sampling;
};

#endif /*DRIFTCTRL_H_*/
