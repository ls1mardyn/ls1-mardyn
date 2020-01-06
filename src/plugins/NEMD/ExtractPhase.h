/*
 * ExtractPhase.h
 *
 *  Created on: 31.12.2019
 *      Author: mheinen
 */

#ifndef EXTRACT_PHASE_H_
#define EXTRACT_PHASE_H_

#include "plugins/PluginBase.h"
#include <cstdint>

class ParticleContainer;
class DomainDecompBase;
class Domain;
class ChemicalPotential;
class CavityEnsemble;
class ExtractPhase : public PluginBase
{
public:
	// constructor and destructor
	ExtractPhase();
	~ExtractPhase();

	/** @brief Read in XML configuration for ExtractPhase and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
		<plugin name="ExtractPhase">
			<change>
				<enabled>1</enabled>
				<cid_ub>2</cid_ub>
			</change>
			<density>
				<value>0.66318</value>
				<percent>0.5</percent>
				<cutoff>2.5</cutoff>
				<range>
					<left>250.0</left>
					<right>350.0</right>
				</range>
			</density>
			<interface>
				<type>vl</type>
				<range>
					<left>190.0</left>
					<right>210.0</right>
				</range>
			</interface>
		</plugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	void init(ParticleContainer *particleContainer,
			  DomainDecompBase *domainDecomp, Domain *domain) override;

	void beforeForces(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) override;

    /** @brief Method afterForces will be called after forcefields have been applied
     *
     * make pure Virtual ?
     */
    void afterForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    ) override;

	void endStep(
			ParticleContainer *particleContainer,
			DomainDecompBase *domainDecomp, Domain *domain,
			unsigned long simstep) override {}

	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain) override {}

	std::string getPluginName() override {return std::string("ExtractPhase");}
	static PluginBase* createInstance() {return new ExtractPhase();}

private:
	struct ActionsDone {
		bool beforeForces, afterForces;
	} _bDone;

	struct ComponentChange {
		bool enabled;
		uint32_t cid_ub;
	} _compChange;

	struct DensityTarget {
		double value;
		double percent;
		double cutoff;
		struct Range {
			double left, right;
			bool enabled;
		} range;
	} _densityTarget;

	enum InterfaceType : uint16_t {
		INTT_VAPOR_LIQUID = 1,
		INTT_LIQUID_VAPOR = 2
	};

	struct InterfaceParams {
		InterfaceType type;
		struct Range {
			double left, right;
		} range;
	} _interface;
};

#endif /* EXTRACT_PHASE_H_ */
