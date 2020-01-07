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
				<enabled>INT</enabled>   <!-- 0:disable | 1:enable component change instead of particle deletion -->
				<cid_ub>INT</cid_ub>     <!-- change component id of overlapping particles to cid_ub=INT, first component: cid_ub=1 (unity based)
			</change>
			<density>
				<value>FLOAT</value>        <!-- use this value as liquid density -->
				<percent>FLOAT</percent>    <!-- percentage of liquid density identifies begin of vapor phase -->
				<cutoff>FLOAT</cutoff>      <!-- cutoff radius for calculating local density -->
				<range>					    <!-- range in which liquid density is calculated -->
					<left>DOUBLE</left>     <!-- left boundary -->
					<right>DOUBLE</right>   <!-- right boundary -->
				</range>
			</density>
			<interface>
				<type>vl</type>            <!-- 1,vl,VL:vapor-liquid or 2,lv,LV:liquid-vapor interface -->
				<range>			           <!-- searching for overlaps in specified range (y axis) -->
					<left>FLOAT</left>     <!-- left boundary -->
					<right>FLOAT</right>   <!-- right boundary -->
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
