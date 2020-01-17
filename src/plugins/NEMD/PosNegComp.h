/*
 * PosNegComp.h
 *
 *  Created on: 03.12.2018
 *      Author: mheinen
 */

#ifndef POSNEGCOMP_H_
#define POSNEGCOMP_H_

#include "plugins/PluginBase.h"
#include <cstdint>

class ParticleContainer;
class DomainDecompBase;
class Domain;
class ChemicalPotential;
class CavityEnsemble;
class PosNegComp : public PluginBase
{
public:
	// constructor and destructor
	PosNegComp();
	~PosNegComp();

	/** @brief Read in XML configuration for PosNegComp and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
		<plugin name="PosNegComp">                                                  <!-- component change according to moving direction of particles -->
			<cid_ub> <pos>INT</pos> <neg>INT</neg> <ignore>INT</ignore> </cid_ub>   <!-- component id of particles moving in pos:positive | neg:negative y-direction, ignore: particles are not affected -->
			<limit_y> <left>FLOAT</left> <right>FLOAT</right> </limit_y>            <!-- component change active in specified range -->
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

	std::string getPluginName() override {return std::string("PosNegComp");}
	static PluginBase* createInstance() {return new PosNegComp();}

private:
	void changeComponents(ParticleContainer* particleContainer);

private:
	struct CompIDS {
		uint32_t pos;
		uint32_t neg;
		uint32_t ignore;
	} _cid_ub;  // ub: unity based
	struct LimitY {
		double left, right;
	} _limitY;
};



#endif /* POSNEGCOMP_H_ */
