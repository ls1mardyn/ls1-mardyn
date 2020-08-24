/*
 * FixRegion.h
 *
 *  Created on: 9 July 2020
 *      Author: Marx
 */

#ifndef MARDYN_TRUNK_FIXREGION_H
#define MARDYN_TRUNK_FIXREGION_H

// class DropalignerTest;
#include "Domain.h"
#include "PluginBase.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"

//#include "parallel/DomainDecompBase.h"

/** @brief
* Plugin: can be enabled via config.xml <br>
*
* Fixes particles in a given region.<br>

* \code{.xml}
* <plugin name="FixRegion">
*			<xmin>50</xmin>
*			<ymin>60</ymin>
*			<zmin>50</zmin>
*			<xmax>50</xmax>
*			<ymax>60</ymax>
*			<zmax>50</zmax>
* </plugin>
* \endcode
*/

class FixRegion : public PluginBase {
private:
	double _xMin;
	double _yMin;
	double _zMin;
	double _xMax;
	double _yMax;
	double _zMax;
	double _boxLength[3];
	unsigned long _molCount;

public:
	/*  Dropaligner(){
		// SETUP
		   for (unsigned d = 0; d < 3; d++) {
			   _balance[d] = 0.0;
			   _motion[d] = 0.0;
		   }
	   };
	   ~Dropaligner(){};*/

	void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
	//        global_log -> debug() << "DropletRealignment enabled" << std::endl;

	//      for(unsigned d = 0; d < 3; d++){
	//           _boxLength[d] = domain->getGlobalLength(d);
	//        }
	//    }

	void readXML(XMLfileUnits& xmlconfig) override;

	void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
					  unsigned long simstep) override;

	void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
					 unsigned long simstep) override;

	void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
				 unsigned long simstep) override;

	void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override{};

	std::string getPluginName() override { return std::string("FixRegion"); }

	static PluginBase* createInstance() { return new FixRegion(); }
	
};

#endif  // MARDYN_TRUNK_FIXREGION_H
