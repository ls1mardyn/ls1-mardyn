/*
 * FixRegion.h
 *
 *  Created on: 9 July 2020
 *      Author: Marx
 */

#pragma once

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
	void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;

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

	void registerCallbacks(std::map<std::string, FunctionWrapper>& callbackMap) override {
		callbackMap["FixRegion::getMoleculesInRegion"] = [this] { return _molCount; };
	}
};
