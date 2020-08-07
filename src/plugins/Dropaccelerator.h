/*
 * Dropaccelerator.h
 *
 *  Created on: 26 June 2020
 *      Author: Koch
 */

#ifndef MARDYN_TRUNK_DROPACCELERATOR_H
#define MARDYN_TRUNK_DROPACCELERATOR_H

// class DropacceleratorTest;
#include "Domain.h"
#include "PluginBase.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"

/** @brief
 * Plugin: can be enabled via config.xml <br>
 *
 * Detects the undirected particle movement of every molecule in a spherical droplet and accelerates the droplet to a
 *velocity in -y direction <br> Acceleratiom happens once every interval-simsteps in determid time steps until the
 *velocity has reached<br> \code{.xml} <plugin name="Dropaccelerator"> <xposition>45.6701</xposition>
 *	<yposition>100</yposition>
 *	<zposition>45.6701</zposition>
 *	<dropradius>10</dropradius>
 *	<velocity>0.5</velocity>
 *   <starttime>100</starttime>
 *	<steps>100</steps>
 * </plugin>
 * \endcode
 */

class Dropaccelerator : public PluginBase {
private:
	// friend DropacceleratorTest;

	bool _enabled = true;

	int _interval = 1;

	double _motion[3];
	double _balance[3];
	double _mass = 0.0;
	double _boxLength[3];
	double _dropRadius;
	double _xPosition;
	double _yPosition;
	double _zPosition;
	double _veloc;
	double _steps;
	double _velocNow;
	double _startSimStep;

	std::vector<char> _particleIsInDroplet;

public:
	/*  Dropaccelerator(){
		// SETUP
		   for (unsigned d = 0; d < 3; d++) {
			   _balance[d] = 0.0;
			   _motion[d] = 0.0;
		   }
	   };
	   ~Dropaccelerator(){};*/

	void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override {
		global_log->debug() << "DropletAccelerator enabled" << std::endl;

		for (unsigned d = 0; d < 3; d++) {
			_boxLength[d] = domain->getGlobalLength(d);
		}
	}

	void readXML(XMLfileUnits& xmlconfig) override;

	void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
					  unsigned long simstep) override;

	void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
				 unsigned long simstep) override;

	void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override{};

	std::string getPluginName() override { return std::string("Dropaccelerator"); }

	static PluginBase* createInstance() { return new Dropaccelerator(); }
};

#endif  // MARDYN_TRUNK_DROPACCELERATOR_H
