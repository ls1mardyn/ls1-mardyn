/*
 * Dropaligner.h
 *
 *  Created on: 22 June 2020
 *      Author: Marx
 */

#pragma once

// class DropalignerTest;
#include "Domain.h"
#include "PluginBase.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"

/** @brief
 * Plugin: can be enabled via config.xml <br>
 *
 * Calculates Center of mass of a spherical droplet and moves all particles to align with the desired center of mass<br>
 * Alignment happens once every interval-simsteps<br>
 * The correction factor can be set from 0-1<br>
 * 1 being full alignment -> 0 no alignment at all<br>
 * <b>HALO must not be present</b> for the alignment. Halo would lead to incorrect alignment.<br>
 * This is guarenteed by calling the alignment in the beforeForces step of the simulation
 * \code{.xml}
 * <plugin name="Dropaligner">
 *			<xpos>50</xpos>
 *			<ypos>60</ypos>
 *			<zpos>50</zpos>
 *           <radius>30</radius>
 *			<interval>1</interval>
 *			<correctionFactor>.5</correctionFactor>
 * </plugin>
 * \endcode
 */

class Dropaligner : public PluginBase {
private:
	bool _enabled = true;

	int _interval = 1;

	double _alignmentCorrection = 1.0;
	double _motion[3];
	double _balance[3];
	double _mass = 0.0;
	double _boxLength[3];
	double _radius;
	double _xPos;
	double _yPos;
	double _zPos;

public:
	void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override {
		Log::global_log->debug() << "DropletRealignment enabled" << std::endl;

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

	std::string getPluginName() override { return std::string("Dropaligner"); }

	static PluginBase* createInstance() { return new Dropaligner(); }
};
