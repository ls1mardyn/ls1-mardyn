/*
 * HardwareInfo.h
 *
 *  Created on: 16 Dec 2025
 *      Author: amartyads
 */

#pragma once

#include "PluginBase.h"

/**
 * @brief Prints the hardware info of current run to stdout or to file
 * @author Amartya Das Sharma
 *
 */
class HardwareInfo : public PluginBase {
public:
	HardwareInfo() = default;
	~HardwareInfo() override = default;

	void readXML(XMLfileUnits& xmlconfig) override {}
	void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
	void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
				 unsigned long simstep) override {};
	void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override {};

	std::string getPluginName() override { return "HardwareInfo"; }

	static PluginBase* createInstance() { return new HardwareInfo(); }

private:
	int _thread, _rank, _totalThreads, _totalRanks;
	unsigned int _openMPCPUID, _openMPNUMA;
	char _processorName[1024];
};