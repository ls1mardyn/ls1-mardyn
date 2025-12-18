/*
 * HardwareInfo.h
 *
 *  Created on: 16 Dec 2025
 *      Author: amartyads
 */

#pragma once

#include <vector>

#include "PluginBase.h"

struct ThreadwiseInfo {
	int thread, totalThreads;
	unsigned int cpuID, numa;
	ThreadwiseInfo(int threadNum = 0, int totalThreadsNum = 1, unsigned int openMPCPUID = 0,
				   unsigned int openMPNUMA = 0)
		: thread(threadNum), totalThreads(totalThreadsNum), cpuID(openMPCPUID), numa(openMPNUMA) {}
};

/**
 * @brief Prints the hardware info of current run to stdout or to file
 * @author Amartya Das Sharma
 *
 */
class HardwareInfo : public PluginBase {
public:
	HardwareInfo() = default;
	~HardwareInfo() override = default;

	void readXML(XMLfileUnits& xmlconfig) override;
	void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;
	void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
				 unsigned long simstep) override {};
	void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override {};

	std::string getPluginName() override { return "HardwareInfo"; }

	static PluginBase* createInstance() { return new HardwareInfo(); }

private:
	void populateData(DomainDecompBase* domainDecomp);
	void printDataToStdout();
	void writeDataToFile();
	std::string convertFullDataToJson();

	std::string _filename;
	int _rank, _totalRanks;
	std::string _processorName;
	std::vector<ThreadwiseInfo> threadData;
	bool _dataPopulated = false;
};