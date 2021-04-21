#pragma once

#include "PluginBase.h"

/**
 * Plugin that produces a load imbalace by sleeping.
 *
 * Current modes:
 * - "varyingSteps" -- introduces need for overlapping globals.
 *     - all even ranks: sleep at afterForces
 *     - all uneven ranks: sleep at endStep
 */
class LoadImbalanceThroughSleepPlugin : public PluginBase {
public:
	/**
	 * \code{.xml}
	 * <plugin name="LoadImbalanceThroughSleepPlugin">
	 *   <varyingSteps>True/False</varyingSteps> <!-- Default: True -->
	 *   <varyingStepsSleepTime>INT: sleep time in ms</varyingStepsSleepTime> <!-- Default: 200 -->
	 * </plugin>
	 * \endcode
	 * @param xmlconfig
	 */
	void readXML(XMLfileUnits &xmlconfig) override;

	void init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) override{};

	void afterForces(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
					 unsigned long simstep) override;

	void endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain,
				 unsigned long simstep) override;

	void finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) override{};

	std::string getPluginName() override{return "LoadImbalanceThroughSleepPlugin";}

	static PluginBase* createInstance() { return new LoadImbalanceThroughSleepPlugin(); }

private:
	bool _varyingSteps{true};
	int _varyingStepsSleepTime{200};
};
