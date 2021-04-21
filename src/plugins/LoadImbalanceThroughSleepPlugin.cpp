#include "LoadImbalanceThroughSleepPlugin.h"

#include "parallel/DomainDecompBase.h"
#include "utils/xmlfileUnits.h"

#include <chrono>
#include <thread>


void LoadImbalanceThroughSleepPlugin::readXML(XMLfileUnits &xmlconfig) {
	xmlconfig.getNodeValue("varyingSteps", _varyingSteps);
}

void LoadImbalanceThroughSleepPlugin::afterForces(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
												  unsigned long simstep) {
	using namespace std::chrono_literals;
	if (_varyingSteps and domainDecomp->getRank() % 2 == 0) {
		std::this_thread::sleep_for(2000ms);
	}
}

void LoadImbalanceThroughSleepPlugin::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
											  Domain *domain, unsigned long simstep) {
	using namespace std::chrono_literals;
	if (_varyingSteps and domainDecomp->getRank() % 2 != 0) {
		std::this_thread::sleep_for(2000ms);
	}
}