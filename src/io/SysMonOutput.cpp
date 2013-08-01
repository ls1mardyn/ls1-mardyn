#include "io/SysMonOutput.h"

#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"

using Log::global_log;

SysMonOutput::SysMonOutput() : _writeFrequency(1) {}


void SysMonOutput::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << std::endl;
}

void SysMonOutput::initOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain){
	SysMon* sysmon = SysMon::getSysMon();
	sysmon->updateExpressionValues();
	//sysmon->writeExpressionValues();
	global_log->info() << *sysmon;
}

void SysMonOutput::doOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep,
			std::list<ChemicalPotential>* lmu) {
	if((simstep % _writeFrequency) == 0) {
		SysMon* sysmon = SysMon::getSysMon();
		sysmon->updateExpressionValues();
		//sysmon->writeExpressionValues();
		global_log->info() << *sysmon;
	}
}

void SysMonOutput::finishOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain){
	SysMon* sysmon = SysMon::getSysMon();
	sysmon->updateExpressionValues();
	//sysmon->writeExpressionValues();
	global_log->info() << *sysmon;
}
