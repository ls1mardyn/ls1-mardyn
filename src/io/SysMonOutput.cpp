#include "io/SysMonOutput.h"

#include <sstream>

#include "utils/Logger.h"
#include "utils/SysMon.h"
#include "utils/xmlfileUnits.h"

using Log::global_log;

SysMonOutput::SysMonOutput() : _writeFrequency(1) {}


void SysMonOutput::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << std::endl;
	SysMon* sysmon = SysMon::getSysMon();
	XMLfile::Query query = xmlconfig.query("expression");
	//string oldpath = xmlconfig.getcurrentnodepath();
	for(XMLfile::Query::const_iterator exprIter = query.begin(); exprIter; exprIter++ )
	{
		xmlconfig.changecurrentnode(exprIter);
		std::string expr(xmlconfig.getNodeValue_string("."));
		std::string label(xmlconfig.getNodeValue_string(std::string("@label")));
		if(label.empty())
		{
			sysmon->addExpression(expr);
		} else {
			sysmon->addExpression(expr,label);
		}
	}
	//xmlconfig.changecurrentnode(oldpath);
	
	//sysmon->updateExpressionValues();
	//global_log->info() << sysmon->InfoString("System Monitor\n","\t");
}

void SysMonOutput::init(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
                        Domain * /*domain*/){
	SysMon* sysmon = SysMon::getSysMon();
	sysmon->updateExpressionValues();
	global_log->info() << sysmon->InfoString("System Monitor (initial)\n","\t");
}

void SysMonOutput::endStep(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
                           Domain * /*domain*/, unsigned long simstep)
{
	if((simstep % _writeFrequency) == 0) {
		SysMon* sysmon = SysMon::getSysMon();
		sysmon->updateExpressionValues();
		std::ostringstream oss;
		oss << "System Monitor (simulation step " << simstep << ")" << std::endl;
		global_log->info() << sysmon->InfoString(oss.str(),"\t");
	}
}

void SysMonOutput::finish(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
						  Domain * /*domain*/){
	SysMon* sysmon = SysMon::getSysMon();
	sysmon->updateExpressionValues();
	global_log->info() << sysmon->InfoString("System Monitor (final)\n","\t");
}
