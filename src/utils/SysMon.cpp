/** \file SysMon.cpp
  * \brief System monitoring class
  * \author Martin Bernreuther <bernreuther@hlrs.de>
*/

#define SYSMON_CPP

#include "SysMon.h"
#include "Logger.h"	// MPI_CHECK, Log::global_log

#include <stack>

#include <fstream>
#include <sstream>

using Log::global_log;

SysMon* SysMon::s_sysmoninstance=NULL;

#ifdef MPI_VERSION
	const MPI_Datatype SysMon::mpiTvalue=MPI_DOUBLE;
	const int SysMon::mpiRootRank=0;
#endif


void SysMon::clear()
{
	for (std::list<Expression>::iterator exprit=_expressions.begin();exprit!=_expressions.end();++exprit)
		exprit->clear();
	_expressions.clear();
	_values.clear();
	_initMinMax.clear();
	_valuesMaxMinPeak.clear();
#ifdef MPI_VERSION
	_valuesMaxMin.clear();
#endif
}

int SysMon::addExpression(const std::string& exprstr)
{
	_expressions.push_back(Expression(std::string(),_variableset));
	Expression& expr=_expressions.back();
	expr.initializeRPN(exprstr,true);	// generate a Label
	//expr.genLabel();	// done at the end of initializeRPN
	
	_values.push_back(0.);
	//_values.resize(_expressions.size(),0.);
	_initMinMax.push_back(true);
	//_initMinMax.resize(_expressions.size(),true);
	_valuesMaxMinPeak.push_back(0.); _valuesMaxMinPeak.push_back(0.);
	//_valuesMaxMinPeak.resize(_values.size()*2,0.);
#ifdef MPI_VERSION
	_valuesMaxMin.push_back(0.); _valuesMaxMin.push_back(0.);
	//_valuesMaxMin.resize(_valuesMaxMinPeak.size(),0.);
#endif
	return _expressions.size()-1;
}


void SysMon::updateExpressionValues(bool resetMinMax)
{
	if(_expressions.empty()) return;	// no expressions?
	
	if(resetMinMax) _initMinMax.assign(numExpressions(),true);
	
	//sync();
	if(_variableset->existVariableGroup("sysconf")) updateVariables_sysconf();
	if(_variableset->existVariableGroup("sysinfo")) updateVariables_sysinfo();
	if(_variableset->existVariableGroup("mallinfo")) updateVariables_mallinfo();
	if(_variableset->existVariableGroup("procmeminfo")) updateVariables_procmeminfo();
	if(_variableset->existVariableGroup("procvmstat")) updateVariables_procvmstat();
	if(_variableset->existVariableGroup("procloadavg")) updateVariables_procloadavg();
	if(_variableset->existVariableGroup("procselfstatm")) updateVariables_procselfstatm();
	if(_variableset->existVariableGroup("procselfschedstat")) updateVariables_procselfschedstat();
	if(_variableset->existVariableGroup("procselfsched")) updateVariables_procselfsched();
	if(_variableset->existVariableGroup("procselfstatus")) updateVariables_procselfstatus();
	
	size_t i=0;
	for(std::list<Expression>::const_iterator exprit=_expressions.begin();exprit!=_expressions.end();++exprit)
	{
		_values[i]=exprit->evaluateFloat();
		++i;
	}
	
	std::vector<Tvalue> valuesMaxMin(2*_values.size());
	for(i=0;i<_values.size();++i)
	{
		valuesMaxMin[2*i]=_values[i];
		valuesMaxMin[2*i+1]=-_values[i];
	}
	
	for(i=0;i<_valuesMaxMinPeak.size();++i)
	{
		if(_initMinMax[i/2] || valuesMaxMin[i]>_valuesMaxMinPeak[i])
			_valuesMaxMinPeak[i]=valuesMaxMin[i];
	}
	
#ifdef MPI_VERSION
	int myrank;
	MPI_CHECK( MPI_Comm_rank(_mpicomm,&myrank) );
	MPI_CHECK( MPI_Reduce(&valuesMaxMin[0],&_valuesMaxMin[0],_valuesMaxMin.size(),mpiTvalue,MPI_MAX,0,MPI_COMM_WORLD) );
	
	if(myrank==mpiRootRank)
	{
		for(i=0;i<_valuesMaxMinPeak.size();++i)
		{
			if(_initMinMax[i/2] || _valuesMaxMin[i]>_valuesMaxMinPeak[i])
				_valuesMaxMinPeak[i]=_valuesMaxMin[i];
		}
	}
#endif
	_initMinMax.assign(_initMinMax.size(),false);
}

int SysMon::getExpressionIndex(const std::string& label) const
{
	int idx=0;
	for(std::list<Expression>::const_iterator it=_expressions.begin();it!=_expressions.end();++it)
	{
		if(it->getLabel()==label) return idx;
		++idx;
	}
	return -1;
}

#ifdef MPI_VERSION
std::pair<SysMon::Tvalue,SysMon::Tvalue> SysMon::getExpressionMinMaxValues(unsigned int index) const
{
	int myrank;
	MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD,&myrank) );
	if(myrank==mpiRootRank&&index*2+1<_valuesMaxMin.size())
		return std::make_pair(_valuesMaxMin[index*2],_valuesMaxMin[index*2+1]);
	else
		return std::make_pair(0,0);
}
#endif


void SysMon::writeExpressionValues(std::ostream& ostrm, std::string header, std::string lineprefix, std::string sep, std::string eol) const
{
	size_t numvalues=_values.size();
	size_t i=0;
#ifdef MPI_VERSION
	// MPI_Comm_get_attr(_mpicomm,MPI_HOST,&mpihost,&flag);
	int myrank;
	MPI_CHECK( MPI_Comm_rank(_mpicomm,&myrank) );
	if(myrank==mpiRootRank) {
		ostrm << header;
		for(std::list<Expression>::const_iterator exprit=_expressions.begin();exprit!=_expressions.end();++exprit)
		{
			ostrm << lineprefix << exprit->getLabel();
			if(i>=numvalues || _initMinMax[i])
			{
				ostrm << sep << "undefined";
			} else {
				ostrm << sep << "[" << -_valuesMaxMin[2*i+1] << "," << _valuesMaxMin[2*i] << "]";
				ostrm << sep << "[" << -_valuesMaxMinPeak[2*i+1] << "," << _valuesMaxMinPeak[2*i] << "]";
				++i;
			}
			ostrm << eol;
		}
	}
#else
	ostrm << header;
	for(std::list<Expression>::const_iterator exprit=_expressions.begin();exprit!=_expressions.end();++exprit)
	{
		ostrm << lineprefix << exprit->getLabel();
		if(i>=numvalues || _initMinMax[i])
		{
			ostrm << sep << "undefined";
		} else {
			ostrm << sep << _values[i];
			ostrm << sep << "[" << -_valuesMaxMinPeak[2*i+1] << "," << _valuesMaxMinPeak[2*i] << "]";
			++i;
		}
		ostrm << eol;
	}
#endif
}


unsigned int SysMon::updateVariables_sysconf()
{
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_SYSCONF
	long val=0;
	
#ifdef __linux__
	val=sysconf(_SC_PHYS_PAGES);
	_variableset->setVariable("sysconf:PHYS_PAGES",val);
	++numvalues;
	val=sysconf(_SC_AVPHYS_PAGES);
	_variableset->setVariable("sysconf:AVPHYS_PAGES",val);
	++numvalues;
#endif
	val=sysconf(_SC_PAGESIZE);
	_variableset->setVariable("sysconf:PAGESIZE",val);
	++numvalues;
#endif
	return numvalues;
}

unsigned int SysMon::updateVariables_sysinfo()
{
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_SYSINFO
	long val=0;
	
	struct sysinfo sysinfodata;
	int rc=sysinfo(&sysinfodata);
	
	if(!rc)
	{
		val=long(sysinfodata.uptime);
		_variableset->setVariable("sysinfo","uptime",val);
		++numvalues;
		val=long(sysinfodata.loads[0]);
		_variableset->setVariable("sysinfo","loads1",val);
		++numvalues;
		val=long(sysinfodata.loads[1]);
		_variableset->setVariable("sysinfo","loads5",val);
		++numvalues;
		val=long(sysinfodata.loads[2]);
		_variableset->setVariable("sysinfo","loads15",val);
		++numvalues;
		val=long(sysinfodata.totalram);
		_variableset->setVariable("sysinfo","totalram",val);
		++numvalues;
		val=long(sysinfodata.freeram);
		_variableset->setVariable("sysinfo","freeram",val);
		++numvalues;
		val=long(sysinfodata.sharedram);
		_variableset->setVariable("sysinfo","sharedram",val);
		++numvalues;
		val=long(sysinfodata.bufferram);
		_variableset->setVariable("sysinfo","bufferram",val);
		++numvalues;
		val=long(sysinfodata.totalswap);
		_variableset->setVariable("sysinfo","totalswap",val);
		++numvalues;
		val=long(sysinfodata.freeswap);
		_variableset->setVariable("sysinfo","freeswap",val);
		++numvalues;
		val=long(sysinfodata.procs);
		_variableset->setVariable("sysinfo","procs",val);
		++numvalues;
		val=long(sysinfodata.totalhigh);
		_variableset->setVariable("sysinfo","totalhigh",val);
		++numvalues;
		val=long(sysinfodata.freehigh);
		_variableset->setVariable("sysinfo","freehigh",val);
		++numvalues;
		val=long(sysinfodata.mem_unit);
		_variableset->setVariable("sysinfo","mem_unit",val);
		++numvalues;
	}
#endif
	return numvalues;
}

unsigned int SysMon::updateVariables_mallinfo()
{
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_MALLINFO
	long val=0;
	
	struct mallinfo mallinfodata=mallinfo();
	
	val=long(mallinfodata.arena);
	_variableset->setVariable("mallinfo","arena",val);
	++numvalues;
	val=long(mallinfodata.ordblks);
	_variableset->setVariable("mallinfo","ordblks",val);
	++numvalues;
	val=long(mallinfodata.smblks);
	_variableset->setVariable("mallinfo","smblks",val);
	++numvalues;
	val=long(mallinfodata.hblks);
	_variableset->setVariable("mallinfo","hblks",val);
	++numvalues;
	val=long(mallinfodata.hblkhd);
	_variableset->setVariable("mallinfo","hblkhd",val);
	++numvalues;
	val=long(mallinfodata.usmblks);
	_variableset->setVariable("mallinfo","usmblks",val);
	++numvalues;
	val=long(mallinfodata.fsmblks);
	_variableset->setVariable("mallinfo","fsmblks",val);
	++numvalues;
	val=long(mallinfodata.uordblks);
	_variableset->setVariable("mallinfo","uordblks",val);
	++numvalues;
	val=long(mallinfodata.fordblks);
	_variableset->setVariable("mallinfo","fordblks",val);
	++numvalues;
	val=long(mallinfodata.keepcost);
	_variableset->setVariable("mallinfo","keepcost",val);
	++numvalues;
#endif
	return numvalues;
}

unsigned int SysMon::updateVariables_procmeminfo()
{	// 
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_PROCMEMINFO
	std::ifstream ifstrm("/proc/meminfo");
	if(!ifstrm) return 0;
	long val=0;
	std::string line,label,valunit;
	unsigned int linenr=0;
	while (getline(ifstrm, line))
	{	// process line
		++linenr;
		std::istringstream iss(line);
		if (!(iss >> label >> val)) { break; }
		size_t i=label.find_first_of(" :");
		while(i!=std::string::npos)
		{
			label.erase(i,1);
			i=label.find_first_of(" :");
		}
		valunit=std::string();	//=""
		if(!iss.eof()) if (!(iss >> valunit)) break;
		//iss.str(std::string());iss.clear();
		if(!valunit.empty())
		{
			if(valunit=="kB")
				val*=1024;
			else
				//cerr << "WARNING /proc/meminfo:" << linenr << ": unknown unit " << valunit << " (no conversion): using " << label << "=" << val << std::endl;
				global_log->warning() << "WARNING /proc/meminfo:" << linenr << ": unknown unit " << valunit << " (no conversion): using " << label << "=" << val << std::endl;
		}
		_variableset->setVariable("procmeminfo",label,val);
		++numvalues;
	}
	ifstrm.close();
#endif
	return numvalues;
}

unsigned int SysMon::updateVariables_procvmstat()
{	//
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_PROCVMSTAT
	std::ifstream ifstrm("/proc/vmstat");
	if(!ifstrm) return 0;
	long val=0;
	std::string line,label;
	unsigned int linenr=0;
	while (getline(ifstrm, line))
	{	// process line
		++linenr;
		std::istringstream iss(line);
		if (!(iss >> label >> val)) { break; }
		//iss.str(std::string());iss.clear();
		_variableset->setVariable("procvmstat",label,val);
		++numvalues;
	}
	ifstrm.close();
#endif
	return numvalues;
}

unsigned int SysMon::updateVariables_procloadavg()
{	// see also `man 5 proc | grep -m 1 -A 12 /proc/loadavg`
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_PROCLOADAVG
	std::ifstream ifstrm("/proc/loadavg");
	if(!ifstrm) return 0;
	std::string line,label;
	if (!getline(ifstrm, line)) return 0;
	ifstrm.close();
	size_t i=line.find("/");
	if(i!=std::string::npos)
	{
		//line[i]=' ';
		line.replace(i,1," ");
	}
	float loadavg1,loadavg5,loadavg15;
	unsigned int numschedentexec,numschedentexist,pidrecent;
	std::istringstream iss(line);
	if (!(iss >> loadavg1 >> loadavg5 >> loadavg15
	          >> numschedentexec >> numschedentexist >> pidrecent)) { return 0; }
	//iss.str(std::string());iss.clear();
	_variableset->setVariable("procloadavg","loadavg1",double(loadavg1));
	++numvalues;
	_variableset->setVariable("procloadavg","loadavg5",double(loadavg5));
	++numvalues;
	_variableset->setVariable("procloadavg","loadavg15",double(loadavg15));
	++numvalues;
	_variableset->setVariable("procloadavg","numschedentexec",long(numschedentexec));
	++numvalues;
	_variableset->setVariable("procloadavg","numschedentexist",long(numschedentexist));
	++numvalues;
	//_variableset->setVariable("procloadavg","pidrecent",long(pidrecent)); ++numvalues;
#endif
	return numvalues;
}

unsigned int SysMon::updateVariables_procselfstatm()
{	// see also `man 5 proc | grep -m 1 -A 12 /statm`
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_PROCSELFSTATM
	std::ifstream ifstrm("/proc/self/statm");
	if(!ifstrm) return 0;
	unsigned long size,resident,share,text,lib,data;
	ifstrm >> size >> resident >> share >> text >> lib >> data;
	ifstrm.close();
	_variableset->setVariable("procselfstatm","size",long(size));
	++numvalues;
	_variableset->setVariable("procselfstatm","resident",long(resident));
	++numvalues;
	_variableset->setVariable("procselfstatm","share",long(share));
	++numvalues;
	_variableset->setVariable("procselfstatm","text",long(text));
	++numvalues;
	_variableset->setVariable("procselfstatm","lib",long(lib));
	++numvalues;
	_variableset->setVariable("procselfstatm","data",long(data));
	++numvalues;
#endif
	return numvalues;
}

unsigned int SysMon::updateVariables_procselfsched()
{
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_PROCSELFSCHED
	std::ifstream ifstrm("/proc/self/sched");
	if(!ifstrm) return 0;
	double val=0;
	std::string line,label,colon;
	unsigned int linenr=0;
	while (getline(ifstrm, line))
	{	// process line
		++linenr;
		std::istringstream iss(line);
		if (!(iss >> label >> colon >> val)) { continue; }
		//iss.str(std::string());iss.clear();
		if(colon!=":") continue;
		_variableset->setVariable("procselfsched",label,val);
		++numvalues;
	}
	ifstrm.close();
#endif
	return numvalues;
}

unsigned int SysMon::updateVariables_procselfschedstat()
{
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_PROCSELFSCHEDSTAT
	std::ifstream ifstrm("/proc/self/schedstat");
	if(!ifstrm) return 0;
	unsigned long runningtime, waitingtime;
	unsigned int numtasks;
	if (!(ifstrm >> runningtime >> waitingtime >> numtasks)) { return 0; }
	ifstrm.close();
	_variableset->setVariable("procselfschedstat:runningtime",long(runningtime));
	++numvalues;
	_variableset->setVariable("procselfschedstat:waitingtime",long(waitingtime));
	++numvalues;
	_variableset->setVariable("procselfschedstat:numtasks",long(numtasks));
	++numvalues;
#endif
	return numvalues;
}

unsigned int SysMon::updateVariables_procselfstatus()
{	// 
	unsigned int numvalues=0;
#ifdef SYSMON_ENABLE_PROCSELFSTATUS
	std::ifstream ifstrm("/proc/self/status");
	if(!ifstrm) return 0;
	long val=0;
	std::string line,label,valunit;
	unsigned int linenr=0;
	while (getline(ifstrm, line))
	{	// process line
		++linenr;
		std::istringstream iss(line);
		if (!(iss >> label)) { continue; }
		//iss.str(std::string());iss.clear();
		size_t i=label.find_first_of(" :");
		while(i!=std::string::npos)
		{
			label.erase(i,1);
			i=label.find_first_of(" :");
		}
		if(label.find("Vm")==0 || label.find("Rss")==0 || label.find("Hugetlb")==0)
		{
			if (!(iss >> val)) break;
			valunit=std::string();	//=""
			if(!iss.eof()) if (!(iss >> valunit)) break;
			if(!valunit.empty())
			{
				if(valunit=="kB")
					val*=1024;
				else
					//cerr << "WARNING /proc/self/status:" << linenr << ": unknown unit " << valunit << " (no conversion): using " << label << "=" << val << std::endl;
					global_log->warning() << "WARNING /proc/self/status:" << linenr << ": unknown unit " << valunit << " (no conversion): using " << label << "=" << val << std::endl;
			}
			_variableset->setVariable("procselfstatus",label,val);
		}
		++numvalues;
	}
	ifstrm.close();
#endif
	return numvalues;
}

