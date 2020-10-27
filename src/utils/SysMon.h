/** \file SysMon.h
  * \brief System monitoring class
  * \author Martin Bernreuther <bernreuther@hlrs.de>
*/

#ifndef SYSMON_H
#define SYSMON_H

/*#define ENABLE_MPI*/

#ifdef _POSIX_VERSION
#define SYSMON_ENABLE_SYSCONF
#endif
#ifdef __linux__
#define SYSMON_ENABLE_SYSINFO
#endif
#define SYSMON_ENABLE_MALLINFO
#define SYSMON_ENABLE_PROCMEMINFO
#define SYSMON_ENABLE_PROCVMSTAT
#define SYSMON_ENABLE_PROCLOADAVG
#define SYSMON_ENABLE_PROCSELFSTATM
#define SYSMON_ENABLE_PROCSELFSCHED
//#define SYSMON_ENABLE_PROCSELFSCHEDSTAT
#define SYSMON_ENABLE_PROCSELFSTATUS

#include <string>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <iostream>
#include <utility>

#include <cmath>


#ifdef ENABLE_MPI
#include "mpi.h"
#endif

#ifdef SYSMON_ENABLE_SYSCONF
#include <unistd.h>	// sysconf, sync; see http://pubs.opengroup.org/onlinepubs/9699919799/basedefs/unistd.h.html
#endif

#ifdef SYSMON_ENABLE_SYSINFO
#include <sys/sysinfo.h>
#endif

#ifdef SYSMON_ENABLE_MALLINFO
#include <malloc.h>	// mallinfo
#endif


#include "Expression.h"

// SysMon ------------------------------------------------------------------------------------------
class SysMon{
public:
	typedef double Tvalue;
#ifdef MPI_VERSION
	static const MPI_Datatype mpiTvalue;
	static const int mpiRootRank;
#endif
	
private:	// Singleton
	static SysMon* s_sysmoninstance;
	SysMon(const SysMon&);
#ifdef MPI_VERSION
	SysMon(MPI_Comm mpicomm=MPI_COMM_WORLD) : _mpicomm(mpicomm) { _variableset=new Expression::VariableSet(); }
#else
	SysMon() { _variableset=new Expression::VariableSet(); }
#endif
	~SysMon() { clear(); delete(_variableset); s_sysmoninstance=NULL; }
	
public:
	static SysMon* getSysMon()
	{
		if(s_sysmoninstance==NULL)
			s_sysmoninstance = new SysMon();
		return s_sysmoninstance;
	}
	
	void clear();
	int addExpression(const std::string& exprstr);
	int addExpression(const std::string& exprstr, const std::string& label)
		{ int idx=addExpression(exprstr); if(idx>=0) _expressions.back().setLabel(label); return idx; }
	unsigned int numExpressions() { return _expressions.size(); }
	void updateExpressionValues(bool resetMinMax=false);
	int getExpressionIndex(const std::string& label) const;
	Tvalue getExpressionValue(unsigned int index) const
		{ if(index<_values.size()) return _values[index]; else return Tvalue(); }
	std::pair<Tvalue,Tvalue> getExpressionMinMaxPeakValues(unsigned int index) const
		{ if(index*2+1<_valuesMaxMinPeak.size()) return std::make_pair(_valuesMaxMinPeak[index*2],_valuesMaxMinPeak[index*2+1]); else return std::make_pair(0,0); }
#ifdef MPI_VERSION
	std::pair<Tvalue,Tvalue> getExpressionMinMaxValues(unsigned int index) const;
#endif
	void writeExpressionValues(std::ostream& ostrm=std::cout
	                         , std::string header=std::string()
	                         , std::string lineprefix=std::string()
	                         , std::string sep=std::string("\t")
	                         , std::string eol=std::string("\n")
	                          ) const;
	std::string InfoString(std::string header=std::string()
	                     , std::string lineprefix=std::string()
	                     , std::string sep=std::string("\t")
	                     , std::string eol=std::string("\n")
	                      ) const
	{
		std::ostringstream oss;
		writeExpressionValues(oss,header,lineprefix,sep,eol);
		return oss.str();
	}
	operator std::string() const { return InfoString(); }
	
private:
	Expression::VariableSet* _variableset;
	std::list<Expression> _expressions;
	std::vector<Tvalue> _values;
	std::vector<bool> _initMinMax;
	// *valuesMaxMin* stores for each value (i) the
	//  maximum (2*i) and negative minimum (=maximum of negative value) (2*i+1)
	std::vector<Tvalue> _valuesMaxMinPeak;
#ifdef MPI_VERSION
	std::vector<Tvalue> _valuesMaxMin;
	MPI_Comm _mpicomm;
#endif
	
private:
	unsigned int updateVariables_sysconf();
	unsigned int updateVariables_sysinfo();
	unsigned int updateVariables_mallinfo();
	unsigned int updateVariables_procmeminfo();
	unsigned int updateVariables_procvmstat();
	unsigned int updateVariables_procloadavg();
	unsigned int updateVariables_procselfstatm();
	unsigned int updateVariables_procselfschedstat();
	unsigned int updateVariables_procselfsched();
	unsigned int updateVariables_procselfstatus();
};
// ------------------------------------------------------------------------------------------ SysMon


inline std::ostream& operator << (std::ostream& ostrm, const SysMon& s)
{
	s.writeExpressionValues(ostrm);
	return ostrm;
}

#endif
