/** \file xmlfileUnits.cpp
  * \brief XML input file with unit handling/conversion support
  * \author Martin Bernreuther <bernreuther@hlrs.de>

  * encapsulates XML main input file with unit attributes
*/

#include "xmlfileUnits.h"

#include <iostream>

#include <cstring>
#include <cmath>

const char *const XMLfileUnits::roottag = { "/mardyn" };
const char *const XMLfileUnits::refunitstag = { "refunits" };
const char *const XMLfileUnits::unitattributetag = { "unit" };

const short int XMLfileUnits::unittypequantifiers[numUnitTypes][8] = { { 0,0,0,0,0,0,0,0 }
                                                                      ,{ 0,1,0,0,0,0,0,0 }
                                                                      ,{ 0,0,1,0,0,0,0,0 }
                                                                      ,{ 0,0,0,1,0,0,0,0 }
                                                                      ,{ 0,0,0,0,1,0,0,0 }
                                                                      ,{ 0,0,0,0,0,1,0,0 }
                                                                      ,{ 0,0,0,0,0,0,1,0 }
                                                                      ,{ 0,0,0,0,0,0,0,1 }
                                                                     };
const char *const XMLfileUnits::unittypesymbols[numUnitTypes] = { ""
                                                                 ,"m"
                                                                 ,"g"
                                                                 ,"s"
                                                                 ,"A"
                                                                 ,"K"
                                                                 ,"mol"
                                                                 ,"cd"
                                                                };
const char *const XMLfileUnits::unittypequalifiers[numUnitTypes] = { "no_unit"
                                                                    ,"length"
                                                                    ,"mass"
                                                                    ,"time"
                                                                    ,"current"
                                                                    ,"temperature"
                                                                    ,"substance"
                                                                    ,"luminous"
                                                                   };

const char *const XMLfileUnits::prefixsymbols[numPrefixes] = { "y"  // "yocto"
                                                              ,"z"  // "zepto"
                                                              ,"a"  // "atto"
                                                              ,"f"  // "femto"
                                                              ,"p"  // "pico"
                                                              ,"n"  // "nano"
                                                              ,"u"  // "micro","µ"
                                                              ,"m"  // "milli"
                                                              ,"c"  // "centi"
                                                              ,"d"  // "deci"
                                                              ,""   // ""
                                                              ,"da" // "deca"
                                                              ,"h"  // "hecto"
                                                              ,"k"  // "kilo"
                                                              ,"M"  // "mega"
                                                              ,"G"  // "giga"
                                                              ,"T"  // "tera"
                                                              ,"P"  // "peta"
                                                              ,"E"  // "exa"
                                                              ,"Z"  // "zetta"
                                                              ,"Y"  // "yotta"
                                                             };
// prefixes are sorted in ascending order
const short XMLfileUnits::prefixquantifiers[numPrefixes] = { -24
                                                            ,-21
                                                            ,-18
                                                            ,-15
                                                            ,-12
                                                            ,-9
                                                            ,-6
                                                            ,-3
                                                            ,-2
                                                            ,-1
                                                            ,0
                                                            ,1
                                                            ,2
                                                            ,3
                                                            ,6
                                                            ,9
                                                            ,12
                                                            ,15
                                                            ,18
                                                            ,21
                                                            ,24
                                                           };

size_t XMLfileUnits::ValueUnit::string_endswith(std::string str, const char* suffix)
{
	unsigned int suffixlength=strlen(suffix);
	if (suffixlength>str.length()) return std::string::npos;
	size_t pos = str.length()-suffixlength;
	//return str.rfind(suffix,pos);
	if (str.substr(pos,suffixlength)==std::string(suffix))
		return pos;
	else
		return std::string::npos;
}

XMLfileUnits::ValueUnit XMLfileUnits::ValueUnit::operator /(const XMLfileUnits::ValueUnit& vu) const
{
	ValueUnit res;
	if(vu.m_value==0.) return res;
	res.m_value=m_value/vu.m_value;
	for(int i=0;i<numUnitTypes;++i) res.m_Q[i]=m_Q[i]-vu.m_Q[i];
	return res;
}

XMLfileUnits::UnitType XMLfileUnits::ValueUnit::unittype() const
{
	for(int i=0;i<numUnitTypes;++i)
	{
		if(  m_Q[1]==unittypequantifiers[i][1]
		  && m_Q[2]==unittypequantifiers[i][2]
		  && m_Q[3]==unittypequantifiers[i][3]
		  && m_Q[4]==unittypequantifiers[i][4]
		  && m_Q[5]==unittypequantifiers[i][5]
		  && m_Q[6]==unittypequantifiers[i][6]
		  && m_Q[7]==unittypequantifiers[i][7])
			return static_cast<UnitType>(i);
	}
	return Unknown_Unit;
}

std::string XMLfileUnits::ValueUnit::unittypesymbol() const
{
	for(int i=0;i<numUnitTypes;++i)
	{
		if(  m_Q[1]==unittypequantifiers[i][1]
		  && m_Q[2]==unittypequantifiers[i][2]
		  && m_Q[3]==unittypequantifiers[i][3]
		  && m_Q[4]==unittypequantifiers[i][4]
		  && m_Q[5]==unittypequantifiers[i][5]
		  && m_Q[6]==unittypequantifiers[i][6]
		  && m_Q[7]==unittypequantifiers[i][7])
			return unittypesymbols[i];
	}
	return std::string();
}

XMLfileUnits::ValueUnit XMLfileUnits::ValueUnit::normalized() const
{
	short int value_exponent = short(log10(m_value));
	short int new_q0 = ((m_Q[0]+value_exponent)/3)*3;
	//int i;
	//for(i=0;i<numPrefixes;++i)
	//	if (prefixquantifiers[i]>m_Q[0]+value_exponent) break;
	//if(i>0) --i;
	//new_q0=prefixquantifiers[i];
	short int q0diff = m_Q[0]-new_q0;
	double new_v = m_value*pow(10., q0diff);
	XMLfileUnits::ValueUnit new_vu(*this);
	new_vu.m_value = new_v;
	new_vu.m_Q[0] = new_q0;
	return new_vu;
}

void XMLfileUnits::ValueUnit::print(std::ostream& ostrm,bool simplify) const
{
	XMLfileUnits::UnitType ut=unittype();
	double v=m_value;
	short int q0=m_Q[0];
	if(simplify)
	{
		ValueUnit vu=normalized();
		v=vu.value();
		q0=vu.m_Q[0];
		if(ut!=Unknown_Unit)
		{
			// determine the prefix
			int i;
			for(i=0;i<numPrefixes;++i)
			{
				if (prefixquantifiers[i]>q0) break;
			}
			if(i>0) --i;
			q0=m_Q[0]-prefixquantifiers[i];
			v=m_value*pow(10., q0);
			ostrm << v << " " << std::string(prefixsymbols[i]) << std::string(unittypesymbols[ut]);
			return;
		}
	}
	// SI units representation
	// g -> kg
	if(m_Q[2]>0) q0-=3*m_Q[2];
	else if(m_Q[2]<0) q0+=3*m_Q[2];
	if(q0) v*=pow(10., q0);
	ostrm << v;
	if(m_Q[1])
	{
		if(m_Q[1]==1) ostrm << " m";
		else ostrm << " m^"<<int(m_Q[1]);
	}
	if(m_Q[2])
	{
		if(m_Q[2]==1) ostrm << " kg";
		else ostrm << " kg^"<<int(m_Q[2]);
	}
	if(m_Q[3])
	{
		if(m_Q[3]==1) ostrm << " s";
		else ostrm << " s^"<<int(m_Q[3]);
	}
	if(m_Q[4])
	{
		if(m_Q[4]==1) ostrm << " A";
		else ostrm << " A^"<<int(m_Q[4]);
	}
	if(m_Q[5])
	{
		if(m_Q[5]==1) ostrm << " K";
		else ostrm << " K^"<<int(m_Q[5]);
	}
	if(m_Q[6])
	{
		if(m_Q[6]==1) ostrm << " mol";
		else ostrm << " mol^"<<int(m_Q[6]);
	}
	if(m_Q[7])
	{
		if(m_Q[7]==1) ostrm << " cd";
		else ostrm << " cd^"<<int(m_Q[7]);
	}
}

void XMLfileUnits::ValueUnit::initialize(double value)
	{
		m_value=value;
		m_Q[0]=m_Q[1]=m_Q[2]=m_Q[3]=m_Q[4]=m_Q[5]=m_Q[6]=m_Q[7]=0;
	}

void XMLfileUnits::ValueUnit::initialize(double value,const std::string& symbol)
{
	initialize(value);
	size_t pos;
	if ((pos=string_endswith(symbol,"mol"))!=std::string::npos)
	{// mole, SubstanceAmount_Unit
		//m_Q[1]=m_Q[2]=m_Q[3]=m_Q[4]=m_Q[5]=m_Q[7]=0;
		m_Q[6]=1;
	}
	else if ((pos=string_endswith(symbol,"cd"))!=std::string::npos)
	{// candela, LuminousIntensity_Unit
		//m_Q[1]=m_Q[2]=m_Q[3]=m_Q[4]=m_Q[5]=m_Q[6]=0;
		m_Q[7]=1;
	}
	else if ((pos=string_endswith(symbol,"m"))!=std::string::npos)
	{// metre, Length_Unit
		//m_Q[2]=m_Q[3]=m_Q[4]=m_Q[5]=m_Q[6]=m_Q[7]=0;
		m_Q[1]=1;
	}
	else if ((pos=string_endswith(symbol,"g"))!=std::string::npos)
	{// gram (used instead of SI base unit "kg" as base unit), Mass_Unit
		//m_Q[1]=m_Q[3]=m_Q[4]=m_Q[5]=m_Q[6]=m_Q[7]=0;
		m_Q[2]=1;
	}
	else if ((pos=string_endswith(symbol,"s"))!=std::string::npos)
	{// second, Time_Unit
		//m_Q[1]=m_Q[2]=m_Q[4]=m_Q[5]=m_Q[6]=m_Q[7]=0;
		m_Q[3]=1;
	}
	else if ((pos=string_endswith(symbol,"A"))!=std::string::npos)
	{// ampere, ElectricCurrent_Unit
		//m_Q[1]=m_Q[2]=m_Q[3]=m_Q[5]=m_Q[6]=m_Q[7]=0;
		m_Q[4]=1;
	}
	else if ((pos=string_endswith(symbol,"K"))!=std::string::npos)
	{// kelvin, ThermodynamicTemperature_Unit
		//m_Q[1]=m_Q[2]=m_Q[3]=m_Q[4]=m_Q[6]=m_Q[7]=0;
		m_Q[5]=1;
	}
	// PREFIX ------------------------------------------------------
	if(pos!=std::string::npos && pos>0)
	{ // could determinate the unit, but prefix is present
		std::string prefix=symbol.substr(0,pos);
		if(prefix==std::string(""))
		{ // no prefix (will never be executed)
			m_Q[0]=0;
		}
		else if(prefix==std::string("da"))
		{ // deca-
			m_Q[0]=1;
		}
		else if(prefix==std::string("h"))
		{ // hecto-
			m_Q[0]=2;
		}
		else if(prefix==std::string("k"))
		{ // kilo-
			m_Q[0]=3;
		}
		else if(prefix==std::string("M"))
		{ // mega-
			m_Q[0]=6;
		}
		else if(prefix==std::string("G"))
		{ // giga-
			m_Q[0]=9;
		}
		else if(prefix==std::string("T"))
		{ // tera-
			m_Q[0]=12;
		}
		else if(prefix==std::string("P"))
		{ // peta-
			m_Q[0]=15;
		}
		else if(prefix==std::string("E"))
		{ // exa-
			m_Q[0]=18;
		}
		else if(prefix==std::string("Z"))
		{ // zetta-
			m_Q[0]=21;
		}
		else if(prefix==std::string("Y"))
		{ // yotta-
			m_Q[0]=24;
		}
		else if(prefix==std::string("d"))
		{ // deci-
			m_Q[0]=-1;
		}
		else if(prefix==std::string("c"))
		{ // centi-
			m_Q[0]=-2;
		}
		else if(prefix==std::string("m"))
		{ // milli-
			m_Q[0]=-3;
		}
		else if(prefix==std::string("u") )
		       // || prefix==std::string("µ"))
		{ // micro-
			m_Q[0]=-6;
		}
		else if(prefix==std::string("n"))
		{ // nano-
			m_Q[0]=-9;
		}
		else if(prefix==std::string("p"))
		{ // pico-
			m_Q[0]=-12;
		}
		else if(prefix==std::string("f"))
		{ // femto-
			m_Q[0]=-15;
		}
		else if(prefix==std::string("a"))
		{ // atto-
			m_Q[0]=-18;
		}
		else if(prefix==std::string("z"))
		{ // zepto-
			m_Q[0]=-21;
		}
		else if(prefix==std::string("y"))
		{ // yocto-
			m_Q[0]=-24;
		}
	}
	/*
	else
		// no prefix
		m_Q[0]=0;
	*/
}

//------------------------------------------------------------------------------

XMLfileUnits::XMLfileUnits(const std::string& filepath)
	: XMLfile(filepath)
{
	if (changecurrentnode(std::string(roottag)+std::string(refunitstag))) {
		double v;
		std::string u;
		for(unsigned int i=1;i<numUnitTypes;++i)
		{
			if(getNodeValue(unittypequalifiers[i],v))
			{
				getNodeValue(std::string(unittypequalifiers[i])+std::string("@")+std::string(unitattributetag),u);
				m_refunits[i]=ValueUnit(v,u);
			}
			else
				m_refunits[i]=ValueUnit(1.,"");
		}
		changecurrentnode("/");
	}
}

unsigned long XMLfileUnits::getNodeValueUnit(const char* nodepath, XMLfileUnits::ValueUnit& value) const
{
	double v;
	unsigned long found=getNodeValue(nodepath,v);
	std::string unitattrpath=std::string(nodepath)+std::string("@")+std::string(unitattributetag);
	std::string u;
	if(getNodeValue(unitattrpath,u))
	{
		value=XMLfileUnits::ValueUnit(v,u);
	}
	else {
		value = v;
	}
	return found;
}

unsigned long XMLfileUnits::getNodeValueReduced(const char* nodepath, double& value, XMLfileUnits::UnitType ut) const
{
	XMLfileUnits::ValueUnit vu;
	unsigned int found=getNodeValueUnit(nodepath,vu);
	if(found)
	{
		if(vu.unittype()==No_Unit)
		{ // read value is already dimensionless
			value=vu.value();
		}
		else
		{
			if(ut==Unknown_Unit) ut=vu.unittype();	// no unittype given - grab from read unit
			if(ut!=Unknown_Unit && ut==vu.unittype())
			{
				XMLfileUnits::ValueUnit vu_red=vu/m_refunits[ut];
				//cout << "reduced value " << nodepath << " from " << vu << " to " << vu_red << std::endl;
				value=vu_red.value();
			}
			else
			{
				//cerr << "ERROR reading " << nodepath << ", found " << found <<" entries (" << vu << ")" << std::endl;
				found=0;
			}
		}
	}
	return found;
}
