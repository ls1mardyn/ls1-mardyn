/** \file xmlfileUnits.h
  * \brief XML input file with unit handling/conversion support
  * \author Martin Bernreuther <bernreuther@hlrs.de>

  * encapsulates XML main input file with unit attributes
*/

#ifndef XMLFILEUNIT_H
#define XMLFILEUNIT_H

#include "xmlfile.h"

#include <string>
#include <map>
#include <iostream>

//+XMLfileUnits=================================================================
/**
* \class XMLfileUnits
* \brief XML file with unit attributes abstraction
*
* DOM representation of an XML file including unit conversion
*/
class XMLfileUnits : public XMLfile
{
public:
	static const char *const roottag;
	static const char *const refunitstag;
	static const char *const unitattributetag;

	enum UnitType { Unknown_Unit=-1, No_Unit=0,
	                Length_Unit=1,
	                Mass_Unit=2,
	                Time_Unit=3,
	                ElectricCurrent_Unit=4,
	                ThermodynamicTemperature_Unit=5,
	                SubstanceAmount_Unit=6,
	                LuminousIntensity_Unit=7,
	                numUnitTypes
	              };
	// definitions corresponding to enum UnitType in xmlfileUnits.cpp
	static const short int unittypequantifiers[numUnitTypes][8]; // Q
	static const char *const unittypesymbols[numUnitTypes]; // corresponding symbols
	static const char *const unittypequalifiers[numUnitTypes]; // XML qualifiers

	static const int numPrefixes=21;
	static const char *const prefixsymbols[numPrefixes]; // corresponding symbols
	static const short int prefixquantifiers[numPrefixes]; // corresponds to Q[0]

//+XMLfileUnits::ValueUnit======================================================
	class ValueUnit
	{
	public:
		static std::size_t string_endswith(std::string str, const char* suffix); // util/helper function

		/// \brief XMLfileUnits::ValueUnit constructor
		/// \param	double	value
		ValueUnit(double value=0.)
			{ initialize(value); }
		/// \brief XMLfileUnits::ValueUnit constructor
		/// \param	double	value
		/// \param	const std::string&	symbol
		ValueUnit(double value, const std::string& symbol)
			{ initialize(value,symbol); }

		/// \brief XMLfileUnits::ValueUnit copy constructor
		/// duplicate a ValueNode
		/// \param	XMLfileUnits::ValueNode&	source ValueUnit
		ValueUnit(const ValueUnit& vu)
			: m_value(vu.m_value)
			{ for(int i=0;i<numUnitTypes;++i) m_Q[i]=vu.m_Q[i];}

		/// \brief assignment operator
		/// copy/duplicate other ValueUnit content to ValueUnit
		/// \param	const XMLfile::ValueUnit& vu	source ValueUnit
		/// \return	XMLfile::ValueUnit&	reference to this node
		ValueUnit& operator =(const ValueUnit& vu)
			{ m_value=vu.m_value; for(int i=0;i<numUnitTypes;++i) m_Q[i]=vu.m_Q[i]; return *this; }

		/// \brief assignment operator for scalar values
		/// assign (dimensionless) double value
		/// \param	double d	source value
		/// \return	XMLfile::ValueUnit&	reference to this node
		ValueUnit& operator =(const double d)
			{ m_value=d; for(int i=0;i<numUnitTypes;++i) m_Q[i]=0; return *this; }

		/// \brief divide operator
		/// divide other ValueUnit
		/// \return XMLfile::ValueUnit	result
		ValueUnit operator /(const ValueUnit& vu) const;

		double value() const
			{ return m_value; }

		UnitType unittype() const;
		std::string unittypesymbol() const;

		bool iscompatible(const ValueUnit& vu) const
			{ return m_Q[1]==vu.m_Q[1]&&m_Q[2]==vu.m_Q[2]&&m_Q[3]==vu.m_Q[3]&&m_Q[4]==vu.m_Q[4]
			       &&m_Q[5]==vu.m_Q[5]&&m_Q[6]==vu.m_Q[6]&&m_Q[7]==vu.m_Q[7]; }
		bool iscompatible(UnitType ut) const
			{ return ut==unittype(); }

		ValueUnit normalized() const;

		/// \brief print data to stream
		/// print the node data
		/// \param	std::ostream&	stream to write to (default: std::cout)
		void print(std::ostream& ostrm=std::cout, bool simplify=true) const;

	private:

		void initialize(double value=0.);
		void initialize(double value,const std::string& symbol);

		double m_value; // value
		short int m_Q[numUnitTypes]; // quantity exponents: 10^Q[0]*m^Q[1]*g^Q[2]*s^Q[3]*A^Q[4]*K^Q[5]*mol^Q[6]*cd^Q[7]
		                             /* In contrast to the SI base units, gram ("g") instead of kilogram ("kg")
		                                is used as base mass unit here!
		                                This seems to be more consistent related to the handling of prefixes */
	};
//-XMLfileUnits::ValueUnit======================================================

	XMLfileUnits(const std::string& filepath);

	unsigned long getNodeValueUnit(const char* nodepath, ValueUnit& value) const;
	unsigned long getNodeValueUnit(const std::string& nodepath, ValueUnit& value) const
		{ return getNodeValueUnit(nodepath.c_str(),value); }
	unsigned long getNodeValueReduced(const char* nodepath, double& value, UnitType ut=Unknown_Unit) const;
	unsigned long getNodeValueReduced(const std::string& nodepath, double& value, UnitType ut=Unknown_Unit) const
		{ return getNodeValueReduced(nodepath.c_str(),value,ut); }

private:
	ValueUnit m_refunits[numUnitTypes];
};
//-XMLfileUnits=================================================================

/// \brief write a ValueUnit to a stream
/// write XML query data to an output stream
/// \param std::ostream&	output stream
/// \param XMLfile::Query&	query
/// \return std::ostream&	output stream
inline std::ostream& operator << (std::ostream& ostrm, const XMLfileUnits::ValueUnit& vu)
{
	vu.print(ostrm,true);
	return ostrm;
}

#endif
