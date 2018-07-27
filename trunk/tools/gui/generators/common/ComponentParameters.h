/*
 * ComponentParameters.h
 *
 * @Date: 30.05.2011
 * @Author: eckhardw
 */

#ifndef COMPONENTPARAMETERS_H_
#define COMPONENTPARAMETERS_H_

#include "Parameters/ParameterCollection.h"
#include "Parameters/ParameterWithDoubleValue.h"
#include "molecules/Component.h"

/**
 * This class represents a collection of parameters corresponding to a
 * component type in Mardyn/ls1.
 */
class ComponentParameters: public ParameterCollection {

public:

	/**
	 * Construct a parameter collection for a component.
	 */
	ComponentParameters(const std::string& id, const std::string& name,
			const std::string& desc, Component& component);

	virtual ~ComponentParameters();

	/**
	 * set the value of the parameter p in the component
	 */
	static void setParameterValue(Component& component, Parameter* p,
			std::string& valueName);

	virtual void load(const std::string& filename, Generator* generator);

private:

	/**
	 * set the value of the parameter p in the charge
	 */
	static void setParameterValue(Charge& charge, ParameterWithDoubleValue* p,
			std::string& valueName);

	/**
	 * set the value of the parameter p in the lennard jones center
	 */
	static void setParameterValue(LJcenter& ljCenter,
			ParameterWithDoubleValue* p, std::string& valueName);

	/**
	 * set the value of the parameter p in the dipole
	 */
	static void setParameterValue(Dipole& dipole, ParameterWithDoubleValue* p,
			std::string& valueName);

	/**
	 * set the value of the parameter p in the quadrupoles
	 */
	static void setParameterValue(Quadrupole& quadrupole,
			ParameterWithDoubleValue* p, std::string& valueName);
};

#endif /* COMPONENTPARAMETERS_H_ */
