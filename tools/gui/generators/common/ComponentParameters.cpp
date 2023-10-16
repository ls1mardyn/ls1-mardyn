/*
 * ComponentParameters.cpp
 *
 * @Date: 30.05.2011
 * @Author: eckhardw
 */

#include "ComponentParameters.h"

#include "Parameters/ParameterCollection.h"
#include "Parameters/ParameterWithIntValue.h"
#include "Parameters/ParameterWithDoubleValue.h"

#include "PMFileReader.h"
#include "Conversions.h"
#include "../MDGenerator.h"

#include <cstdlib>

//#define DEBUG


ComponentParameters::ComponentParameters(const std::string& id,
		const std::string& name, const std::string& desc, Component& component) :
	ParameterCollection(id, name, desc, Parameter::BUTTON, true, false) {
//	parameters.push_back(new ParameterWithDoubleValue(name + ".Temperature", "Temperature", "Temperature of this component", Parameter::LINE_EDIT,
//			false,component.T()));
	parameters.push_back(
			new ParameterWithIntValue(name + ".NumberOfLJCenters",
					"Number of LJCenters", "number of LJCenters",Parameter::SPINBOX,
					true, component.numLJcenters()));
	parameters.push_back(
			new ParameterWithIntValue(name + ".NumberOfCharges",
					"Number of charges", "Number of charges", Parameter::SPINBOX, true,
					component.numCharges()));
	parameters.push_back(
			new ParameterWithIntValue(name + ".NumberOfDipoles",
					"Number of dipoles", "Number of Dipoles",Parameter::SPINBOX, true,
					component.numDipoles()));
	parameters.push_back(
			new ParameterWithIntValue(name + ".NumberOfQuadrupoles","Number of Quadrupoles",
					"Number of Quadrupoles",Parameter::SPINBOX,true, component.numQuadrupoles()));

	for (size_t i = 0; i < component.numCharges(); i++) {
		std::stringstream baseNameStream;
		baseNameStream << name << ".Charge" << i;
		std::string baseName = baseNameStream.str();
		ParameterCollection* chargeCollection = new ParameterCollection(
				baseName, baseName, baseName, Parameter::BUTTON);
		Charge charge = component.charge(i);
		chargeCollection->addParameter(
				new ParameterWithDoubleValue(baseName + ".x", baseName + ".x [Angstrom]",
						baseName + ".x",Parameter::LINE_EDIT, false, charge.rx() / MDGenerator::angstroem_2_atomicUnitLength));
		chargeCollection->addParameter(
				new ParameterWithDoubleValue(baseName + ".y", baseName + ".y [Angstrom]",
						baseName + ".y",Parameter::LINE_EDIT, false, charge.ry() / MDGenerator::angstroem_2_atomicUnitLength));
		chargeCollection->addParameter(
				new ParameterWithDoubleValue(baseName + ".z", baseName + ".z [Angstrom]",
						baseName + ".z",Parameter::LINE_EDIT, false, charge.rz() / MDGenerator::angstroem_2_atomicUnitLength));
		chargeCollection->addParameter(
				new ParameterWithDoubleValue(baseName + ".mass",
						baseName + ".mass [u]", baseName + ".mass",Parameter::LINE_EDIT, false,
						charge.m() / MDGenerator::unitMass_2_mardyn));
		chargeCollection->addParameter(
				new ParameterWithDoubleValue(baseName + ".charge", "charge",
						baseName + ".charge [el. charge e]", Parameter::LINE_EDIT,false, charge.q()));
		addParameter(chargeCollection);
	}
	for (unsigned int i = 0; i < component.numLJcenters(); i++) {
		std::stringstream baseNameStream;
		baseNameStream << name << ".LJCenter" << i;
		std::string baseName = baseNameStream.str();
		ParameterCollection* ljCenterCollection = new ParameterCollection(
				baseName, baseName, baseName, Parameter::BUTTON);
		LJcenter ljCenter = component.ljcenter(i);
		ljCenterCollection->addParameter(
				new ParameterWithDoubleValue(baseName + ".x", baseName + ".x [Angstrom]",
						baseName + ".x", Parameter::LINE_EDIT,false, ljCenter.rx() / MDGenerator::angstroem_2_atomicUnitLength));
		ljCenterCollection->addParameter(
				new ParameterWithDoubleValue(baseName + ".y", baseName + ".y [Angstrom]",
						baseName + ".y",Parameter::LINE_EDIT, false, ljCenter.ry() / MDGenerator::angstroem_2_atomicUnitLength));
		ljCenterCollection->addParameter(
				new ParameterWithDoubleValue(baseName + ".z", baseName + ".z [Angstrom]",
						baseName + ".z",Parameter::LINE_EDIT, false, ljCenter.rz() / MDGenerator::angstroem_2_atomicUnitLength));
		ljCenterCollection->addParameter(
				new ParameterWithDoubleValue(baseName + ".mass",
						baseName + ".mass [u]", baseName + ".mass",Parameter::LINE_EDIT, false,
						ljCenter.m() / MDGenerator::unitMass_2_mardyn));
		ljCenterCollection->addParameter(
				new ParameterWithDoubleValue(baseName + ".epsilon",
						baseName + ".epsilon [K]", "Epsilon normalized to the Boltzmann constant", Parameter::LINE_EDIT,false,
						ljCenter.eps() / MDGenerator::kelvin_2_mardyn));
		ljCenterCollection-> addParameter(
				new ParameterWithDoubleValue(baseName + ".sigma",
						baseName + ".sigma [Angstrom]", baseName + ".sigma", Parameter::LINE_EDIT,false,
						ljCenter.sigma() / MDGenerator::angstroem_2_atomicUnitLength));
		addParameter(ljCenterCollection);
	}

	for (unsigned int i = 0; i < component.numDipoles(); i++) {
		std::stringstream baseNameStream;
		baseNameStream << name << ".Dipole" << i;
		std::string baseName = baseNameStream.str();
		ParameterCollection* dipoleCollection = new ParameterCollection(
				baseName, baseName, baseName, Parameter::BUTTON);
		Dipole dipole = component.dipole(i);
		dipoleCollection->addParameter(new ParameterWithDoubleValue(baseName + ".x", baseName + ".x [Angstrom]", baseName + ".x",
				Parameter::LINE_EDIT,false, dipole.rx() / MDGenerator::angstroem_2_atomicUnitLength));
		dipoleCollection->addParameter(new ParameterWithDoubleValue(baseName + ".y", baseName + ".y [Angstrom]", baseName + ".y",
				Parameter::LINE_EDIT,false, dipole.ry() / MDGenerator::angstroem_2_atomicUnitLength));
		dipoleCollection->addParameter(new ParameterWithDoubleValue(baseName + ".z", baseName + ".z [Angstrom]", baseName + ".z",
				Parameter::LINE_EDIT,false, dipole.rz() / MDGenerator::angstroem_2_atomicUnitLength));
		dipoleCollection->addParameter(new ParameterWithDoubleValue(baseName + ".eMyx", baseName + ".eMyx", baseName + ".eMyx",
				Parameter::LINE_EDIT,false, dipole.ex()));
		dipoleCollection->addParameter(new ParameterWithDoubleValue(baseName + ".eMyy", baseName + ".eMyy", baseName + ".eMyy",
				Parameter::LINE_EDIT,false, dipole.ey()));
		dipoleCollection->addParameter(new ParameterWithDoubleValue(baseName + ".eMyz", baseName + ".eMyz", baseName + ".eMyz",
				Parameter::LINE_EDIT,false, dipole.ez()));
		dipoleCollection->addParameter(new ParameterWithDoubleValue(baseName + ".dipole", baseName + ".dipole [Debye]", baseName + ".dipole",
				Parameter::LINE_EDIT,false, dipole.absMy() / MDGenerator::debye_2_mardyn));
		addParameter(dipoleCollection);

	}
	for (unsigned int i = 0; i < component.numQuadrupoles(); i++) {
		std::stringstream baseNameStream;
		baseNameStream << name << ".Quadrupole" << i;
		std::string baseName = baseNameStream.str();
		Quadrupole quad = component.quadrupole(i);
		ParameterCollection* quadrupoleCollection = new ParameterCollection(baseName, baseName, baseName, Parameter::BUTTON);
		quadrupoleCollection->addParameter(new ParameterWithDoubleValue(baseName + ".x", baseName + ".x [Angstrom]", baseName + ".x",
				Parameter::LINE_EDIT,false, quad.rx() / MDGenerator::angstroem_2_atomicUnitLength));
		quadrupoleCollection->addParameter(new ParameterWithDoubleValue(baseName + ".y", baseName + ".y [Angstrom]", baseName + ".y",
				Parameter::LINE_EDIT,false, quad.ry() / MDGenerator::angstroem_2_atomicUnitLength));
		quadrupoleCollection->addParameter(new ParameterWithDoubleValue(baseName + ".z", baseName + ".z [Angstrom]", baseName + ".z",
				Parameter::LINE_EDIT,false, quad.rz() / MDGenerator::angstroem_2_atomicUnitLength));
		quadrupoleCollection->addParameter(new ParameterWithDoubleValue(baseName + ".eQx", baseName + ".eQx", baseName + ".eQx",
				Parameter::LINE_EDIT,false, quad.ex()));
		quadrupoleCollection->addParameter(new ParameterWithDoubleValue(baseName + ".eQy", baseName + ".eQy", baseName + ".eQy",
				Parameter::LINE_EDIT,false, quad.ey()));
		quadrupoleCollection->addParameter(new ParameterWithDoubleValue(baseName + ".eQz", baseName + ".eQz", baseName + ".eQz",
				Parameter::LINE_EDIT,false, quad.ez()));
		quadrupoleCollection->addParameter(new ParameterWithDoubleValue(baseName + ".quadrupole", baseName + ".quadrupole [Buckingham]", baseName + ".quadrupole",
				Parameter::LINE_EDIT,false, quad.absQ() / MDGenerator::buckingham_2_mardyn));
		addParameter(quadrupoleCollection);
	}

}

ComponentParameters::~ComponentParameters() {
}

void ComponentParameters::setParameterValue(Component& component, Parameter* p,
		std::string& valueName) {
#ifdef DEBUG
	std::cout
			<< "ComponentParameters::setParameterValue() from Parameter Object to component: name is "
			<< valueName << std::endl;
#endif

	if (valueName == "NumberOfCharges") {
		size_t numCharges =
				dynamic_cast<const ParameterWithIntValue*> (p)->getValue();
		while (component.numCharges() < numCharges) {
			component.addCharge(0.0, 0.0, 0.0, 0.0, 0);
		}
		while (component.numCharges() > numCharges) {
			component.deleteCharge();
		}
	} else if (valueName == "NumberOfLJCenters") {
		size_t numLJCenters =
				dynamic_cast<const ParameterWithIntValue*> (p)->getValue();
		while (component.numLJcenters() < numLJCenters) {
			component.addLJcenter(0, 0, 0, 0, 0, 0, 0, 0);
		}
		while (component.numLJcenters() > numLJCenters) {
			component.deleteLJCenter();
		}
	} else if (valueName == "NumberOfDipoles"){
		size_t numDipoles = static_cast<const ParameterWithIntValue*> (p)->getValue();
		while (component.numDipoles() < numDipoles){
			component.addDipole(0,0,0,0,0,0,0);
		}

		while (component.numDipoles() > numDipoles) {
			component.deleteDipole();
		}
 	} else if (valueName == "NumberOfQuadrupoles") {
 		size_t numQuadrupoles = static_cast<const ParameterWithIntValue*> (p)->getValue();
 		while (component.numQuadrupoles() < numQuadrupoles) {
 			component.addQuadrupole(0,0,0,0,0,0,0);
 		}
 		while (component.numQuadrupoles() > numQuadrupoles) {
 			component.deleteQuadrupole();
 		}
 	} else if (valueName == "Temperature") {
 		component.setT(static_cast<ParameterWithDoubleValue*> (p)->getValue());
#ifdef DEBUG
 		std::cout<<"Temperature set to "<<static_cast<ParameterWithDoubleValue*> (p)->getValue()<<std::endl;
#endif
 	} else if (valueName.find("Charge") != std::string::npos) {
		int chargeIndex = convertFromChar<int> (valueName.at(6));
		std::string part = valueName.substr(8);
#ifdef DEBUG
		std::cout << "Index of charge: " << chargeIndex << std::endl;
		std::cout << "part is: " << part << std::endl;
#endif
		setParameterValue(component.charge(chargeIndex),
				static_cast<ParameterWithDoubleValue*> (p), part);
	} else if (valueName.find("LJCenter") != std::string::npos) {
		int ljCenterIndex = convertFromChar<int> (valueName.at(8));
		std::string part = valueName.substr(10);
#ifdef DEBUG
		std::cout << "Index of LJ Center: " << ljCenterIndex << std::endl;
		std::cout << "part is: " << part << std::endl;
#endif
		setParameterValue(component.ljcenter(ljCenterIndex),
				static_cast<ParameterWithDoubleValue*> (p), part);
	} else if (valueName.find("Dipole") != std::string::npos) {
		int dipoleIndex = convertFromChar<int> (valueName.at(6));
		std::string part = valueName.substr(8);
#ifdef DEBUG
		std::cout << "Index of Dipole: " << dipoleIndex << std::endl;
		std::cout << "part is: " << part << std::endl;
#endif
		setParameterValue(component.dipole(dipoleIndex), static_cast<ParameterWithDoubleValue*> (p), part);
	} else if (valueName.find("Quadrupole") != std::string::npos){
		int quadIndex = convertFromChar<int> (valueName.at(10));
		std::string part = valueName.substr(12);
#ifdef DEBUG
		std::cout<< "Index of Quadrupole: " << quadIndex << std::endl;
		std::cout << "part is: " << part << std::endl;
#endif
		setParameterValue(component.quadrupole(quadIndex), static_cast<ParameterWithDoubleValue*> (p), part);
	} else {
		std::cout << "ComponentParameters::setParameterValue(): unkown parameter! (valueName="
				<< valueName << ")" << std::endl;
		exit(-1);
	}
}

void ComponentParameters::setParameterValue(Charge& charge,
		ParameterWithDoubleValue* p, std::string& valueName) {
#ifdef DEBUG
	std::cout << "SetParameterValue in Charge! valueName=" << valueName
			<< "value: " << p->getValue() << std::endl;
#endif

	if (valueName == "x") {
		charge.setR(0, p->getValue() * MDGenerator::angstroem_2_atomicUnitLength);
	} else if (valueName == "y") {
		charge.setR(1, p->getValue() * MDGenerator::angstroem_2_atomicUnitLength);
	} else if (valueName == "z") {
		charge.setR(2, p->getValue()  * MDGenerator::angstroem_2_atomicUnitLength);
	} else if (valueName == "mass") {
		charge.setM(p->getValue() * MDGenerator::unitMass_2_mardyn);
	} else if (valueName == "charge") {
		charge.setQ(p->getValue());
	} else {
		std::cout << "ComponentParameters::setParameterValue(Charge): unkown parameter! (valueName="
				<< valueName << ")" << std::endl;
		exit(-1);
	}

}

void ComponentParameters::setParameterValue(LJcenter& ljCenter,
		ParameterWithDoubleValue* p, std::string& valueName) {
#ifdef DEBUG
	std::cout << "SetParameterValue in LJCenter! valueName =" << valueName
			<< "value: " << p->getValue() << std::endl;
#endif

	if (valueName == "x") {
		ljCenter.setR(0, p->getValue() * MDGenerator::angstroem_2_atomicUnitLength);
	} else if (valueName == "y") {
		ljCenter.setR(1, p->getValue() * MDGenerator::angstroem_2_atomicUnitLength);
	} else if (valueName == "z") {
		ljCenter.setR(2, p->getValue() * MDGenerator::angstroem_2_atomicUnitLength);
	} else if (valueName == "mass") {
		ljCenter.setM(p->getValue()  * MDGenerator::unitMass_2_mardyn);
	} else if (valueName == "epsilon") {
		ljCenter.setEps(p->getValue() * MDGenerator::kelvin_2_mardyn);
	} else if (valueName == "sigma") {
		ljCenter.setSigma(p->getValue() * MDGenerator::angstroem_2_atomicUnitLength);
	} else {
		std::cout << "ComponentParameters::setParameterValue(LJCenter): unkown parameter! (valueName="
				<< valueName << ")" << std::endl;
		exit(-1);
	}
}

void ComponentParameters::setParameterValue(Dipole& dipole, ParameterWithDoubleValue* p, std::string& valueName){
#ifdef DEBUG
	std::cout << "SetParameterValue in Dipole! valueName =" << valueName << "value: " << p->getValue() << std::endl;
#endif

	if (valueName == "x") {
		dipole.setR(0,p->getValue() * MDGenerator::angstroem_2_atomicUnitLength);
	} else if (valueName == "y") {
		dipole.setR(1,p->getValue() * MDGenerator::angstroem_2_atomicUnitLength);
	} else if (valueName == "z") {
		dipole.setR(2,p->getValue() * MDGenerator::angstroem_2_atomicUnitLength);
	} else if (valueName == "eMyx") {
		dipole.setE(0,p->getValue());
	} else if (valueName == "eMyy") {
		dipole.setE(1,p->getValue());
	} else if (valueName == "eMyz") {
		dipole.setE(2,p->getValue());
	} else if (valueName == "dipole") {
		dipole.setAbyMy(p->getValue() * MDGenerator::debye_2_mardyn);
	} else {
		std::cout << "ComponentParameters::setParameterValue(Dipole): unkown parameter! (valueName="
				<< valueName << ")" << std::endl;
		exit(-1);
	}
}

void ComponentParameters::setParameterValue(Quadrupole& quadrupole, ParameterWithDoubleValue* p, std::string& valueName){
#ifdef DEBUG
	std::cout << "SetParameterValue in Quadrupole! valueName =" << valueName << "value: " << p->getValue() << std::endl;
#endif

	if (valueName == "x") {
		quadrupole.setR(0,p->getValue() * MDGenerator::angstroem_2_atomicUnitLength);
	} else if (valueName == "y") {
		quadrupole.setR(1,p->getValue() * MDGenerator::angstroem_2_atomicUnitLength);
	} else if (valueName == "z") {
		quadrupole.setR(2,p->getValue() * MDGenerator::angstroem_2_atomicUnitLength);
	} else if (valueName == "eQx") {
		quadrupole.setE(0,p->getValue());
	} else if (valueName == "eQy") {
		quadrupole.setE(1,p->getValue());
	} else if (valueName == "eQz") {
		quadrupole.setE(2,p->getValue());
	} else if (valueName == "quadrupole") {
		quadrupole.setAbsQ(p->getValue() * MDGenerator::buckingham_2_mardyn);
	} else {
		std::cout << "ComponentParameters::setParameterValue(Quadrupole): unkown parameter! (valueName="
				<< valueName << ")" << std::endl;
		exit(-1);
	}
}

void ComponentParameters::load(const std::string& filename, Generator* generator) {
	PMFileReader::readPMFile(filename, generator, this);
}
