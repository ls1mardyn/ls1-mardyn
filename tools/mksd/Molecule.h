/*
 * Molecule.h
 *
 *  Created on: 09.01.2012
 *      Author: Stefan Becker <stefan.becker@mv.uni-kl.de>
 *
 * The instance of this class writes the body of the phase space file.
 * It receives the component specific information (from the object "Component")
 * and the position of every single particle (from the object "StartGeometry"). The body of the phase space file contains:
 * 		 - molecule number (_moleculeId)
 * 		 - component type (_componentId)
 * 		 - position
 * 		 - velocity
 * 		 - quaterion oreintation
 * 		 - angular momentum
 */

#ifndef MOLECULE_H_
#define MOLECULE_H_

#include<iostream>
#include<cmath>
#include"RandomNumber.h"
#include<map>

extern const double PI;


class Molecule
{
public:
	// constructor and methods
	Molecule(double compT, double mass);
	// the get-methods, self-explanatory
	unsigned gMoleculeCID(unsigned moleculeCounter);
	unsigned gComponentId();
	double gXPos(unsigned id);
	double gYPos(unsigned id);
	double gZPos(unsigned id);
	double gXVelocity(unsigned id);
	double gYVelocity(unsigned id);
	double gZVelocity(unsigned id);
	double gNumberOfMolecules();
	//Quaternion gOrientationQuaternion; 	// @todo: to be implemented when needed.
	//double Molecule::gAngularMomentum();// @todo: to be implemented when needed.

	// the set-methods
	void sIdNumber(unsigned MoleculeId);
	void sComponentId(unsigned cId);
	void sPosition(double pos[3]);
	void sVelocity(double vel[3]);
	void sOrientationQuaternion(double quatOrient[4]);	//@todo: to be implemented when needed
	void sAngularMomentum(double angmom[3]);			//@todo: to be implemented when needed

	/*void calculateCoordinatesOfWallMolecule(unsigned numberOfLayers,double xLength, double yLength, double zLength,
											double off0, double off1, double off2, double latticeConstant);*/
	void calculateCoordinatesOfWallMolecule(double xLength, double zLength, double off1,
											double latticeConstant, double shielding);
	void calculateCoordinatesOfWallMolecule(double xLength, double zLength, double off1,
													  double latticeConstant, double shielding,
													  unsigned numberOfStripes);
	void calculateVelocities();

	// the operator writing to the stream, i.e. the phase space file.
//	ostream& operator << (ostrean& ofstr, WallMolecule wallMol);

private:
	double _componentTemperature;
	double _componentMass;
//	unsigned _moleculeId;
	unsigned _numberOfMolecules;
	//unsigned _componentId;
	//double _position[3];
	//double _velocity[3];
	//Quaternion orientationQuaternion; // orientation in quaternions (q0,q1,q2,q3), where q = q0 + i*q1+j*q2 + k*q3 => class to be implemented
	//double _angularMomentum[3];

	// @brief: the key element of the map is always the internal (within this class only) molecule-ID
	// the value element of the map is denoted by second part of the name (e.g. position, velocity, component-ID, etc.) and corresponds to the
	// particular molecule (of any component)
	std::map<unsigned, double> _moleculePosition[3];  // for the following three maps: the key value starts at 0.
	std::map<unsigned, double> _moleculeVelocity[3];
	std::map<unsigned, unsigned> _moleculeCID;
	//map<unsigned,unsigned> _moleculeComponentID;

};

#endif /* MOLECULE_H_ */
