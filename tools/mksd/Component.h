/*
 * Component.h
 *
 *  Created on: 09.01.2012
 *      Author: Stefan Becker <stefan.becker@mv.uni-kl.de>
 *
 * This class contains the models of the substances to be simulated, e.g. sigma, eps, Tersoff-parameters, etc.
 */

#ifndef COMPONENT_H_
#define COMPONENT_H_

#include<iostream> // for cerr <<
#include<cstdlib> // for exit()
#include<string>
#include<vector>
#include<cmath>  // for pow()


class Component{
private:
	std::string _substance;
	unsigned _numberLJCenters;
	unsigned _numberCharges;
	unsigned _numberQuadrupoles;
	unsigned _numberDipoles;
	unsigned _numberTersoff;
	// LJ centers
	// LJ center position, vector quantities to account for multicenter LJ potentials
	std::vector<double> _vecLJx;
	std::vector<double> _vecLJy;
	std::vector<double> _vecLJz;
	// mass of each LJ center
	std::vector<double> _vecLJMass;
	// energy of each LJ center => eps
	std::vector<double> _vecLJEps;
	// length of each LJ center =>sigma
	std::vector<double> _vecLJSigma;
	// Lennard-Jones cut off radius
	double _TSLJCutOff;



	// moment of inertia
	double _ixx;
	double _iyy;
	double _izz;

	// reference quantities
	double _refLength;
	double _refEnergy;
	double _refMass;
	//double _refTime;


public:
	// the constructor and destructor, respectively
	// formal arguments of the constructor: substance, reference energy, reference length, reference mass
	//Component(std::string substance);
	// second Constructor, in case reference quantities different from atomic units are allowed
//	Component(std::string in_substance, double refNRG, double refLgth, double refM);
	// another constructor, called if the reference quantities are atomic units
	Component(std::string in_substance, bool in_LJunits);
	~Component();

	// the methods
	// !!!@todo: if refence quantities different from atomic units (e.g. LJ units) are allowed
	// most of the methods have to be adjusted!
	// (i) get methods
	unsigned gNumberLJCenters();
	unsigned gNumberCharges();
	unsigned gNumberQuadrupoles();
	unsigned gNumberDipoles();
	unsigned gNumberTersoff();

	//@todo: @implementation: so far, the first
	double gSigma(unsigned i);
	double gEps(unsigned i);
	double gMass(unsigned i);
	double gRCutLJ();

	// for use in unit reduction by LJ units => the reference length corresponds to the smallest sigma
	double gSigmaMin();
	double gRefTime(bool in_LJunits);

	// (ii) the set methods
	/*void sNumberLJCenters(unsigned nCLJ);
	void sNumberCharges(unsigned NChg);
	void sNumberQuadrupoles (unsigned NQdpl);
	void sNumberDipoles(unsigned NDpl);
	void sNumberTersoff(unsigned NTers);*/

	// (ii) initializing a 1CLJ fluid, parameters depend on the specific fluid (e.g. Argon, CH4, etc.)
	void init1CLJ(std::string substance);

	// calculating the liquid and vapor densities, respectively
	double calculateLiquidDensity(double T);
	double calculateVaporDensity(double T, double factor);


};


#endif /* COMPONENT_H_ */
