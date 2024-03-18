/*
 * Component.cpp
 *
 *  Created on: 09.01.2012
 *      Author: becker
 */

#include "Component.h"

// fluids
extern const std::string FLUID_AR = "Ar";
extern const std::string FLUID_CH4 = "CH4";

const double EPS_AR = 4.36704e-04;
const double SIGMA_AR = 6.40920;
const double AR_MASS = 0.03994;

const double EPS_CH4 = 5.54383e-04;
const double SIGMA_CH4 = 7.03753;
const double CH4_MASS = 0.016042;

extern const double RHO_CRITICAL_1CLJ = 0.3190;
const double T_CRITICAL_1CLJ = 1.0779;



// solids
extern const std::string WALL_CU_LJ = "CU_LJ"; // copper by Lennard-Jones
extern const std::string WALL_TERSOFF = "Tersoff";
extern const std::string WALL_MAT_AKA_FIT = "MatAkaFit";

// According to Phillpot, S. R. in: "Reconstruction of grain boundaries in copper and gold by simulation"
// Journal of Materials Ressearch, Vol 9, No.3 1994;
//const double EPS_CU = 6.13686e-03; // not used
const double EPS_CU = 4.36704e-02;  // a hundred times the fluid epsilon (in this case of argon)
const double SIGMA_CU = SIGMA_AR;//4.37472;
const double CU_MASS = 0.063546;
// for SIGMA_CU = SIGMA_AR
// for SIGMA_CU = 0.8*SIGMA_AR
double LATTICE_CONST_WALL_LJTS;



Component::Component(std::string in_substance, bool in_LJunits){
	_substance = in_substance;
	double facLatConst = 1.0; //0.79852 for sigma_ss = 0.8*sigma_ff and 1.0 for sigma_ss = sigma_ff
	if (in_LJunits){
	_refEnergy = EPS_AR;
	_refLength = SIGMA_AR;
	_refMass = AR_MASS;
	LATTICE_CONST_WALL_LJTS = facLatConst *9.92234/ SIGMA_AR;
	}
	else{	// i.e. atomic units
		_refEnergy = 1.0;
		_refLength = 1.0;
		_refMass = 1.0;
		LATTICE_CONST_WALL_LJTS = facLatConst * 9.92234;
	}
	//_refTime = _refLength*sqrt(_refMass/_refEnergy);

	if( _substance == FLUID_AR || _substance == FLUID_CH4 || _substance == WALL_CU_LJ){
		init1CLJ(_substance);
	}
	else{
		std::cerr << "No other interaction models implemented.";
	}
}

Component::~Component(){

}

void Component::init1CLJ(std::string substance){

			// parameters all the 1C LJ fluid models have in common
			_numberLJCenters = 1;
			_numberCharges = 0;
			_numberQuadrupoles = 0;
			_numberDipoles = 0;
			_numberTersoff = 0;
			// Vektoren auf passende Länge bringen (Anzahl Vektorelemente, in diesem Fall: 2 Elemente. Eines pro LJ-center)
			// alternativ: Bsp.: _vecLJx(2,0) => 2-elementiger Vektor, initialisiert mit "0", i.e. _vecLJx == [0,0]
			_vecLJx.resize(_numberLJCenters);
			_vecLJy.resize(_numberLJCenters);
			_vecLJz.resize(_numberLJCenters);
			_vecLJMass.resize(_numberLJCenters);
			_vecLJEps.resize(_numberLJCenters);
			_vecLJSigma.resize(_numberLJCenters);
			// setting the model parameters
			_vecLJx[0] = 0;   	// x-position of the LJ-center, with respect to the local particle fixed coordinate system
			_vecLJy[0] = 0;		// y-position ...
			_vecLJz[0] = 0;		// z-position ...
			// moments of inertia
			_ixx = 0;
			_iyy = 0;
			_izz = 0;
			// specific fluid parameters mass, sigma, epsilon
			if(substance == FLUID_AR)
			{
				_vecLJMass[0] = AR_MASS;
				_vecLJEps[0] = EPS_AR;
				_vecLJSigma[0] = SIGMA_AR;
				_TSLJCutOff = 2.5*SIGMA_AR;
			}
			else if(substance == FLUID_CH4)
			{
				_vecLJMass[0] = CH4_MASS;
				_vecLJEps[0] = EPS_CH4;
				_vecLJSigma[0] = SIGMA_CH4;
				_TSLJCutOff = 2.5*SIGMA_CH4;
			}
			else if(substance == WALL_CU_LJ)
			{
				_vecLJMass[0] = CU_MASS;
				_vecLJEps[0] = EPS_CU;
				//_vecLJSigma[0] = SIGMA_AR;
				_vecLJSigma[0] = SIGMA_CU;
				_TSLJCutOff = 2.5*SIGMA_AR;
			}
}

unsigned Component::gNumberLJCenters(){
	return _numberLJCenters;
}

unsigned Component::gNumberCharges(){
	return _numberCharges;
}

unsigned Component::gNumberQuadrupoles(){
	return _numberQuadrupoles;
}

unsigned Component::gNumberDipoles(){
	return _numberDipoles;
}

unsigned Component::gNumberTersoff(){
	return _numberTersoff;
}

//@brief: multicenter LJ, method returns the smallest sigma; to be used for unit reduction
double Component::gSigmaMin()
{
	double temp;
	temp = _vecLJSigma[0];
	for(unsigned i = 0; i < _vecLJSigma.size(); i++)
	{
		if(temp > _vecLJSigma[i])
		{
			temp = _vecLJSigma[i];
		}
	}
	return temp/_refLength;
}

double Component::gSigma(unsigned i){
	return _vecLJSigma.at(i)/_refLength;
}

double Component::gEps(unsigned i){
	return _vecLJEps.at(i)/_refEnergy;
}

double Component::gMass(unsigned i){
	return _vecLJMass.at(i)/_refMass;
}

double Component::gRCutLJ(){
	return _TSLJCutOff/_refLength;
}

double Component::calculateLiquidDensity(double T){
	double rhoLiq;
	if(_substance == FLUID_AR || _substance == FLUID_CH4){
		// @brief: calculation of the bulk densities in the liquid and vapor phase of a 1CLJ fluid, respectively;
		// 		   according to Kedia et al. in:  Molecular Physics, vol. 104, Issue 9, p.1509-1527
		rhoLiq =(RHO_CRITICAL_1CLJ + 0.5649* pow(T_CRITICAL_1CLJ - T, (1.0/3.0))
				+ 0.1314 * (T_CRITICAL_1CLJ - T) + 0.0413 * pow(T_CRITICAL_1CLJ - T, 1.5) )
				/ (gSigma(0)*gSigma(0)*gSigma(0));
	}
	else{
		std::cerr << "Error in Component class: Claculation of the liquid density. Liquid density of the 1C LJ model not calculated!";
		exit(-201);
	}
	//cout << "Calculating rhoLiq in Component:\ngSigma(0) = " << (gSigma(0)) << "\n"<< "temperature = "<< T <<"\n";
	std::cout << "rhoLiq = " << rhoLiq << "\n";
	return rhoLiq;
}

double Component::calculateVaporDensity(double T, double factor){
	double rhoVap;
	if(_substance == FLUID_AR || _substance == FLUID_CH4){
		// @brief: calculation of the bulk densities in the liquid and vapor phase of a 1CLJ fluid, respectively;
		// 		   according to Kedia et al. in:  Molecular Physics, vol. 104, Issue 9, p.1509-1527
		rhoVap = factor * ( RHO_CRITICAL_1CLJ - 0.5649* pow(T_CRITICAL_1CLJ - T, 1.0/3.0)
				  + 0.2128 * (T_CRITICAL_1CLJ - T) + 0.0702 * pow(T_CRITICAL_1CLJ - T, 1.5) )
				 / (gSigma(0)*gSigma(0)*gSigma(0));
	}
	else{
			std::cerr << "Error in Component class: Claculation of the liquid density. Liquid density of the 1C LJ model not calculated!";
			exit(-202);
		}
	//cout << "Calculating rhoVap in Component: \ngSigma(0) = " << (gSigma(0)) << "\n"<< "temperature = "<< T <<"\n";
	std::cout << "rhoVap = " << rhoVap << "\n";
	return rhoVap;
}


double Component::gRefTime(bool in_LJunits){
	double refTime;
	if(in_LJunits){
		refTime = _refLength *sqrt(_refMass/_refEnergy);
	}
	else{
		refTime = 1.0;
	}
	return refTime;
}
