/*
 * PhaseSpaceWriter.h
 *
 *  Created on: 21.01.2012
 *      Author: Stefan Becker <stefan.becker@mv.uni-kl.de>
 */

#ifndef PHASESPACEWRITER_H_
#define PHASESPACEWRITER_H_


#include<string>
#include<iostream>
#include<fstream>

#include"GlobalStartGeometry.h"
#include"Molecule.h"
#include"Component.h"
#include"RandomNumber.h"

using namespace std;


class PhaseSpaceWriter{

public:
	// constructor + destructor
	PhaseSpaceWriter(string in_prefix, double in_Temperature, double in_densFac, unsigned in_nFluid, string in_fluidComponent,
					string in_wallComponent, unsigned in_wallLayers, double in_xi12, double in_xi13, double in_eta,  double in_alpha, double in_beta, double in_gamma, double in_edgeProp, bool in_stripes,
					unsigned in_numberOfStripes, bool in_LJShifted, bool in_LJunits);

	~PhaseSpaceWriter();
	// methods
	void write();
	double gBoxLengthY();
	double gTemperature();
	double gAverageMassPerParticle();
	unsigned gNTotal();



private:

	string _fileName;
	string _fileNameXyz;
	string _fluidComponentName;
	string _wallComponentName;
	unsigned _nFluid;
	unsigned _nTotal;
	unsigned _wallLayers;
	unsigned _numberOfStripes;
	double _temperature;
	double _densFac;
	double _xi12, _xi13;
	double _eta12;
	double _alpha;
	double _beta;
	double _gamma;
	double _edgeProp;
	double _boxLengthY;
	double _avMass;
	bool _stripes;
	bool _LJShifted;
	bool _LJunits;

	// ofstream: phase space stream
//	ofstream _psstrm;
	// stringstream
	//stringstream _strstrm;

};

#endif /* PHASESPACEWRITER_H_ */
