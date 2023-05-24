/*
 * Molecule.cpp
 *
 *  Created on: 09.01.2012
 *      Author: becker
 */
#include "Molecule.h"
Molecule::Molecule(double compT, double mass){
_componentTemperature = compT;
_componentMass = mass;
_numberOfMolecules = 0;
}


unsigned Molecule::gMoleculeCID(unsigned moleculeCounter){
	return _moleculeCID[moleculeCounter];
}

double Molecule::gXPos(unsigned id){
	return _moleculePosition[0][id];
}

double Molecule::gYPos(unsigned id){
	return _moleculePosition[1][id];
}

double Molecule::gZPos(unsigned id){
	return _moleculePosition[2][id];
}

void Molecule::calculateVelocities(){
	double phi, omega;
	double absoluteVelocity = sqrt(3.0 * _componentTemperature / _componentMass);
	RandomNumber rdm;
	for(unsigned i = 0; i < _numberOfMolecules; i++){
		phi = 2*PI*rdm.randNum();
		omega = 2*PI*rdm.randNum();
		_moleculeVelocity[0][i] = absoluteVelocity*cos(phi)*cos(omega);
		_moleculeVelocity[1][i] = absoluteVelocity*cos(phi)*sin(omega);
		_moleculeVelocity[2][i] = absoluteVelocity*sin(phi);
	}
}

double Molecule::gXVelocity(unsigned id){
	return _moleculeVelocity[0][id];
}

double Molecule::gYVelocity(unsigned id){
	return _moleculeVelocity[1][id];
}

double Molecule::gZVelocity(unsigned id){
	return _moleculeVelocity[2][id];
}

double Molecule::gNumberOfMolecules(){
	return _numberOfMolecules;
}


// used if there's a single set of fluid wall interaction parameters, in contrast to the stripes shaped wall, e.g.
void Molecule::calculateCoordinatesOfWallMolecule(double xLength, double zLength, double off1,
												  double latticeConstant, double shielding){
	double yLength = off1 - shielding;
	double pos[3][3]; 	// 3 vectors which are added to the current position in order to set three new particles
	for(unsigned i = 0; i<3; i++){    // !!!this initialisation of the 3 vectors is valid only for equilateral lattices!!! => see sketch!
		for(unsigned j= 0; j<3; j++){
			if(i!= j) pos[i][j] = 0.5*latticeConstant;
			else pos[i][j] = 0;
			//cout << "pos["<<i<<"]["<< j<<"] = " << pos[i][j] << "\n";
		}
	}
	_numberOfMolecules = 0;
	double currentZeroPos[3];
	currentZeroPos[0] = -0.95*latticeConstant; // 1-0.95 = 0.05 => offset so that all the particles are within the box!

	for(unsigned i = 0; i*latticeConstant < xLength; i++){
		currentZeroPos[0] = currentZeroPos[0] +latticeConstant;
		currentZeroPos[1] = off1-shielding + 1.05*latticeConstant;  // little offset so that there is no risk that particles are placed out of the box (i.e. y-coordinate < 0)
		currentZeroPos[2] = -0.95*latticeConstant;
		for(unsigned j = 0; j*latticeConstant < yLength; j++){
			currentZeroPos[1] = currentZeroPos[1] - latticeConstant; // minus because the start is at the top of the wall going to the bottom of the box afterwards
			currentZeroPos[2] =  -0.95*latticeConstant; //_moleculePosition[2][1] - latticeConstant;
			for(unsigned k = 0; k*latticeConstant < zLength; k++){
				currentZeroPos[2] = currentZeroPos[2] +latticeConstant;
				_moleculePosition[0][_numberOfMolecules] = currentZeroPos[0];
				_moleculePosition[1][_numberOfMolecules] = currentZeroPos[1];
				_moleculePosition[2][_numberOfMolecules] = currentZeroPos[2];
				_moleculeCID[_numberOfMolecules] = 2;
				_numberOfMolecules++;
				for(unsigned l = 0; l < 3; l++){
					_moleculePosition[0][_numberOfMolecules] = currentZeroPos[0] + pos[l][0];
					_moleculePosition[1][_numberOfMolecules] = currentZeroPos[1] - pos[l][1];
					_moleculePosition[2][_numberOfMolecules] = currentZeroPos[2] + pos[l][2];
					_moleculeCID[_numberOfMolecules] = 2;
					_numberOfMolecules++;
				} // end for(l)
			}// end for (k)
		} // end for(j)
	} // end for(i)

} // end of method


// overloaded form of the above method: used if there's a stripes shaped wall with 2 different sets of fluid wall interaction parameters.
void Molecule::calculateCoordinatesOfWallMolecule(double xLength, double zLength, double off1,
												  double latticeConstant, double shielding, unsigned numberOfStripes){
	double yLength = off1 - shielding;
	double pos[3][3]; 	// 3 vectors which are added to the current position in order to set three new particles
	for(unsigned i = 0; i<3; i++){    // !!!this initialisation of the 3 vectors is valid only for equilateral lattices!!! => see sketch!
		for(unsigned j= 0; j<3; j++){
			if(i!= j) pos[i][j] = 0.5*latticeConstant;
			else pos[i][j] = 0;
			//cout << "pos["<<i<<"]["<< j<<"] = " << pos[i][j] << "\n";
		}
	}
	_numberOfMolecules = 0;
	double stripeWidth = xLength / numberOfStripes; // width of a single stripe in x-direction => delta_x
	if (stripeWidth < 0.5*latticeConstant){
		std::cerr << "Too many stripes: width would be smaller than half the lattice constant! Meaningless => width adapted to half the lattice constant!";
		stripeWidth = 0.5 * latticeConstant;
	}
	double currentZeroPos[3];
	currentZeroPos[0] = -0.95*latticeConstant; // 1-0.95 = 0.05 => offset so that all the particles are within the box!

	for(unsigned i = 0; i*latticeConstant < xLength; i++){
		currentZeroPos[0] = currentZeroPos[0] +latticeConstant;
		currentZeroPos[1] = off1-shielding + 1.05*latticeConstant;
		currentZeroPos[2] = -0.95*latticeConstant;
		for(unsigned j = 0; j*latticeConstant < yLength; j++){
			currentZeroPos[1] = currentZeroPos[1] - latticeConstant; // minus because the start is at the top of the wall going to the bottom of the box afterwards
			currentZeroPos[2] =  -0.95*latticeConstant; //_moleculePosition[2][1] - latticeConstant;
			for(unsigned k = 0; k*latticeConstant < zLength; k++){
				currentZeroPos[2] = currentZeroPos[2] +latticeConstant;
				_moleculePosition[0][_numberOfMolecules] = currentZeroPos[0];
				_moleculePosition[1][_numberOfMolecules] = currentZeroPos[1];
				_moleculePosition[2][_numberOfMolecules] = currentZeroPos[2];
				_moleculeCID[_numberOfMolecules] = ((int)floor( _moleculePosition[0][_numberOfMolecules]/ stripeWidth) % 2) +2; // if modulo returns zero => cid == 2, otherwise cid == 3
				_numberOfMolecules++;
				for(unsigned l = 0; l < 3; l++){
					_moleculePosition[0][_numberOfMolecules] = currentZeroPos[0] + pos[l][0];
					_moleculePosition[1][_numberOfMolecules] = currentZeroPos[1] - pos[l][1];
					_moleculePosition[2][_numberOfMolecules] = currentZeroPos[2] + pos[l][2];
					_moleculeCID[_numberOfMolecules] = ((int)floor( _moleculePosition[0][_numberOfMolecules]/ stripeWidth) % 2) +2;
					_numberOfMolecules++;
				} // end for(l)
			}// end for (k)
		} // end for(j)
	} // end for(i)

} // end of method

/*
void Molecule::calculateCoordinatesOfWallMolecule(double xLength, double zLength, double off1,
                                                 double latticeConstant, double shielding, unsigned numberOfCircles, unsigned deltaRCircle){
       double yLength = off1 - shielding;
       double pos[3][3];       // 3 vectors which are added to the current position in order to set three new particles
       for(unsigned i = 0; i<3; i++){    // !!!this initialisation of the 3 vectors is valid only for equilateral lattices!!! => see sketch!
               for(unsigned j= 0; j<3; j++){
                       if(i!= j) pos[i][j] = 0.5*latticeConstant;
                       else pos[i][j] = 0;
                       //cout << "pos["<<i<<"]["<< j<<"] = " << pos[i][j] << "\n";
               }
       }
       _numberOfMolecules = 0;
       if (deltaRCircle < 0.5*latticeConstant){
               std::cerr << "Circle size chosen too small! delta_r has to be at least 1*lattice constant! delta_r SWITCHED to 1*lattice constant.";
               deltaRCircle = latticeConstant;
       }
       double currentZeroPos[3];
       currentZeroPos[0] = -0.95*latticeConstant; // 1-0.95 = 0.05 => offset so that all the particles are within the box!
       for(unsigned i = 0; i*latticeConstant < xLength; i++){
               currentZeroPos[0] = currentZeroPos[0] +latticeConstant;
               currentZeroPos[1] = off1-shielding + 1.05*latticeConstant;
               currentZeroPos[2] = -0.95*latticeConstant;
               for(unsigned j = 0; j*latticeConstant < yLength; j++){
                       currentZeroPos[1] = currentZeroPos[1] - latticeConstant; // minus because the start is at the top of the wall going to the bottom of the box afterwards
                       currentZeroPos[2] =  -0.95*latticeConstant; //_moleculePosition[2][1] - latticeConstant;
                       for(unsigned k = 0; k*latticeConstant < zLength; k++){
                               currentZeroPos[2] = currentZeroPos[2] +latticeConstant;
                               _moleculePosition[0][_numberOfMolecules] = currentZeroPos[0];
                               _moleculePosition[1][_numberOfMolecules] = currentZeroPos[1];
                               _moleculePosition[2][_numberOfMolecules] = currentZeroPos[2];
                               _moleculeCID[_numberOfMolecules] = ((int)floor( _moleculePosition[0][_numberOfMolecules]/ stripeWidth) % 2) +2; // if modulo returns zero => cid == 2, otherwise cid == 3
                               _numberOfMolecules++;
                               for(unsigned l = 0; l < 3; l++){
                                       _moleculePosition[0][_numberOfMolecules] = currentZeroPos[0] + pos[l][0];
                                       _moleculePosition[1][_numberOfMolecules] = currentZeroPos[1] - pos[l][1];
                                       _moleculePosition[2][_numberOfMolecules] = currentZeroPos[2] + pos[l][2];
                                       _moleculeCID[_numberOfMolecules] = ((int)floor( _moleculePosition[0][_numberOfMolecules]/ stripeWidth) % 2) +2;
                                       _numberOfMolecules++;
                               } // end for(l)
                       }// end for (k)
               } // end for(j)
       } // end for(i)

  
} // end of method
*/

