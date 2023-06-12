/*
 * GlobalStartGeometry.cpp
 *
 *  Created on: 12.01.2012
 *      Author: Stefan Becker <stefan.becker@mv.uni-kl.de>
 */


#include"GlobalStartGeometry.h"



const unsigned PRECISION = 10; // the number of times the loops are run determining the fill array _fill[i][j][k][l]

extern double LATTICE_CONST_WALL_LJTS;
extern double RHO_CRITICAL_1CLJ;
//extern const double LATTICE_CONST_CU;

//@brief: zero temperature lattice constant of the LJ model used for copper, needed to build up the wall
//@todo: more robust way necessary: example of a problem: if sigma of the LJ wall is changed => Lattice Constant will change, too!
//@todo: CU_LC to be shifted to the class PhaseSpaceWriter
//const double CU_LC = 6.83325;

GlobalStartGeometry::GlobalStartGeometry(unsigned in_nFluid, double in_rhoLiq, double in_rhoVap, double in_alpha, double in_beta, double in_gamma)
{
	_nFluid = in_nFluid;
	_rhoLiq = in_rhoLiq;
	_rhoVap = in_rhoVap;
	_alpha = in_alpha;
	_beta = in_beta;
	_gamma = in_gamma;
	_nLiq = _nFluid / (1+ (_alpha*_beta*_gamma-1.0) *_rhoVap/_rhoLiq );
	_nVap = _nFluid - _nLiq;
	std::cout << "\n**********************************\n";
	std::cout << "GloablStartGeometry:\n";
	//cout << "rhoVap = "<< _rhoVap << "\t rhoLiq = " << _rhoLiq << "\n";
	std::cout << "N liquid: " << _nLiq << "\n";
	std::cout << "N vapor: " << _nVap << "\n";
}

GlobalStartGeometry::~GlobalStartGeometry(){

}

void GlobalStartGeometry::calculateBoxFluidOffset(double hWall, double shielding)
{
	double aApprox;
	unsigned nEkX;	// number of elementary lattice boxes of the solid wall in x-direction
	unsigned nEkZ;	// "		"			"			"		"		"	 "	z-direction
	double effVap[3];

	// @brief: 2nd step: determining the simulation box, offset and fluid cuboid dimensions
	aApprox = pow(_nLiq/_rhoLiq, 1.0/3.0);
	nEkX = round(_alpha*aApprox / LATTICE_CONST_WALL_LJTS); 	// the number of lattices in the x-direction (Ek == elementary cristal)
	nEkZ = ceil(_gamma*aApprox /  LATTICE_CONST_WALL_LJTS);	// the number of lattices in the z-direction (Ek == elementary cristal)
	//cout << "nEKx = "<< nEkX << " \t nEkz = "<< nEkZ << "\n";

	_box[0] = nEkX * LATTICE_CONST_WALL_LJTS;
	_box[2] = nEkZ * LATTICE_CONST_WALL_LJTS;

	// introducing the factor pow(2.0, 1.0/3.0) in order to allow the liquid cuboid to have different edge length:
	// the edge lengths in x- and z-direction are twice the edge length in y-direction
	_effLiq[0] = pow(1.0 , 1.0/3.0) * _box[0] / _alpha;
	_effLiq[2] = pow(1.0 , 1.0/3.0) * _box[2] / _gamma;
	_effLiq[1] = (_nLiq / _rhoLiq) / (_effLiq[0]*_effLiq[2]);

	_box[1] = _nVap/_rhoVap/_box[0]/_box[2] + (_effLiq[0]*_effLiq[1]*_effLiq[2])/_box[0]/_box[2] + hWall + shielding;
	std::cout << "hWall: " << hWall << "\n";
	std::cout << "shielding: " << shielding <<"\n";
	std::cout << "box[0]: " << _box[0] << "\n";
	std::cout << "box[1]: " << _box[1] << "\n";
	std::cout << "box[2]: " << _box[2] << "\n";

	effVap[0] = 0.97*_box[0];
	effVap[1] = 0.97*(_box[1]-hWall-shielding) ;
	effVap[2] = 0.97*_box[2];
	_grossFluidDens = _nFluid / (_box[0] * _box[2] * (_box[1] - hWall));
	std::cout << "gross fluid density: " <<  _grossFluidDens << "\n";
	std::cout << "critical density 1CLJ: " << RHO_CRITICAL_1CLJ << "\n";
	std::cout << "vapour fraction x = " << (1.0/_grossFluidDens - 1.0/_rhoLiq) / (1.0/_rhoVap - 1.0/_rhoLiq) << "\n";

	_offLiq[0] = 0.5*(_box[0] - _effLiq[0]);
	_offLiq[1] = hWall + shielding;
	_offLiq[2] = 0.5*(_box[2]-_effLiq[2]);

	// @brief: 3rd step: setting up a cuboid containing all the fluid particles at liquid phase density
	double nLiqBoxes;
	nLiqBoxes = (_nLiq )/3.0;
	_liqUnits[0] = round( pow( _effLiq[0]*_effLiq[0]*nLiqBoxes/_effLiq[1]/_effLiq[2], 1.0/3.0 ) );
	//double fluidUnitsYz = (double)nFluidBoxes/_fluidUnits[0];
	_liqUnits[1] = round( sqrt(_effLiq[1]*nLiqBoxes/_liqUnits[0]/ _effLiq[2]) );
	_liqUnits[2] = ceil( nLiqBoxes/_liqUnits[0]/_liqUnits[1] );

	_liqUnit[0] =1.1* (double)_effLiq[0]/_liqUnits[0];
	_liqUnit[1] =1.1* (double)_effLiq[1]/_liqUnits[1];
	_liqUnit[2] =1.1* (double)_effLiq[2]/_liqUnits[2];


	_liqFillProbability = (_nLiq)/3.0 / ( _liqUnits[0] * _liqUnits[1] * _liqUnits[2] );

	// calculating the number of elementary vapour lattice boxes per direction
	// the number calculated assumes the entire fluid volume being available for the vapor
	// Since the liquid cuboid demands a part of this volume the number of lattice boxes is higher
	// making sure that enough lattice boxes are availabe for the vapour particles.
	// In fact 4 times more lattice slots than actual vapor particles are available.
	// In order to preserve the desired vapour density (number of vapour particles) the filling of the vapour lattice is controlled
	// by the bool-array _fillVap

	double nVapBoxes;
	nVapBoxes = 4.0*_nVap/3.0;
	_vapUnits[0] = floor(pow(effVap[0]*effVap[0]*nVapBoxes/effVap[1]/effVap[2],1.0/3.0) + 0.5);
	_vapUnits[1] = floor(sqrt(effVap[1]*nVapBoxes/_vapUnits[0]/effVap[2])  + 0.5);
	_vapUnits[2] = ceil(nVapBoxes/_vapUnits[0]/_vapUnits[1]);

	_vapUnit[0] = (double)effVap[0]/_vapUnits[0];
	_vapUnit[1] = (double)effVap[1]/_vapUnits[1];
	_vapUnit[2] = (double)effVap[2]/_vapUnits[2];

	std::cout << "_vapUnit[0]  = "<< _vapUnit[0] << "\t _vapUnit[1] = " << _vapUnit[1] << "\t _vapUnit[2] = " << _vapUnit[2] << "\n";
	std::cout << "_vapUnits[0]  = "<< _vapUnits[0] << "\t _vapUnits[1] = " << _vapUnits[1] << "\t _vapUnits[2] = " << _vapUnits[2] << "\n";
	std::cout << "effVap[0] = " << effVap[0] << "\t effVap[1] = " << effVap[1] << "\teffVap[2] = " << effVap[2] << "\n";

	std::cout << "_effLiq[0]  = " << _effLiq[0] << "\t _effLiq[1] = " << _effLiq[1] << "\t _effLiq[2] = " << _effLiq[2] <<"\n";
	std::cout << "upper edge of the liquid cuboid: _effLiq[1] + _offLiq[1] = " << _effLiq[1] + _offLiq[1] << "blub\n";


	_vapFillProbability = _nVap/3.0 / ( _vapUnits[0] * _vapUnits[1] * _vapUnits[2] );
	std::cout << "_vapFillProbability \t "<< _vapFillProbability << "\n";

	_offVap[0] = 0.1 * _vapUnit[0];
	_offVap[1] = _offLiq[1];
	_offVap[2] = 0.1 * _vapUnit[2];


}

void GlobalStartGeometry::calculateBoxFluidOffset(double hWall, double shielding, double edgeProportion)
{
	double aApprox;
	unsigned nEkX;	// number of elementary lattice boxes of the solid wall in x-direction
	unsigned nEkZ;	// "		"			"			"		"		"	 "	z-direction
	double effVap[3];

	// @brief: 2nd step: determining the simulation box, offset and fluid cuboid dimensions
	aApprox = pow(_nLiq/_rhoLiq, 1.0/3.0);
	nEkX = round(_alpha*aApprox / LATTICE_CONST_WALL_LJTS); 	// the number of lattices in the x-direction (Ek == elementary cristal)
	nEkZ = ceil(_gamma*aApprox /  LATTICE_CONST_WALL_LJTS);	// the number of lattices in the z-direction (Ek == elementary cristal)
	//cout << "nEKx = "<< nEkX << " \t nEkz = "<< nEkZ << "\n";

	_box[0] = nEkX * LATTICE_CONST_WALL_LJTS;
	_box[2] = nEkZ * LATTICE_CONST_WALL_LJTS;

	// introducing the factor pow(2.0, 1.0/3.0) in order to allow the liquid cuboid to have different edge length:
	// the edge lengths in x- and z-direction are twice the edge length in y-direction
	_effLiq[0] = pow(edgeProportion , 1.0/3.0) * _box[0] / _alpha;
	_effLiq[2] = pow(edgeProportion , 1.0/3.0) * _box[2] / _gamma;
	_effLiq[1] = (_nLiq / _rhoLiq) / (_effLiq[0]*_effLiq[2]);

	_box[1] = _nVap/_rhoVap/_box[0]/_box[2] + (_effLiq[0]*_effLiq[1]*_effLiq[2])/_box[0]/_box[2] + hWall + shielding;
	std::cout << "hWall: " << hWall << "\n";
	std::cout << "shielding: " << shielding <<"\n";
	std::cout << "box[0]: " << _box[0] << "\n";
	std::cout << "box[1]: " << _box[1] << "\n";
	std::cout << "box[2]: " << _box[2] << "\n";

	effVap[0] = 0.97*_box[0];
	effVap[1] = 0.97*(_box[1]-hWall-shielding) ;
	effVap[2] = 0.97*_box[2];
	_grossFluidDens = _nFluid / (_box[0] * _box[2] * (_box[1] - hWall));
	std::cout << "gross fluid density: " <<  _grossFluidDens << "\n";
	std::cout << "critical density 1CLJ: " << RHO_CRITICAL_1CLJ << "\n";
	std::cout << "vapour fraction x = " << (1.0/_grossFluidDens - 1.0/_rhoLiq) / (1.0/_rhoVap - 1.0/_rhoLiq) << "\n";

	_offLiq[0] = 0.5*(_box[0] - _effLiq[0]);
	_offLiq[1] = hWall + shielding;
	_offLiq[2] = 0.5*(_box[2]-_effLiq[2]);

	// @brief: 3rd step: setting up a cuboid containing all the fluid particles at liquid phase density
	double nLiqBoxes;
	nLiqBoxes = (_nLiq )/3.0;
	_liqUnits[0] = round( pow( _effLiq[0]*_effLiq[0]*nLiqBoxes/_effLiq[1]/_effLiq[2], 1.0/3.0 ) );
	//double fluidUnitsYz = (double)nFluidBoxes/_fluidUnits[0];
	_liqUnits[1] = round( sqrt(_effLiq[1]*nLiqBoxes/_liqUnits[0]/ _effLiq[2]) );
	_liqUnits[2] = ceil( nLiqBoxes/_liqUnits[0]/_liqUnits[1] );

	_liqUnit[0] = (double)_effLiq[0]/_liqUnits[0];
	_liqUnit[1] = (double)_effLiq[1]/_liqUnits[1];
	_liqUnit[2] = (double)_effLiq[2]/_liqUnits[2];


	_liqFillProbability = (_nLiq)/3.0 / ( _liqUnits[0] * _liqUnits[1] * _liqUnits[2] );

	// calculating the number of elementary vapour lattice boxes per direction
	// the number calculated assumes the entire fluid volume being available for the vapor
	// Since the liquid cuboid demands a part of this volume the number of lattice boxes is higher
	// making sure that enough lattice boxes are availabe for the vapour particles.
	// In fact 4 times more lattice slots than actual vapor particles are available.
	// In order to preserve the desired vapour density (number of vapour particles) the filling of the vapour lattice is controlled
	// by the bool-array _fillVap

	double nVapBoxes;
	nVapBoxes = 4.0*_nVap/3.0;
	_vapUnits[0] = floor(pow(effVap[0]*effVap[0]*nVapBoxes/effVap[1]/effVap[2],1.0/3.0) + 0.5);
	_vapUnits[1] = floor(sqrt(effVap[1]*nVapBoxes/_vapUnits[0]/effVap[2])  + 0.5);
	_vapUnits[2] = ceil(nVapBoxes/_vapUnits[0]/_vapUnits[1]);

	_vapUnit[0] = (double)effVap[0]/_vapUnits[0];
	_vapUnit[1] = (double)effVap[1]/_vapUnits[1];
	_vapUnit[2] = (double)effVap[2]/_vapUnits[2];

	std::cout << "_vapUnit[0]  = "<< _vapUnit[0] << "\t _vapUnit[1] = " << _vapUnit[1] << "\t _vapUnit[2] = " << _vapUnit[2] << "\n";
	std::cout << "_vapUnits[0]  = "<< _vapUnits[0] << "\t _vapUnits[1] = " << _vapUnits[1] << "\t _vapUnits[2] = " << _vapUnits[2] << "\n";
	std::cout << "effVap[0] = " << effVap[0] << "\t effVap[1] = " << effVap[1] << "\teffVap[2] = " << effVap[2] << "\n";

	std::cout << "_effLiq[0]  = " << _effLiq[0] << "\t _effLiq[1] = " << _effLiq[1] << "\t _effLiq[2] = " << _effLiq[2] <<"\n";
	std::cout << "upper edge of the liquid cuboid: _effLiq[1] + _offLiq[1] = " << _effLiq[1] + _offLiq[1] << "\n";


	_vapFillProbability = _nVap/3.0 / ( _vapUnits[0] * _vapUnits[1] * _vapUnits[2] );

	_offVap[0] = 0.1 * _vapUnit[0];
	_offVap[1] = _offLiq[1];
	_offVap[2] = 0.1 * _vapUnit[2];


}

//@todo: the tensor accounting for the filling of a single slot, i.e. _fill[][][][], should be part of an extra class that also determines the
// posistion of the particles(both fluid and wall)
void GlobalStartGeometry::calculateLiqFillProbabilityArray()
{
	//@brief: initializing each element of _fill[][][][] by true
	_fill[_liqUnits[0]][_liqUnits[1]][_liqUnits[2]][3] = new bool;
	for(unsigned i = 0; i< _liqUnits[0]; i++){
		for(unsigned j = 0; j< _liqUnits[1]; j++){
			for(unsigned k = 0; k< _liqUnits[2]; k++){
				for(unsigned l = 0; l < 3; l++){
					_fill[i][j][k][l] = true;
				}
			}
		}
	}

	// number of filled slots, initialized by the total number of slots (both filled and not filled)
	//unsigned nFilledSlots;
	// total number of slots
	double totalNSlots;
	// a measure of the probability by which an element of _fill is swapped
	double pSwap;
	// number of slots beeing ideally filled
	double nIdeallyFilled;
	// bool valued quantity checking if particles have to be added or removed respectively
	bool tSwap;

	_nFilledLiqSlots = 3*_liqUnits[0]*_liqUnits[1]*_liqUnits[2];
//	std::cout << "number of filled slots at the beginning of Gloablstartgeometry:" << _nFilledSlots <<"\n";
	totalNSlots = _nFilledLiqSlots; // slots is and "always" will be the total number of slots
	nIdeallyFilled = _liqFillProbability * totalNSlots;
	RandomNumber rdm;

	for(unsigned m = 0; m < PRECISION; m++){
		tSwap = (_nFilledLiqSlots < nIdeallyFilled);
		pSwap = (nIdeallyFilled - (double)_nFilledLiqSlots)/( (tSwap ? totalNSlots : 0.0) - _nFilledLiqSlots );
		for(unsigned i=0; i < _liqUnits[0]; i++){
			for(unsigned j=0; j < _liqUnits[1]; j++){
				for(unsigned k=0; k < _liqUnits[2]; k++){
					for(unsigned short l=0; l < 3; l++){
						if(pSwap >= rdm.randNum()){
							if(_fill[i][j][k][l]) _nFilledLiqSlots--;
							_fill[i][j][k][l] = tSwap;
							if(tSwap) _nFilledLiqSlots++;
						}
					}
				}
			}
		}
	}
/*	std::cout << "Filling" << _nFilledSlots << " out of a total number of " << totalNSlots << "fluid slots. Ideally a number of "
		<< nIdeallyFilled << "fluid slots was to be filled.\n" << "fluidFillProbability = " << _fluidFillProbability << "\n";*/
}

void GlobalStartGeometry::calculateVapFillProbabilityArray(){

  //@brief: initializing each element of _fillVap[][][][] by true
	for(unsigned i = 0; i< _vapUnits[0]; i++){
		for(unsigned j = 0; j< _vapUnits[1]; j++){
			for(unsigned k = 0; k< _vapUnits[2]; k++){
				for(short l = 0; l < 3; l++){
					_fillVap[i][j][k][l] = true;
				}
			}
		}
	}


	// number of filled slots, initialized by the total number of slots (both filled and not filled)
	// total number of slots
	double totalNSlots;
	// a measure of the probability by which an element of _fillVap is swapped
	double pSwap;
	// number of slots beeing ideally filled
	double nIdeallyFilled;
	// bool valued quantity checking if particles have to be added or removed respectively
	bool tSwap;

	_nFilledVapSlots = 3*_vapUnits[0]*_vapUnits[1]*_vapUnits[2];
//	std::cout << "number of filled slots at the beginning of Gloablstartgeometry:" << _nFilledSlots <<"\n";
	totalNSlots = _nFilledVapSlots;
	nIdeallyFilled = _vapFillProbability * totalNSlots;
	RandomNumber rdm;

	for(unsigned m = 0; m < PRECISION; m++){
		tSwap = (_nFilledVapSlots < nIdeallyFilled);
		pSwap = (nIdeallyFilled - (double)_nFilledVapSlots)/( (tSwap ? totalNSlots : 0.0) - _nFilledVapSlots );
		for(long int i=0; i < _vapUnits[0]; i++){
			for(long int j=0; j < _vapUnits[1]; j++){
				for(long int k=0; k < _vapUnits[2]; k++){
				  // if(...) if the position of the slots is overlapping with the liquid cuboid => no slot beeing filled!!!
				  if ((i+2)* _vapUnit[0] >= 0.5*(_box[0]-_effLiq[0]) && (i-1) * _vapUnit[0] <= 0.5 * (_box[0]+_effLiq[0])
				    && (j-1)* _vapUnit[1] <= _effLiq[1]
				    && (k+2)* _vapUnit[2] >= 0.5*(_box[2]-_effLiq[2]) && (k-1)* _vapUnit[2] <= 0.5 * (_box[2] +_effLiq[2]) )
				  {
				    for(unsigned short l = 0; l < 3; l++){
				      if(_fillVap[i][j][k][l]) _nFilledVapSlots--;
				       _fillVap[i][j][k][l] = false;
				    }// end for
				  } // end if
				  else{
				    for(unsigned short l=0; l < 3; l++){
						  if(pSwap >= rdm.randNum()){
							  if(_fillVap[i][j][k][l]) _nFilledVapSlots--;
							  _fillVap[i][j][k][l] = tSwap;
							  if(tSwap) _nFilledVapSlots++;
						}
					}
				  } // end else
				}
			}
		}
	}
/*	std::cout << "Filling" << _nFilledSlots << " out of a total number of " << totalNSlots << "fluid slots. Ideally a number of "
		<< nIdeallyFilled << "fluid slots was to be filled.\n" << "fluidFillProbability = " << _fluidFillProbability << "\n";*/


}

bool GlobalStartGeometry::gFillLiqArray(unsigned fluidUnits0, unsigned fluidUnits1, unsigned fluidUnits2, unsigned particleInElementaryBox){
	return _fill[fluidUnits0][fluidUnits1][fluidUnits2][particleInElementaryBox];
}

bool GlobalStartGeometry::gFillVapArray(unsigned vapUnits0, unsigned vapUnits1, unsigned vapUnits2, unsigned particleInElementaryBox){
  return _fillVap[vapUnits0][vapUnits1][vapUnits2][particleInElementaryBox];
}

double GlobalStartGeometry::gBoxLength(unsigned direction){
	return _box[direction];
}

unsigned GlobalStartGeometry::gLiqUnits(unsigned direction){
	return _liqUnits[direction];
}

unsigned GlobalStartGeometry::gVapUnits(unsigned short direction){
	return _vapUnits[direction];
}

double GlobalStartGeometry::gOffsetLiq(unsigned direction){
	return _offLiq[direction];
}

double GlobalStartGeometry::gOffsetVap(unsigned short direction){
	return _offVap[direction];
}

double GlobalStartGeometry::gLiqUnit(unsigned short direction){
	return _liqUnit[direction];
}

double GlobalStartGeometry::gVapUnit(unsigned short direction){
  return _vapUnit[direction];
}

unsigned GlobalStartGeometry::gNFilledLiqSlots(){
	return _nFilledLiqSlots;
}

unsigned GlobalStartGeometry::gNFilledVapSlots(){
	return _nFilledVapSlots;
}

double GlobalStartGeometry::gGrossFluidDens(){
  return _grossFluidDens;
}
