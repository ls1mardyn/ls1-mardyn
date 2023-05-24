/*
 * GlobalStartGeometry.h
 *
 *  Created on: 12.01.2012
 *      Author:Stefan Becker <stefan.becker@mv.uni-kl.de>
 *
 * The instance of this class defines the global geometry of the start configuration, e.g. the dimensions of the simulation box
 * or the length of the liquid cuboid in relation to the box length, amongst others.
 */

#ifndef GLOBALSTARTGEOMETRY_H_
#define GLOBALSTARTGEOMETRY_H_

#include<iostream>
#include<string>
#include<math.h>
#include<map>
#include"RandomNumber.h"


class GlobalStartGeometry
{
public:
	GlobalStartGeometry(unsigned in_nFluid, double in_rhoLiq, double in_rhoVap, double in_alpha, double in_beta, double in_gamma);
	~GlobalStartGeometry();

	void calculateBoxFluidOffset(double hWall, double shielding);
	void calculateBoxFluidOffset(double hWall, double shielding, double edgeProportion);
	void calculateLiqFillProbabilityArray();
	void calculateVapFillProbabilityArray();

	double gBoxLength(unsigned direction);
	unsigned gLiqUnits(unsigned direction);
	unsigned gVapUnits(unsigned short direction);
	double gOffsetLiq(unsigned direction);
	double gOffsetVap(unsigned short direction);
	double gLiqUnit(unsigned short direction);
	double gVapUnit(unsigned short direction);
	bool gFillLiqArray(unsigned fluidUnits0, unsigned fluidUnits1, unsigned flunidUnits2, unsigned particleInElementaryBox);
	bool gFillVapArray(unsigned vapUnits0, unsigned vapUnits1, unsigned vapUnits2, unsigned particleInElementaryBox);
	unsigned gNFilledLiqSlots();
	unsigned gNFilledVapSlots();
	double gGrossFluidDens();

private:
	unsigned _nFluid;
	//unsigned _wallThick; // thickness of the wall as a multiple of the lattice constant in y-direction
	double _alpha;
	double _beta;
	double _gamma;

	//@brief: the next three variables needed to compute the density of the fluid (so far: 1CLJ)
	double _rhoLiq;
	double _rhoVap;
	double _grossFluidDens;
	double _sigmaFluid1CLJ;

	double _nLiq;
	double _nVap;
	unsigned _nFilledLiqSlots;
	unsigned _nFilledVapSlots;

	//@brief: lattice contant of the solid phase (wall) in the three dimensions
	//double _latticeConstSolid[3];
	//@brief: dimensions of the simulation box
	double _box[3];
	//@brief: offset of the liquid particles posistion with respect to the simulation box (x,z direction)
	// and with respect to the wall (y-direction)
	double _offLiq[3];
	//@brief: offset of the vapour particles posistion with respect to the simulation box (x,z direction)
	// and with respect to the wall (y-direction)
	double _offVap[3];
	// @brief: ectual length of the liquid cuboid in each direction
	double _effLiq[3];
	//@brief: number of elementary liquid lattices in three directions
	unsigned _liqUnits[3];
	//@brief: lengths of a single elementary liquid lattice in three directions
	double _liqUnit[3];
	//@brief: number of elementary vapour lattice boxes per direction
	long int _vapUnits[3];
	//@brief: lengths of a single elementary vapour lattice with respect to three directions
	double _vapUnit[3];
	//@brief: gross probability of a liquid elementary box to be filled
	double _liqFillProbability;
	//@brief: gross probability of a vapour elementary box to be filled
	double _vapFillProbability;
	//@brief: 4-dimensional array addressing single slots that may be filled with a particle
	// => 3 diemsions for a fluid elementary box (due to 3 directions in space)
	// and one dimension addressing one of three slots within an elementary box
	std::map< unsigned, std::map<unsigned, std::map<unsigned, std::map<unsigned, bool> > > > _fill;
	//@brief: 4-dimensional array addressing single slots that may be filled with a vapour particle
	// => 3 diemsions for a vapour elementary box (due to 3 directions in space)
	// and one dimension addressing one of three slots within an elementary box
	std::map< unsigned, std::map<unsigned, std::map<unsigned, std::map<unsigned, bool> > > > _fillVap;
	
};
#endif /* GLOBALSTARTGEOMETRY_H_ */
