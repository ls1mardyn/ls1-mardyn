/*
 * TemperatureControl.h
 *
 *  Created on: 27.05.2015
 *      Author: mheinen
 *
 *  Accumulate kinetic energy dependent on which translatoric directions should be thermostated
 */
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "molecules/Molecule.h"

class Accumulator {
private:
	bool _accumulateX, _accumulateY, _accumulateZ;

public:
	Accumulator(bool accX, bool accY, bool accZ) : _accumulateX(accX), _accumulateY(accY), _accumulateZ(accZ) {}

	double CalcKineticEnergyContribution(Molecule* mol) {
		double vx = _accumulateX ? mol->v(0) : 0.0;
		double vy = _accumulateY ? mol->v(1) : 0.0;
		double vz = _accumulateZ ? mol->v(2) : 0.0;
		double m = mol->mass();

		return m * (vx * vx + vy * vy + vz * vz);
	}
	void ScaleVelocityComponents(Molecule* mol, double vcorr) {
		if (_accumulateX) mol->setv(0, mol->v(0) * vcorr);
		if (_accumulateY) mol->setv(1, mol->v(1) * vcorr);
		if (_accumulateZ) mol->setv(2, mol->v(2) * vcorr);
	}
};