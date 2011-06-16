/*
 * ReceiveCouplingInfoCommand.h
 *
 *  Created on: Apr 24, 2009
 *      Author: hpcdjenz
 */

#ifdef STEEREO
#ifndef RECEIVECOUPLINGMDCOMMAND_H_
#define RECEIVECOUPLINGMDCOMMAND_H_

#include <steereo/steereoCouplingCommand.h>
#include <queue>
#include "../molecules/Molecule.h"

class CouplingBoundary;

typedef struct
{
	double translatedPos[3];
	//CouplingBoundary boundary;
} DataPackage;

class ReceiveCouplingMDCommand: public SteereoCouplingCommand
{
public:
	ReceiveCouplingMDCommand();
	virtual ~ReceiveCouplingMDCommand();

  virtual ReturnType executeProcessing();
  virtual ReturnType executeTransfer ();
  void setParameters (std::list<std::string> params);

  static SteereoCommand* generateNewInstance ();
  static void setReferenceMolecule (Molecule* mol) {referenceMolecule = new Molecule(*mol);};
  bool condition ();

private:
  double getPotentialEnergyAndForce (Molecule* victim);
  void updatePosition (Molecule* victim, double stepSize, double* forceModulus);
  double getNewStepSize (double currentSize, double potentialEnergy, double targetEnergy, double force);
  int insertParticle(double startPos[3], double maxDist, double targetEnergy, CouplingBoundary* boundary);

  std::vector<double*> receiveBuffers;
  std::queue<DataPackage> workPackages;
  double stepOverlap;
  double potentialEnergyOverlap;
  double minimalStepSize;
  double accuracy;
  double maximumDistance;
  double targetEnergy;
  int maximumIterations;
  int maximumTries;
  //int borderToLook;
  int numberOfWorkpackagesPerCycle;

	static Molecule* referenceMolecule;

};

#endif /* RECEIVECOUPLINGMDCOMMAND_H_ */
#endif /* STEEREO */

