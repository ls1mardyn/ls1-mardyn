/*
 * SendCouplingInfoCommand.h
 *
 *  Created on: Apr 23, 2009
 *      Author: hpcdjenz
 */

#ifdef STEEREO

#ifndef __SENDCOUPLINGINFOCOMMAND_H__
#define __SENDCOUPLINGINFOCOMMAND_H__

#include <steereo/steereoCouplingCommand.h>
#include <vector>
#include "../molecules/Molecule.h"

class CouplingBoundary;

class SendCouplingMDCommand: public SteereoCouplingCommand
{
public:
	SendCouplingMDCommand();
	virtual ~SendCouplingMDCommand();

  virtual ReturnType executeProcessing ();
  virtual ReturnType executeTransfer ();
  void setParameters (std::list<std::string> params);

  static SteereoCommand* generateNewInstance ();

  bool condition ();

private:
  void transferMolecules ();
  bool isInBounds (double ro1, double ro2, CouplingBoundary* boundary);

  int borderToLook;
  int outmin, outmax;

  static std::vector<Molecule>* transferContainer;
};

#endif /* SENDCOUPLINGINFOCOMMAND_H_ */
#endif /* STEEREO */
