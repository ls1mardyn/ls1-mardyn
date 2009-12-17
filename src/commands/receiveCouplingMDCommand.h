/*
 * ReceiveCouplingInfoCommand.h
 *
 *  Created on: Apr 24, 2009
 *      Author: hpcdjenz
 */

#ifdef STEEREO
#ifndef RECEIVECOUPLINGMDCOMMAND_H_
#define RECEIVECOUPLINGMDCOMMAND_H_

#include <steereoCommand.h>

class ReceiveCouplingMDCommand: public SteereoCommand
{
public:
	ReceiveCouplingMDCommand();
	virtual ~ReceiveCouplingMDCommand();

  virtual ReturnType execute ();
  void setParameters (std::list<std::string> params);

  static SteereoCommand* generateNewInstance ();

  bool condition ();
  void setStepInterval (int interval) {stepInterval = interval;};

private:
  static int startStep;
  int stepInterval;
};

#endif /* RECEIVECOUPLINGMDCOMMAND_H_ */
#endif /* STEEREO */

