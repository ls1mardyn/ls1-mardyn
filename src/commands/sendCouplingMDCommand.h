/*
 * SendCouplingInfoCommand.h
 *
 *  Created on: Apr 23, 2009
 *      Author: hpcdjenz
 */

#ifdef STEEREO

#ifndef SENDCOUPLINGINFOCOMMAND_H_
#define SENDCOUPLINGINFOCOMMAND_H_

#include <steereoCommand.h>

class SendCouplingMDCommand: public SteereoCommand
{
public:
	SendCouplingMDCommand();
	virtual ~SendCouplingMDCommand();

  virtual ReturnType execute ();
  void setParameters (std::list<std::string> params);

  static SteereoCommand* generateNewInstance ();

  bool condition ();
  void setStepInterval (int interval) {stepInterval = interval;};

private:
  static int startStep;
  static int borderToLook;
  static int outmin, outmax;
  int stepInterval;
};

#endif /* SENDCOUPLINGINFOCOMMAND_H_ */
#endif /* STEEREO */
