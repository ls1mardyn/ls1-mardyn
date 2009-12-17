/*
 * megaMolSnapshotCommand.h
 *
 *  Created on: Mar 19, 2009
 *      Author: hpcdjenz
 */

#ifdef STEEREO
#ifndef MEGAMOLSNAPSHOTCOMMAND_H_
#define MEGAMOLSNAPSHOTCOMMAND_H_

#include <steereoCommand.h>
#include <Simulation.h>

class MegaMolSnapshotCommand : public SteereoCommand
{
	public:
  MegaMolSnapshotCommand ();
  ~MegaMolSnapshotCommand ();
  virtual ReturnType execute ();
  void setParameters (std::list<std::string> params);
  static void setSimData (Simulation* simu) {sim = simu;};
  static SteereoCommand* generateNewInstance ();

  bool condition ();
  void setStepInterval (int interval) {stepInterval = interval;};


 private:
  // parameters needed for execution
  int sockfd;
  static Simulation* sim;
  static int startStep;
  int stepInterval;

};

#endif /* MEGAMOLSNAPSHOTCOMMAND_H_ */
#endif /* STEEREO */
