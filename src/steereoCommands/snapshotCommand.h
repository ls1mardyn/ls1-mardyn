#ifndef __SNAPSHOTCOMMAND_H__
#define __SNAPSHOTCOMMAND_H__

#ifdef STEEREO
#include "Simulation.h"
#include <steereo/steereoCommand.h>


class SnapshotCommand : public SteereoCommand
{
 public:
  SnapshotCommand ();
  ~SnapshotCommand ();
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
  bool sendVelocity;
  bool sendForces;
  bool sendV2;
  bool sendV2Max;

};

#endif

#endif
