#ifdef STEEREO
/*
 * SteereoIntegration.h
 *
 *  Created on: Jun 29, 2010
 *      Author: hpcdjenz
 */

#ifndef STEEREOINTEGRATION_H_
#define STEEREOINTEGRATION_H_

class SimSteering;

SimSteering* initSteereo(int startPort, int ownRank);
void registerSteereoCommands (SimSteering* simSteer, Simulation* sim);
void startListeningSteereo (SimSteering* simSteer, int ownRank);


#endif /* STEEREOINTEGRATION_H_ */
#endif /* STEEREO */
