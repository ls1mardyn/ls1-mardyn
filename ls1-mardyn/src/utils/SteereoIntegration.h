#ifdef STEEREO
/*
 * SteereoIntegration.h
 *
 *  Created on: Jun 29, 2010
 *      Author: hpcdjenz
 */

#ifndef STEEREOINTEGRATION_H_
#define STEEREOINTEGRATION_H_

class SteereoSimSteering;
class ParticleContainer;
class SteereoCouplingSim;

SteereoSimSteering* initSteereo(int ownRank, int numProcs);
#ifdef STEEREO_COUPLING
SteereoCouplingSim* initCoupling (SteereoSimSteering* sim, long* simStep);
#endif
void registerSteereoCommands (SteereoSimSteering* simSteer, Simulation* sim);
void startListeningSteereo (SteereoSimSteering* simSteer);
void checkMoleculeContainer (ParticleContainer* pc);


#endif /* STEEREOINTEGRATION_H_ */
#endif /* STEEREO */
