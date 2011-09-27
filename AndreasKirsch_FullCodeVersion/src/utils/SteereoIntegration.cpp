/*
 * SteereoIntegration.cpp
 *
 *  Created on: Jun 29, 2010
 *      Author: hpcdjenz
 */

#ifdef STEEREO

#include "utils/Logger.h"
#include "commands/snapshotCommand.h"
#include "commands/megaMolSnapshotCommand.h"
#include "commands/sendCouplingMDCommand.h"
#include "commands/receiveCouplingMDCommand.h"
#include "Domain.h"
#ifdef PARALLEL
#include <mpi.h>
#include <steereoMPIIntraCommunicator.h>
#endif //PARALLEL
#include <steerParameterCommand.h>
#include <steereoSocketCommunicator.h>
#include <simSteering.h>

using Log::global_log;

SimSteering* initSteereo(int startPort, int ownRank) {

	SteereoLogger::setOutputLevel(2);
	SimSteering* _steer = new SimSteering();
	_steer->setNumberOfQueues(1);
	int portNumber = 44445;
#ifdef PARALLEL
	int ownSize;
	MPI_Comm_size(MPI_COMM_WORLD, &ownSize);
	SteereoMPIIntraCommunicator* mpiIntraComm = new SteereoMPIIntraCommunicator();
	int partNum = 1;
#ifdef STEEREO_PARTITIONS
	partNum = STEEREO_PARTITIONS;
#endif /* STEEREO_PARTITIONS */
	char* partVal = getenv("STEEREO_PARTITIONS");
	if (partVal != NULL) {
		partNum = atoi(partVal);
	}
	global_log->debug() << "going to divide the " << ownSize << " processes into " << partNum << " partitions" << std::endl;
	mpiIntraComm->generateEqually(ownRank, partNum, ownSize);
	global_log->debug() << "equally generated" << std::endl;
	_steer->setIntraCommunicator(mpiIntraComm);
	if (mpiIntraComm->amIRoot()) {
		portNumber += (ownRank * partNum) / ownSize;
#endif /* PARALLEL */
	std::stringstream strstr;
	strstr << portNumber;
	_steer->setCommunicator(new SteereoSocketCommunicator(strstr.str()));
#ifdef PARALLEL
}
#endif /* PARALLEL */
	return _steer;
}

void registerSteereoCommands(SimSteering* simSteer, Simulation* sim) {

	simSteer->registerCommand(MegaMolSnapshotCommand::generateNewInstance, "getMegaMolSnapshot");
	simSteer->registerCommand(SendCouplingMDCommand::generateNewInstance, "sendCouplingMD");
	simSteer->registerCommand(ReceiveCouplingMDCommand::generateNewInstance, "receiveCouplingMD");
	simSteer->registerCommand(SnapshotCommand::generateNewInstance, "getSnapshot");
	MegaMolSnapshotCommand::setSimData(sim);
	SnapshotCommand::setSimData(sim);
	SteerParameterCommand::registerScalarParameter("temp", sim->getDomain(), &Domain::getGlobalCurrentTemperature,
			&Domain::setGlobalTemperature);
	SendCouplingMDCommand::addData("sendCouplingMD", sim);
	ReceiveCouplingMDCommand::addData("receiveCouplingMD", sim);
}

void startListeningSteereo (SimSteering* simSteer, int ownrank)
{
#ifdef PARALLEL
	if (simSteer->getIntraCommunicator()->amIRoot()) {
#endif
	simSteer->startListening();
#ifdef PARALLEL
}
#endif
}
#endif //STEEREO
