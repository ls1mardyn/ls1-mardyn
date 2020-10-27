/*
 * SteereoIntegration.cpp
 *
 *  Created on: Jun 29, 2010
 *      Author: hpcdjenz
 */

#ifdef STEEREO

#include "utils/Logger.h"
#include "Domain.h"
#include "steereoCommands/snapshotCommand.h"
#include "steereoCommands/megaMolSnapshotCommand.h"
#include "steereoCommands/sendCouplingMDCommand.h"
#include "steereoCommands/receiveCouplingMDCommand.h"
#include "steereoCommands/estimateRemainingTimeCommand.h"
//#include "steereoCommands/getVisDataCommand.h"
#include <steereo/steerParameterCommand.h>
#include <steereo/steereoSocketCommunicator.h>
#include <steereo/steereoSimSteering.h>
#include <steereo/steereoCouplingSim.h>
#include <steereo/steereoXMLReader.h>
#include "particleContainer/ParticleContainer.h"

#ifdef ENABLE_MPI
#include "mpi.h"
#include <steereoMPIIntraCommunicator.h>
#endif //ENABLE_MPI

using Log::global_log;

SteereoSimSteering* initSteereo(int ownRank, int numProcs) {

	SteereoLogger::setOutputLevel(2);
	global_log->info() << "In initSteereo" << std::endl;
	SteereoSimSteering* _steer = new SteereoSimSteering();
/*	SteereoXMLReader* xmlReader;
#ifdef ENABLE_MPI
	global_log->debug() << "I am parallel" << std::endl;
	xmlReader = new SteereoXMLReader ("./utils/steereoParallelConfig.xml");
#else
	global_log->debug() << "I am sequential" << std::endl;
	xmlReader = new SteereoXMLReader ("./utils/steereoSequentialConfig.xml");
#endif
	xmlReader->configSimulation(_steer, ownRank, numProcs);
	global_log->info() << "configured simulation" << std::endl;*/
#ifdef STEEREO_COUPLING
	_steer->setNumberOfQueues(2);
#else
	_steer->setNumberOfQueues(1);
#endif
	int portNumber = 44445;
#ifdef ENABLE_MPI
	int ownSize;
	MPI_Comm_size(MPI_COMM_WORLD, &ownSize);
	SteereoMPIIntraCommunicator* mpiIntraComm = new SteereoMPIIntraCommunicator();
	int partNum = 1;
#ifdef STEEREO_PARTITIONS
	partNum = STEEREO_PARTITIONS;
#endif // STEEREO_PARTITIONS
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
#endif // ENABLE_MPI

	std::stringstream strstr;
	strstr << portNumber;
	_steer->setCommunicator(new SteereoSocketCommunicator(strstr.str()));
#ifdef ENABLE_MPI
}
#endif // ENABLE_MPI
	global_log->info() << "done init_steereo" << std::endl;
	return _steer;

}

#ifdef STEEREO_COUPLING
SteereoCouplingSim* initCoupling(SteereoSimSteering* simSteer, long* stepNum) {
	SteereoCouplingSim* coupling = new SteereoCouplingSim (simSteer);
	coupling->setQueueID(0);
	coupling->registerReceiveCommand(ReceiveCouplingMDCommand::generateNewInstance, 1, 1);
	coupling->registerSendCommand(SendCouplingMDCommand::generateNewInstance, 1, 1);

	coupling->registerStepCounter(stepNum);
	return coupling;
}
#endif

void registerSteereoCommands(SteereoSimSteering* simSteer, Simulation* sim) {

	simSteer->registerCommand(MegaMolSnapshotCommand::generateNewInstance, "getMegaMolSnapshot");
	simSteer->registerCommand(SnapshotCommand::generateNewInstance, "getSnapshot");
	//simSteer->registerCommand(GetVisDataCommand::generateNewInstance, "getVisData");
	// register estimateRemainingTimeCommand
    simSteer->registerCommand(EstimateRemainingTimeCommand::generateNewInstance, "estimateRemainingTime");
	MegaMolSnapshotCommand::setSimData(sim);
	//simSteer->registerSignalHandler(EstimateRemainingTimeCommand::generateNewInstance, 15);
	SnapshotCommand::setSimData(sim);

	SteerParameterCommand::registerScalarParameter(
	    "temp",
	    sim->getDomain(),
	    &Domain::getGlobalCurrentTemperature,
	    &Domain::setGlobalTemperature);
	//GetVisDataCommand::addData("getVisData", sim);
	global_log->info() << "add Data to EstimateRemainingTimeCommand" << std::endl;
	EstimateRemainingTimeCommand::addData ("estimateRemainingTime", sim);
	// add data for the EstimateRemainingTimeCommand

#ifdef STEEREO_COUPLING
	SendCouplingMDCommand::addData("sendCouplingMD", sim);
	ReceiveCouplingMDCommand::addData("receiveCouplingMD", sim);
#endif
}

void startListeningSteereo(SteereoSimSteering* simSteer) {
#ifdef ENABLE_MPI
	if (simSteer->getIntraCommunicator()->amIRoot()) {
#endif
	simSteer->startListening();
#ifdef ENABLE_MPI
}
#endif
}

void checkMoleculeContainer(ParticleContainer* pc) {
	std::cout << "In checkMoleculeContainer" << std::endl;

	double rmin; // lower corner of the process-specific domain //ENABLE_MPI
	double rmax;
	rmin = pc->getBoundingBoxMin(0);
	rmax = pc->getBoundingBoxMax(0);

	Molecule* currentMolecule;

	double low_limit; // particles below this limit have to be copied or moved to the lower process
	double high_limit; // particles above(or equal) this limit have to be copied or moved to the higher process

	// set limits (outside "inner" region)
	low_limit = rmin;
	high_limit = rmax;
	currentMolecule = pc->begin();
	int outOfBounds = 0;
	std::cout << "low_limit: " << low_limit << " / high_limit: " << high_limit << std::endl;
	while (currentMolecule != pc->end()) {
		for (int i = 0; i < 3; i++)
		{
			if ((currentMolecule->r(i) > high_limit) || (currentMolecule->r(i) < low_limit))
			{
				outOfBounds++;
				std::cout << "current molecule " << currentMolecule->getID() << " is out of bounds: "
						<< currentMolecule->r(0) << "," << currentMolecule->r(1) << "," << currentMolecule->r(2) << std::endl;
			}
		}
		currentMolecule = pc->next();
	}
	std::cout << outOfBounds << " molecules out of bounds " << std::endl;
	std::cout << "-------------------- end of checkMolecules " << std::endl;
}
#endif //STEEREO
