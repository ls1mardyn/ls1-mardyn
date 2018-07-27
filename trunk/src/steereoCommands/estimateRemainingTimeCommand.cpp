/*
 * estimateRemainingTimeCommand.cpp
 *
 *  Created on: Jul 1, 2010
 *      Author: hpcdjenz
 */
#ifdef STEEREO

#include "estimateRemainingTimeCommand.h"
#include <Simulation.h>

EstimateRemainingTimeCommand::EstimateRemainingTimeCommand() : SteereoCommand (false, "estimateRemainingTime", 0) {}

EstimateRemainingTimeCommand::~EstimateRemainingTimeCommand() {}

ReturnType EstimateRemainingTimeCommand::execute () {
	Simulation* sim = (Simulation*) (this->getData(0));
	sim->stopTimer("SIMULATION_LOOP");
	double elapsedTime = sim->getTime("SIMULATION_LOOP");
	unsigned long currentStep = sim->getSimulationStep();
	unsigned long numberOfTimesteps = sim->getNumTimesteps();
	double estimation = ((double) numberOfTimesteps / currentStep) * elapsedTime - elapsedTime;
	SteereoStream stream;
	stream << currentStep << numberOfTimesteps << elapsedTime << estimation;
	logger->debug() << "estimating " << estimation << std::endl;
	if (this->getIntraCommunicator()->amIRoot ()) {
		this->getCommunicator()->sendStream (this->getConnection()->getDataConnection(), stream);
	}
	sim->startTimer("SIMULATION_LOOP");
	return EXECUTED;
}

void EstimateRemainingTimeCommand::setParameters (std::list<std::string> params) {}

SteereoCommand* EstimateRemainingTimeCommand::generateNewInstance () {
	return new EstimateRemainingTimeCommand;
}

bool EstimateRemainingTimeCommand::condition () {
	return true;
}

#endif /* STEEREO */
