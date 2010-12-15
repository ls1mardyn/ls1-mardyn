/*
 * estimateRemainingTimeCommand.cpp
 *
 *  Created on: Jul 1, 2010
 *      Author: hpcdjenz
 */
#ifdef STEEREO

#include "estimateRemainingTimeCommand.h"
#include <Simulation.h>
#include <utils/Timer.h>

EstimateRemainingTimeCommand::EstimateRemainingTimeCommand() : SteereoCommand (false, "estimateRemainingTime", 0) {
	// TODO Auto-generated constructor stub
}

EstimateRemainingTimeCommand::~EstimateRemainingTimeCommand() {
	// TODO Auto-generated destructor stub
}

ReturnType EstimateRemainingTimeCommand::execute ()
{
	Simulation* sim = (Simulation*) (this->getData(0));
	Timer* loopTimer = sim->getLoopTimer();
	loopTimer->stop();
	double elapsedTime = loopTimer->get_etime();
	unsigned long currentStep = sim->getSimStep();
	unsigned long numberOfTimesteps = sim->getNumberOfTimeSteps();
	double estimation = ((double) numberOfTimesteps / currentStep) * elapsedTime - elapsedTime;
	SteereoStream stream;
	stream << currentStep << numberOfTimesteps << elapsedTime << estimation;
	logger->debug() << "estimating " << estimation << std::endl;
	if (this->getIntraCommunicator()->amIRoot ())
	{
		this->getCommunicator()->sendStream (this->getConnection()->getDataConnection(), stream);
	}
	loopTimer->start();
	return EXECUTED;

}

void EstimateRemainingTimeCommand::setParameters (std::list<std::string> params)
{

}

SteereoCommand* EstimateRemainingTimeCommand::generateNewInstance ()
{
	return new EstimateRemainingTimeCommand;
}

bool EstimateRemainingTimeCommand::condition ()
{
	return true;
}

#endif /* STEEREO */
