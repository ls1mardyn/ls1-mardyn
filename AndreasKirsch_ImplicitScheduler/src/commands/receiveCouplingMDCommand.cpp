/*
 * ReceiveCouplingMDCommand.cpp
 *
 *  Created on: Apr 24, 2009
 *      Author: hpcdjenz
 */

#ifdef STEEREO
#include "receiveCouplingMDCommand.h"
#include "../Simulation.h"

int ReceiveCouplingMDCommand::startStep = 0;

ReceiveCouplingMDCommand::ReceiveCouplingMDCommand() : SteereoCommand(false, "receiveCouplingMD")
{
	// TODO Auto-generated constructor stub

}

ReceiveCouplingMDCommand::~ReceiveCouplingMDCommand()
{
	// TODO Auto-generated destructor stub
}

ReturnType ReceiveCouplingMDCommand::execute ()
{
	return EXECUTED;
}

void ReceiveCouplingMDCommand::setParameters (std::list<std::string> params)
{

}

SteereoCommand* ReceiveCouplingMDCommand::generateNewInstance ()
{
	return new ReceiveCouplingMDCommand;
}

bool ReceiveCouplingMDCommand::condition ()
{
	if ((stepInterval > 0) && (this->getData() != NULL))
	{
		unsigned long actStep = ((Simulation*) this->getData())-> getSimStep ();
		return (((actStep - startStep) % stepInterval) == 0);
	}
	else
	{
		return true;
	}
}

COMMAND_AUTOREG(ReceiveCouplingMDCommand)
#endif /* STEEREO */
