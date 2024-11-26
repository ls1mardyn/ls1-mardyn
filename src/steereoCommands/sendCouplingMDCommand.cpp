/*
 * SendCouplingMDCommand.cpp
 *
 *  Created on: Apr 23, 2009
 *      Author: hpcdjenz
 */

#if defined(STEEREO) && defined(STEEREO_COUPLING)
#include "sendCouplingMDCommand.h"
#include <cstdio>
#include <cstdlib>
#include <map>
#include "Simulation.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"
#include "Domain.h"

std::vector<Molecule>* SendCouplingMDCommand::transferContainer=NULL;

SendCouplingMDCommand::SendCouplingMDCommand() : SteereoCouplingCommand (false, "sendCouplingMD")
{
	//startStep = 1;
	outmin = 0;
	outmax = 0;
	//borderToLook = 1;
	borderToLook = 0;
	transferContainer = NULL;
}

SendCouplingMDCommand::~SendCouplingMDCommand()
{
}

/** transfer the accumulated molecules to the coupled simulation.
 * The data format is:
 * //bounding box MD - 3 doubles (first of outflow direction, then the other 2 directions)
 * (Both simulation should have gotten the same boundaries)
 * For each boundary:
 * 	 Number of Particles - 1 int
 *   Per Particle:
 *     mass     - 1 double
 *     position - 3 doubles
 *     velocity - 3 doubles
 *	@return returns EXECUTED, so the command will not be put in the queue again (in transfer mode)
 */

ReturnType SendCouplingMDCommand::executeTransfer ()
{
	SteereoStream resultStream;
	CouplingInformationType* couplingInfo = (CouplingInformationType*) this->getData(0);
	Simulation* theSim = (Simulation*) this->getData(1);
	//ParticleContainer* moleculeContainer = theSim->getMolecules();

	/*int direction = borderToLook / 2;
	int otherdir1 = (direction + 1) % 3;
	int otherdir2 = (direction + 2) % 3;
	resultStream << moleculeContainer->getBoundingBoxMax(direction);
	resultStream << moleculeContainer->getBoundingBoxMax(otherdir1) << moleculeContainer->getBoundingBoxMax(otherdir2);
	logger->debug() << "I will now execute the transfer of the SendCouplingMDCommand" << std::endl;
	logger->debug() << "borderToLook is " << borderToLook << std::endl;*/
	for (int i = 0; i < couplingInfo->numberOfBoundaries; i++)
	{
		logger->debug() << "For boundary " << i << " I will transfer " << transferContainer[i].size() << " molecules" << std::endl;
		std::vector<Molecule>::iterator molIt = transferContainer[i].begin();
		resultStream.allocateMemory(7 * transferContainer[i].size() * sizeof(double) + sizeof(int));
		resultStream << (int) transferContainer[i].size();
		while (molIt != transferContainer[i].end()) {
			resultStream << (theSim->getDomain()->getComponents())[molIt->componentid()].m();
			resultStream << molIt->r(0) << molIt->r(1) << molIt->r(2) << molIt->v(0) << molIt->v(1) << molIt->v(2);
			molIt++;
		}
		transferContainer[i].clear();
	}
	this->getCommunicator()->sendStream (this->getConnection()->getDataConnection(), resultStream);

	return EXECUTED;
}

ReturnType SendCouplingMDCommand::executeProcessing()
{
	logger->debug() << "starting SendCouplingMDCommand::executeProcessing" << std::endl;
	CouplingInformationType* couplingInfo = (CouplingInformationType*) this->getData(0);
	Simulation* theSim = (Simulation*) this->getData(1);
	ParticleContainer* moleculeContainer = theSim->getMoleculeContainer();

	if (transferContainer == NULL)
	{
		transferContainer = new std::vector<Molecule>[couplingInfo->numberOfBoundaries];
	}

/*	int dim = borderToLook / 2;
	int dir = borderToLook % 2;*/
	int dim, dir;
	outmin = 0;
	outmax = 0;
	double rmin =  moleculeContainer->getBoundingBoxMin(dim);
	double rmax =  moleculeContainer->getBoundingBoxMax(dim);

	logger->debug() << "dim is " << dim << ", dir is " << dir << std::endl;
	logger->debug() << "halo is " << moleculeContainer->get_halo_L(dim) << std::endl;
	Molecule* currentMolecule;

	double low_limit = rmin; // particles below this limit have to be copied or moved to the lower process
	double high_limit = rmax; // particles above(or equal) this limit have to be copied or moved to the higher process

	currentMolecule = moleculeContainer->begin();

	logger->debug() << "low_limit: " << low_limit << " / high_limit: " << high_limit << std::endl;
	while(currentMolecule!=moleculeContainer->end()){
		for (int i = 0; i < couplingInfo->numberOfBoundaries; i++)
		{
			CouplingBoundary currentBoundary = couplingInfo->boundaries[i];
			//const double& rd=currentMolecule->r(dim);
			const double& rd = currentMolecule->r(currentBoundary.outFlowDirection);
			const double& ro1 = currentMolecule->r(currentBoundary.otherDirection[0]);
			const double& ro2 = currentMolecule->r(currentBoundary.otherDirection[1]);
			if ((currentBoundary.lowerHigher == 1) && (rd > high_limit) && isInBounds(ro1,ro2, &currentBoundary))
			{
				outmax++;
				transferContainer[i].push_back(*currentMolecule);
				currentMolecule = moleculeContainer->deleteCurrent ();
				break;
			}
			else if ((currentBoundary.lowerHigher == 0) && (rd < low_limit) && isInBounds(ro1,ro2, &currentBoundary))
			{
				transferContainer[i].push_back(*currentMolecule);
				currentMolecule = moleculeContainer->deleteCurrent();
				outmin++;
				break;
			}
			else
			{
				currentMolecule = moleculeContainer->next();
			}
		}
	}
	logger->debug() << "outmin["<< dim << "] = " << outmin << std::endl;
	logger->debug() << "outmax["<< dim << "] = " << outmax << std::endl;
	//logger->debug() << "there are now " << transferContainer.size() << " molecules to transfer " << std::endl;
	if (getStepInterval() > 0)
	{
		return REPETITION_REQUESTED;
	}
	else
	{
		return EXECUTED;
	}

}

void SendCouplingMDCommand::setParameters (std::list<std::string> params)
{
	int listSize = params.size ();
	if (listSize > 0)
	{
		borderToLook = atoi (params.front().c_str ());
		params.pop_front();
		logger->debug() << "borderToLook is set to " << borderToLook << std::endl;
	}
}

SteereoCommand* SendCouplingMDCommand::generateNewInstance ()
{
	return new SendCouplingMDCommand;
}

bool SendCouplingMDCommand::condition ()
{
	return this->isReady ();
}

bool SendCouplingMDCommand::isInBounds (double ro1, double ro2, CouplingBoundary* boundary)
{
	bool returnValue;
	double borderO1[2], borderO2[2];
	borderO1[0] = boundary->startMD[boundary->otherDirection[0]];
	borderO1[1] = borderO1[0] + boundary->dimMD[boundary->otherDirection[0]];
	borderO2[0] = boundary->startMD[boundary->otherDirection[1]];
	borderO2[1] = borderO2[0] + boundary->dimMD[boundary->otherDirection[1]];
	returnValue = (ro1 >= borderO1[0]) && (ro1 <= borderO1[1]);
	returnValue = returnValue && (ro2 >= borderO2[0]) && (ro2 <= borderO2[1]);
	return returnValue;
}



COMMAND_AUTOREG(SendCouplingMDCommand)

#endif
