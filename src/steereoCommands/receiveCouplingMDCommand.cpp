/*
 * ReceiveCouplingMDCommand.cpp
 *
 *  Created on: Apr 24, 2009
 *      Author: hpcdjenz
 */

#if defined(STEEREO) && defined(STEEREO_COUPLING)
#include "receiveCouplingMDCommand.h"
#include "../Simulation.h"
#include "../Domain.h"
#include "../particleContainer/ParticleContainer.h"
#include "../molecules/Molecule.h"
#include <cmath>

Molecule* ReceiveCouplingMDCommand::referenceMolecule = NULL;

ReceiveCouplingMDCommand::ReceiveCouplingMDCommand() : SteereoCouplingCommand(false, "receiveCouplingMD")
{
	 srand( time(NULL) );
	 //receiveBuffers = NULL;
	 stepOverlap = 0.3;
	 potentialEnergyOverlap = 1.5;
	 minimalStepSize = 0.01;
	 accuracy = 0.05;
	 maximumIterations = 30;
	 maximumTries = 5;
	 maximumDistance = 1.;
	 numberOfWorkpackagesPerCycle = 10;
	 targetEnergy = 0.001;
}

ReceiveCouplingMDCommand::~ReceiveCouplingMDCommand()
{
	// TODO Auto-generated destructor stub
}

ReturnType ReceiveCouplingMDCommand::executeProcessing()
{
	logger->debug() << "I have " << workPackages.size() << " molecules ready for insertion "<< std::endl;
	logger->info() << "I try to insert " << numberOfWorkpackagesPerCycle << std::endl;
	for (int i = 0; i < numberOfWorkpackagesPerCycle; i++)
	{
		if (workPackages.size() > 0 )
		{
			DataPackage currentPackage = workPackages.front();
			workPackages.pop();
			int success = insertParticle(currentPackage.translatedPos, maximumDistance, targetEnergy, &currentPackage.boundary);
			if (!success)
			{
				logger->info() << "failed to insert particle" << std::endl;
				workPackages.push (currentPackage);
			}
			else
			{
				logger->info() << "succesfully inserted particle " << std::endl;
			}
		}
	}

	if (getStepInterval() > 0)
	{
		return REPETITION_REQUESTED;
	}
	else
	{
		return EXECUTED;
	}
}


/** receive the accumulated distribution data to the coupled simulation.
 * The data format is:
 * Number of vertices in outflowing direction - 1 int (usually the value is 5)
 * Per boundary (the receiving command should already know the boundary data):
 *   The buffer:
 *	   per cell:
 *	 	   Distribution values in outflowing direction - NumberOfVertices (s.a.) doubles
 * @return returns EXECUTED, so the command will not be put in the queue again (in transfer mode)
 */
ReturnType ReceiveCouplingMDCommand::executeTransfer()
{
	CouplingInformationType* couplingInfo = (CouplingInformationType*) this->getData(0);
	Simulation* theSim = (Simulation*) this->getData(1);
	Domain* theDomain = theSim->getDomain();
	//ParticleContainer* particles = theSim->getMolecules ();
	SteereoStream recvStream;
	logger->debug() << "receive data from coupling partner" << std::endl;
	this->getCommunicator()->receiveStream(this->getConnection()->getDataConnection(), &recvStream);
	logger->debug() << "I received " << recvStream.getStreamSize() << " byte " << std::endl;
	int verticesInDirection;
	//int otherDir[2];

	//int direction = borderToLook / 2;
	//int maxMin = borderToLook % 2;
	//otherDir[0] = (direction + 1) % 3;
	//otherDir[1] = (direction + 1) % 3;
	/*recvStream >> directionDimension;
	recvStream >> borderDimensions[0];
	recvStream >> borderDimensions[1];*/
	recvStream >> verticesInDirection;
	//int bufferSize = borderDimensions[0] * borderDimensions[1];
	int bufferSize = 0;
	//logger->debug() << "my dimensions: " << directionDimension << "; " <<
	//			borderDimensions[0] << " " << borderDimensions[1] << std::endl;
	int dataSize = recvStream.getRestStreamSize() / (verticesInDirection * sizeof(double));

	if (receiveBuffers.size() == 0)
	{
		for (int i = 0; i < couplingInfo->numberOfBoundaries; i++)
		{
			CouplingBoundary tempBoundary= couplingInfo->boundaries[i];
			bufferSize = tempBoundary.dimLB[tempBoundary.otherDirection[0]] *
					tempBoundary.dimLB[tempBoundary.otherDirection[1]];
			double* tempBuffer = new double[bufferSize];
			memset (tempBuffer, '\0', bufferSize * sizeof(double));
			receiveBuffers.push_back (tempBuffer);
		}
	}

	for (int i = 0; i < couplingInfo->numberOfBoundaries; i++)
	{
		CouplingBoundary currentBoundary = couplingInfo->boundaries[i];
		double* myReceiveBuffer = receiveBuffers[i];
		//std::cout << "I received " << dataSize << " (partial) cells" << std::endl;
		/*rhoTranslation = 1.0;
		for (int i = 0; i < 2; i++) {
			unitTranslation[i] = myBoxSize[otherDir[i]] / borderDimensions[i];
			rhoTranslation *= unitTranslation[i];
		}
		rhoTranslation *= myBoxSize[direction] / directionDimension;*/
		DataPackage tempPackage;
		double referenceMass = (theDomain->getComponents())[referenceMolecule->componentid()].m();
		//referenceMolecule->componentid()

		// loop over the border cells (without halo)

		int otherDir[2];
		otherDir[0] = currentBoundary.otherDirection[0];
		otherDir[1] = currentBoundary.otherDirection[1];
		tempPackage.boundary = currentBoundary;
		tempPackage.translatedPos[currentBoundary.outFlowDirection] = (currentBoundary.lowerHigher == 0) ? 0 : theDomain->getGlobalLength(currentBoundary.outFlowDirection);
		int bufferIndex = 0;
		for (int o1 = 0; o1 < currentBoundary.dimLB[otherDir[0]]; o1++)
		{
			tempPackage.translatedPos[otherDir[0]] = o1 / currentBoundary.lengthScaleMDToLBM[otherDir[0]];
			for (int o2 = 0; o2 < currentBoundary.dimLB[otherDir[1]]; o2++) {
				tempPackage.translatedPos[otherDir[1]] = o2 / currentBoundary.lengthScaleMDToLBM[otherDir[1]];
				//tempPackage.dist = new double[verticesInDirection];
				double tempDist;
				for (int v = 0; v < verticesInDirection; v++) {
					//				recvStream >> tempPackage.dist[i];
					recvStream >> tempDist;
					receiveBuffers[i][bufferIndex] += tempDist;
				}
				int numberOfPackages = (receiveBuffers[i][bufferIndex] / currentBoundary.rhoScalingMDToLBM) / referenceMass;
				receiveBuffers[i][bufferIndex] -= (numberOfPackages * referenceMass) * currentBoundary.rhoScalingMDToLBM;
				for (int i = 0; i < numberOfPackages; i++) {
					workPackages.push(tempPackage);
				}
				bufferIndex++;
			}
		}
	}
	return EXECUTED;
}

void ReceiveCouplingMDCommand::setParameters (std::list<std::string> params)
{
	int listSize = params.size ();
	if (listSize > 0)
	{
/*		borderToLook = atoi (params.front ().c_str ());
		params.pop_front ();
		logger->debug() << "borderToLook is set to " << borderToLook << std::endl;*/
	}
}



SteereoCommand* ReceiveCouplingMDCommand::generateNewInstance ()
{
	return new ReceiveCouplingMDCommand;
}

bool ReceiveCouplingMDCommand::condition ()
{
	return this->isReady();
}

// ----------------------------- private functions -----------------------------

double ReceiveCouplingMDCommand::getPotentialEnergyAndForce (Molecule* victim)
{
	Simulation* theSim = (Simulation*) this->getData();
	//Domain* domain = (Domain*) theSim->getDomain();
	ParticleContainer* molecules = theSim->getMolecules();
	double potEnergy = molecules->getEnergy(victim);
	victim->calcFM();
	return potEnergy;
}

double ReceiveCouplingMDCommand::getNewStepSize (double currentStepSize, double potentialEnergy, double targetEnergy, double force)
{
	if (potentialEnergy > potentialEnergyOverlap)
	{
		return stepOverlap;
	}
	else
	{
		double proposedSize = (potentialEnergy - targetEnergy) / force;
		if (proposedSize < minimalStepSize)
		{
			proposedSize = minimalStepSize;
		}
		return proposedSize;
	}
}


void ReceiveCouplingMDCommand::updatePosition (Molecule* victim, double stepSize, double* forceModulus)
{
	double tempForce[3];
	*forceModulus = 0;
	for (int i = 0; i < 3; i++)
	{
		tempForce[i] = victim->F(i);
		*forceModulus += tempForce[i] * tempForce[i];
	}
	*forceModulus = sqrt (*forceModulus);
	for (int i = 0; i < 3; i++)
	{
		victim->move (i, (tempForce[i] / *forceModulus) * stepSize);
//		position[i] += (tempForce[i] / forceModulus) * stepSize;
	}
}

void setToRandomPosition (Molecule* victim, int direction, int maxMin, double maxDist)
{
	int otherDir[2] = {(direction + 1) % 3, (direction + 2) % 3};
	double randomDist = ((rand() / (double) RAND_MAX) * maxDist);
	// phi \in [0, 2*pi]
	double randomPhi = ((rand() / (double) RAND_MAX) * 2 * M_PI);
	// theta \in [0, pi/2]
	double randomTheta = ((rand() / (double) RAND_MAX) * (M_PI/2.));
	// place randomly in the half-sphere in direction around the current position of victim
	victim->move (otherDir[0], randomDist * sin (randomTheta) * cos(randomPhi));
	victim->move (otherDir[1], randomDist * sin (randomTheta) * sin(randomPhi));
	victim->move (direction, (maxMin == 0 ? 1 : -1) * randomDist * cos(randomTheta));
}

double distFromStartSqr (double p1[3], Molecule* victim)
{
	double result = 0.0;
	double temp = 0.0;
	for (int i = 0; i < 3; i++)
	{
		temp = victim->r(i) - p1[i];
		result += temp*temp;
	}
	return result;
}

int ReceiveCouplingMDCommand::insertParticle (double startPos[3], double maxDist, double targetEnergy, CouplingBoundary* boundary)
{
	Simulation* theSim = (Simulation*) this->getData(1);
	//int direction = borderToLook / 2;
	int direction = boundary->outFlowDirection;
	//int maxMin = borderToLook % 2;
	double stepLength = 1.0;
	double relativeDifference = 99.9;
	double forceModulus;
	double posDiff = 0.0;
	double maxDistSqr = maximumDistance * maximumDistance;
	Molecule* newMol = new Molecule(*referenceMolecule);
	unsigned long molID = theSim->getMaxID() + 1;
	theSim->setMaxID(molID);
	newMol->setid(molID);
	double nullArray[3] = {0,0,0};
	newMol->setF (nullArray);
	newMol->setM (nullArray);
	for (int i = 0; i < 3; i++)
	{
		newMol->setr(i, startPos[i]);
	}
	int step = 0;
	int currentTry = 0;
	double potentialEnergy = 0.0;

	while ((currentTry < maximumTries) && (relativeDifference > accuracy))
	{
		newMol->setF (nullArray);
		newMol->setM (nullArray);
		setToRandomPosition(newMol, direction, boundary->lowerHigher, maximumDistance / 2.);
		while ((posDiff < maxDistSqr) && (step < maximumIterations) && (relativeDifference > accuracy))
		{
			potentialEnergy = getPotentialEnergyAndForce(newMol);
			updatePosition(newMol, stepLength, &forceModulus);
			stepLength = getNewStepSize(stepLength, potentialEnergy, targetEnergy, forceModulus);
			relativeDifference = abs ((potentialEnergy - targetEnergy) / targetEnergy);
			posDiff = distFromStartSqr(startPos, newMol);
			step++;
		}
		currentTry++;
	}
	if (relativeDifference > accuracy)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}



COMMAND_AUTOREG(ReceiveCouplingMDCommand)
#endif /* STEEREO */
