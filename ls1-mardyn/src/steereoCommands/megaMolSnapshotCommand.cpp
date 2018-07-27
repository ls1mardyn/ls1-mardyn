/*
 * megaMolSnapshotCommand.cpp
 *
 *  Created on: Mar 19, 2009
 *      Author: hpcdjenz
 */

#ifdef STEEREO

#include "megaMolSnapshotCommand.h"

#include <iostream>
#include "snapshotCommand.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"
#include <steereo/steereoSteeringBase.h>
#include <steereo/steereoStream.h>
#include <steereo/steereoDefinitions.h>

#ifdef ENABLE_MPI
#include <mpi.h>
#endif



Simulation* MegaMolSnapshotCommand::sim = NULL;
int MegaMolSnapshotCommand::startStep = 0;

MegaMolSnapshotCommand::MegaMolSnapshotCommand () : SteereoCommand (false, "getMegaMolSnapshot", 0)
{
  //this->setCommandName ("getSnapshot");
  stepInterval = 0;
  colouringVal = 0;
}

MegaMolSnapshotCommand::~MegaMolSnapshotCommand ()
{

}

/**
 * Format:
 *
 */

ReturnType MegaMolSnapshotCommand::execute ()
{
	//SteereoLogger::setOutputLevel(4);

  Domain* dom = sim -> getDomain ();
  int numberOfComponents = dom -> getComponents().size();
  SteereoStream processStream;
  SteereoStream* wholeStream;
  SteereoStream partStream[numberOfComponents];

  int molCompNumber[numberOfComponents];


  int counter = 0;
  // at the moment the simulation area extents need 6 floats and minVal and maxVal need 2;
    // for colouring this gives the minimum and maximum value
  float minVal = 1e100, maxVal = 0;
  ParticleContainer* m_molecules = sim -> getMolecules ();

  if (this->getIntraCommunicator()->amIRoot())
  {
  	processStream.allocateMemory(8 * sizeof(float) + 3 * sizeof(int));

  	processStream << (float) 0.0 << (float) 0.0 << (float) 0.0;
  	processStream << (float) dom->getGlobalLength(0) << (float) dom->getGlobalLength(1) << (float) dom->getGlobalLength(2);
  	processStream << (int) this->getIntraCommunicator()->getSize();
  	processStream << (int) numberOfComponents;
  	processStream << (int) colouringVal;
  } else
  {
  	processStream.allocateMemory(2 * sizeof(float));
  }

  logger->debug() << "sending " << numberOfComponents << " components " << std::endl;

  Molecule* pos = NULL;

  for (int i = 0; i < numberOfComponents; i++)
  {
  	molCompNumber[i] = 0;
  }

  for (pos = m_molecules->begin (); pos != m_molecules->end (); pos = m_molecules->next ()) {
  	molCompNumber[pos->componentid()]++;
  }

  for (int i = 0; i < numberOfComponents; i++)
  {
  	logger->debug() << "allocating for " << molCompNumber[i] << " molecules " << std::endl;
  	partStream[i].allocateMemory(molCompNumber[i] * sizeof(float) * (3 + (colouringVal > 1)) + 1 * sizeof(int));
  	logger->debug() << "inserting number of particles in component: " << molCompNumber[i] << std::endl;
  	partStream[i] << molCompNumber[i];
  }

  for (pos = m_molecules->begin (); pos != m_molecules->end (); pos = m_molecules->next ()) {
  	partStream[pos->componentid()] << (float) pos->r(0);
    partStream[pos->componentid()] << (float) pos->r(1);
    partStream[pos->componentid()] << (float) pos->r(2);
    counter++;
  }

  logger->debug() << "colouringVal is " << colouringVal << std::endl;
  if (colouringVal > 1)
  {
  	float colVal = 0;
  	for (pos = m_molecules->begin (); pos != m_molecules->end(); pos = m_molecules->next())
  	{
  		colVal = 0;
  		// Colour by Force
  		if (colouringVal == 2)
  		{
  			colVal = pos->F(0) * pos->F(0) + pos->F(1) * pos->F(1) + pos->F(2) * pos->F(2);
  			colVal = sqrt(colVal);
  		}
			else if (colouringVal == 3)  // colour by velocity
			{
				colVal = sqrt (pos->v2());
			}
  		if (colVal > maxVal)
  		{
  			maxVal = colVal;
  		}
  		if (colVal < minVal)
  		{
  			minVal = colVal;
  		}
  		partStream[pos->componentid()] << (float) colVal;
  	}
  }
  else
  {
  	minVal = 0;
  	maxVal = (float) numberOfComponents - 1;
  }
  processStream << minVal << maxVal;
  for (int i = 0; i < numberOfComponents; i++)
  {
  	logger->debug() << "inserting partStream[" << i << "] of size " << partStream[i].getStreamSize() << std::endl;
  	processStream << partStream[i];
  	logger->debug() << "size of processStream is now " << processStream.getStreamSize() << std::endl;
  }

	this->getIntraCommunicator()->gatherOnRoot (&processStream, &wholeStream);


  if (this->getIntraCommunicator()->amIRoot ())
	{
  	logger->debug() << "snapshotCommand: assembled the data to be sent in a stream: " << wholeStream->getStreamSize()
  	  							  << std::endl;
		this->getCommunicator()->sendStream (this->getConnection()->getDataConnection(), *wholeStream);
	}
  logger->debug() << "megMolSnapshotCommand: sent the stream" << std::endl;

  if (stepInterval > 0)
  {
    return REPETITION_REQUESTED;
  }
  else
  {
    return EXECUTED;
  }


}

void MegaMolSnapshotCommand::setParameters (std::list <std::string> params)
{
  int listSize = params.size ();
  if (listSize > 0)
  {
    stepInterval = atoi (params.front().c_str ());
    params.pop_front();
  }
  if (listSize > 1)
  {
  	colouringVal = atoi (params.front().c_str());
  	params.pop_front();
  }

}

SteereoCommand* MegaMolSnapshotCommand::generateNewInstance ()
{
  SteereoCommand* comm = new MegaMolSnapshotCommand ();
  startStep = sim -> getSimStep ();
  ((MegaMolSnapshotCommand*) comm) -> setSimData (sim);
  return comm;
}

bool MegaMolSnapshotCommand::condition ()
{
  if (stepInterval > 0)
  {
    unsigned long actStep = sim -> getSimStep ();
    return (((actStep - startStep) % stepInterval) == 0);
  }
  else
  {
    return true;
  }
}
#endif
