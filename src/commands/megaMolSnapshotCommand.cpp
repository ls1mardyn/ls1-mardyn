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
#include "../Domain.h"
#include "../datastructures/ParticleContainer.h"
#include "../molecules/Molecule.h"
#include <baseSimSteering.h>
#include <steereoStream.h>
#include <steereoDefinitions.h>
#ifdef PARALLEL
  #include <mpi.h>
#endif



Simulation* MegaMolSnapshotCommand::sim = NULL;
int MegaMolSnapshotCommand::startStep = 0;

MegaMolSnapshotCommand::MegaMolSnapshotCommand () : SteereoCommand (false, "getMegaMolSnapshot", 0)
{
  //this->setCommandName ("getSnapshot");
  stepInterval = 0;
}

MegaMolSnapshotCommand::~MegaMolSnapshotCommand ()
{

}

ReturnType MegaMolSnapshotCommand::execute ()
{

  std::cout << "----------------->Not Parallel snapshot " << std::endl;
  Domain* dom = sim -> getDomain ();
  int numberOfComponents = dom -> getComponents().size();
  SteereoStream daStream;
  SteereoStream outStream;
  SteereoStream compStream[numberOfComponents];
  int molCompNumber[numberOfComponents];

  int counter = 0;
 /* int* molNumbers = NULL;
  int* dataSizes = NULL;
  int* displ = NULL;
  int baseDispl = 0;
  int factor = 3;*/
  // at the moment the simulation area extents need 6 floats
  int offsetFloat = 6;

  // 1 "float" for the number of components
  int offsetInt = 1;

  ParticleContainer* m_molecules = sim -> getMolecules ();

  //int anzahl_mol = sim -> getDomain () -> getglobalNumMolecules ();
  //int local_mol = m_molecules->getNumberOfParticles();

  daStream.allocateMemory(offsetFloat * sizeof(float) + offsetInt * sizeof(int));

  //send the extents of the domain
  daStream << (float) 0.0 << (float) 0.0 << (float) 0.0;
  daStream << (float) dom->getGlobalLength(0) << (float) dom->getGlobalLength(1) << (float) dom->getGlobalLength(2);
  daStream << (int) numberOfComponents;
  std::cout << "sending " << numberOfComponents << " components " << std::endl;
  //std::cout << "in getSnapshot we have "<< local_mol << " particles from " << anzahl_mol  << std::endl;

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
  	std::cout << "allocating " << molCompNumber[i] << " floats " << std::endl;
  	compStream[i].allocateMemory(molCompNumber[i] * sizeof(float) * 3);
  }

  for (pos = m_molecules->begin (); pos != m_molecules->end (); pos = m_molecules->next ()) {
  	compStream[pos->componentid()] << (float) pos->r(0);
    compStream[pos->componentid()] << (float) pos->r(1);
    compStream[pos->componentid()] << (float) pos->r(2);
    counter++;
  }

  std::cout << "snapshotCommand: assembled the data to be sent in a stream: " << daStream.getStreamSize() <<  " "<< counter <<std::endl;
  BaseSimSteering::comm->sendStream (this->getConnection()->getDataConnection(), daStream);
	std::cout << "sending " << numberOfComponents << " component streams" << std::endl;
  for (int i = 0; i < numberOfComponents; i++) {
  	BaseSimSteering::comm->sendStream (this->getConnection()->getDataConnection(), compStream[i]);
  }
  std::cout << "snapshotCommand: sent the stream" << std::endl;

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
  std::cout << "set parameters has " << params.size () << " parameters" << std::endl;
  if (listSize > 0)
  {
    stepInterval = atoi (params.front().c_str ());
    params.pop_front();
  }

}

SteereoCommand* MegaMolSnapshotCommand::generateNewInstance ()
{
	std::cout << "generate new MegaMolSnapshotCommand " << std::endl;
  SteereoCommand* comm = new MegaMolSnapshotCommand ();
  startStep = sim -> getSimStep ();
  ((MegaMolSnapshotCommand*) comm) -> setSimData (sim);
  return comm;
}

bool MegaMolSnapshotCommand::condition ()
{
  //int actStep = sd -> getSteps ();
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
