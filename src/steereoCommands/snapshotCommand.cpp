#ifdef STEEREO


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



Simulation* SnapshotCommand::sim = NULL;
int SnapshotCommand::startStep = 0;

SnapshotCommand::SnapshotCommand () : SteereoCommand (false, "getSnapshot", 0)
{
  //this->setCommandName ("getSnapshot");
  stepInterval = 0;
  sendForces = false;
  sendVelocity = false;
  sendV2 = false;
}

SnapshotCommand::~SnapshotCommand ()
{

}

ReturnType SnapshotCommand::execute ()
{
#ifdef ENABLE_MPI
  int myRank, mySize;
  MPI_Comm_size (MPI_COMM_WORLD, &mySize);
  MPI_Comm_rank (MPI_COMM_WORLD, &myRank);

#else
  logger->debug() << "----------------->Not Parallel snapshot " << std::endl;
  //return EXECUTED;
#endif


  Domain* dom = sim -> getDomain ();
  int numberOfComponents = dom -> getComponents().size();
  SteereoStream daStream;
  SteereoStream outStream;
  int counter = 0;
#ifdef ENABLE_MPI
  int* molNumbers = NULL;
  int* dataSizes = NULL;
  int* displ = NULL;
  int baseDispl = 0;
#endif
  int factor = 3 +  sendForces + sendVelocity + sendV2;
  //int offset = sendV2Max;
  // at the moment the simulation area extents need 6 floats
  int offset = 7;

  ParticleContainer* m_molecules = sim -> getMolecules ();
#ifdef ENABLE_MPI
  int anzahl_mol = sim -> getDomain () -> getglobalNumMolecules ();
#endif
  int local_mol = m_molecules->getNumberOfParticles();
#ifdef ENABLE_MPI
  if (myRank == 0)
  {
    outStream.allocateMemory (factor * anzahl_mol * sizeof(float) + offset * sizeof(float));
  }

  // The stream will only be used for sending. And we will send maximally 3*sizeof(float) * local_mol at once
  daStream.allocateMemory (3 * local_mol * sizeof(float));

#else
  daStream.allocateMemory(factor * local_mol * sizeof(float) + offset * sizeof(float));
#endif

  //send the extents of the domain
  daStream << (float) 0.0 << (float) 0.0 << (float) 0.0;
  daStream << (float) dom->getGlobalLength(0) << (float) dom->getGlobalLength(1) << (float) dom->getGlobalLength(2);
  daStream << (int) numberOfComponents;

#ifdef ENABLE_MPI
  if (myRank == 0)
  {
    molNumbers = new int [mySize];
    dataSizes = new int [mySize];
  }
  MPI_Gather (&local_mol, 1, MPI_INT, molNumbers, 1, MPI_INT, 0, MPI_COMM_WORLD );

  if (myRank == 0)
  {
    displ = new int[mySize];
    displ[0] = 0;
    dataSizes[0] = molNumbers[0] * 3 * sizeof(float);
    for (int i = 1; i < mySize; i++)
    {
      displ[i] = displ[i-1] + dataSizes[i-1];
      dataSizes[i] = molNumbers[i] * 3 * sizeof(float);
      std::cout << "MolNumber ["<< i-1 << "]: " << molNumbers[i-1] << std::endl;
    }
    std::cout << "MolNumber ["<< mySize-1 << "]: " << molNumbers[mySize-1] << std::endl;

  }
#endif

  Molecule* pos = NULL;
  double crit = 0.0;
  for (pos = m_molecules->begin (); pos != m_molecules->end (); pos = m_molecules->next ()) {
    //crit = pos->r(0) * pos->r(0) + pos->r(1) * pos->r(1) + pos->r(2) * pos->r(2);
    daStream << (float) pos->r(0);
    daStream << (float) pos->r(1);
    daStream << (float) pos->r(2);
    counter++;
  }

#ifdef ENABLE_MPI
  MPI_Gatherv (daStream.getStream(), local_mol * 3 * sizeof(float), MPI_CHAR, outStream.getStream (), dataSizes, displ, MPI_CHAR, 0, MPI_COMM_WORLD);

  daStream.resetActPos();
  daStream.resetReadPos();
#endif
  float tempF = 0.0;
  if (sendForces)
  {
    for (pos = m_molecules->begin (); pos != m_molecules->end (); pos = m_molecules->next ())
    {
      tempF = pos->F(0) * pos->F(0) + pos->F(1) * pos->F(1) + pos->F(2) * pos->F(2);
      tempF = sqrt(tempF);
      /*daStream << (float) pos->F(0);
      daStream << (float) pos->F(1);
      daStream << (float) pos->F(2);*/
      daStream << tempF;
    }
#ifdef ENABLE_MPI
    if (myRank == 0)
    {
      baseDispl = displ[mySize-1] + dataSizes[mySize-1];
      displ[0] = baseDispl;
      dataSizes[0] = molNumbers[0] * sizeof(float);
      for (int i = 1; i < mySize; i++)
      {
        displ[i] = displ[i-1] + dataSizes[i-1];
        dataSizes[i] = molNumbers[i] * sizeof(float);
      }
    }
    MPI_Gatherv (daStream.getStream(), local_mol * sizeof(float), MPI_CHAR, outStream.getStream (), dataSizes, displ, MPI_CHAR, 0, MPI_COMM_WORLD);
    daStream.resetActPos();
    daStream.resetReadPos();
#endif
  }
  if (sendVelocity)
  {
    float vel = 0;
    for (pos = m_molecules->begin (); pos != m_molecules->end (); pos = m_molecules->next ())
    {
      vel = sqrt (pos->v2());
      /*daStream << (float) pos->v(0);
      daStream << (float) pos->v(1);
      daStream << (float) pos->v(2);*/
      daStream << vel;
    }
#ifdef ENABLE_MPI
    if (myRank == 0)
    {
      baseDispl = displ[mySize-1] + dataSizes[mySize-1];
      displ[0] = baseDispl;
      dataSizes[0] = molNumbers[0] * sizeof(float);
      for (int i = 1; i < mySize; i++)
      {
        displ[i] = displ[i-1] + dataSizes[i-1];
        dataSizes[i] = molNumbers[i] * sizeof(float);
      }
    }
    MPI_Gatherv (daStream.getStream(), local_mol * sizeof(float), MPI_CHAR, outStream.getStream(), dataSizes, displ, MPI_CHAR, 0, MPI_COMM_WORLD);
    daStream.resetActPos();
    daStream.resetReadPos();
#endif
  }
  float v2max = 0;
  if (sendV2)
  {
    float v2 = 0;
    for (pos = m_molecules->begin (); pos != m_molecules->end (); pos = m_molecules->next ())
    {
      v2 = pos->v2();
      v2max = v2 > v2max ? v2 : v2max;
      daStream << v2;

    }
#ifdef ENABLE_MPI
    if (myRank == 0)
    {
      baseDispl = displ[mySize-1] + dataSizes[mySize-1];
      displ[0] = baseDispl;
      std::cout << "0: " << displ[0] << std::endl;
      dataSizes[0] = molNumbers[0] * sizeof(float);
      std::cout << "datasize[0]: " << dataSizes[0] << std::endl;
      for (int i = 1; i < mySize; i++)
      {
        displ[i] = displ[i-1] + dataSizes[i-1];
        dataSizes[i] = molNumbers[i] * sizeof(float);
        std::cout << "i: " << displ[i] << std::endl;
      }
    }
    MPI_Gatherv (daStream.getStream(), local_mol * sizeof(float), MPI_CHAR, outStream.getStream(), dataSizes, displ, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
    /*if (sendV2Max)
    {
      daStream << v2max;
    }*/
  }
#ifdef ENABLE_MPI
  if (myRank == 0)
  {
 //   std::cout << "snapshotCommand: assembled the data to be sent in a stream: " << outStream.getStreamSize() <<  " "<< counter <<std::endl;
    this->getCommunicator()->sendStream (this->getConnection()->getDataConnection(), outStream);
   // std::cout << "snapshotCommand: sent the stream" << std::endl;
  }
#else
 // std::cout << "snapshotCommand: assembled the data to be sent in a stream: " << daStream.getStreamSize() <<  " "<< counter <<std::endl;
  this->getCommunicator()->sendStream (this->getConnection()->getDataConnection(), daStream);
 // std::cout << "snapshotCommand: sent the stream" << std::endl;
#endif

  if (stepInterval > 0)
  {
    return REPETITION_REQUESTED;
  }
  else
  {
    return EXECUTED;
  }


}

void SnapshotCommand::setParameters (std::list <std::string> params)
{
  int listSize = params.size ();
  logger->debug() << "set parameters has " << params.size () << " parameters" << std::endl;
  if (listSize > 0)
  {
    stepInterval = atoi (params.front().c_str ());
    params.pop_front();
  }
  if (listSize > 1)
  {
    sendForces = atoi (params.front().c_str()) > 0 ? true : false;
    params.pop_front();
  }
  if (listSize > 2)
  {
    sendVelocity = atoi (params.front().c_str()) > 0 ? true : false;
    params.pop_front ();
  }
  if (listSize > 3)
  {
    sendV2 = atoi (params.front().c_str()) > 0 ? true : false;
    params.pop_front ();
  }
  if (listSize > 4)
  {
    sendV2Max = atoi (params.front().c_str()) > 0 ? true : false;
    params.pop_front ();
  }

  if (listSize > 4)
  {
    logger->debug() << "got Parameters: " << stepInterval << " " << sendForces << " " << sendVelocity << " " << sendV2  << " "
    << sendV2Max <<  std::endl;
  }
}

SteereoCommand* SnapshotCommand::generateNewInstance ()
{
  SteereoCommand* comm = new SnapshotCommand ();
  startStep = sim -> getSimStep ();
  ((SnapshotCommand*) comm) -> setSimData (sim);
  return comm;
}

bool SnapshotCommand::condition ()
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

COMMAND_AUTOREG(SnapshotCommand)

#endif
