/*
 * SendCouplingMDCommand.cpp
 *
 *  Created on: Apr 23, 2009
 *      Author: hpcdjenz
 */

#ifdef STEEREO
#include "sendCouplingMDCommand.h"
#include <cstdio>
#include <cstdlib>
#include <map>
#include "../Simulation.h"
#include "../particleContainer/ParticleContainer.h"
#include "../molecules/Molecule.h"
#include "../Domain.h"

int SendCouplingMDCommand::startStep = 0;
int SendCouplingMDCommand::outmin;
int SendCouplingMDCommand::outmax;
int SendCouplingMDCommand::borderToLook = 0; // {xmax=0, xmin, ymax, ymin, zmax, zmin}

SendCouplingMDCommand::SendCouplingMDCommand() : SteereoCommand (false, "sendCouplingMD")
{
		outmin = 0;
		outmax = 0;
}

SendCouplingMDCommand::~SendCouplingMDCommand()
{
	// TODO Auto-generated destructor stub
}

ReturnType SendCouplingMDCommand::execute ()
{
	std::map<int,int> deleteMap;
	std::cout << "I will now execute the SendCouplingMDCommand" << std::endl;
	std::cout << "borderToLook is " << borderToLook << std::endl;
	ParticleContainer* moleculeContainer = ((Simulation*) this->getData())->getMolecules();
	std::vector<Component> components = ((Simulation*) this->getData())->getDomain()->getComponents();

	double rmin; // lower corner of the process-specific domain //ENABLE_MPI
	double rmax;
	int dim = borderToLook / 2;
	int dir = borderToLook % 2;
	outmin = 0;
	outmax = 0;
	rmin =  moleculeContainer->getBoundingBoxMin(dim);
	rmax =  moleculeContainer->getBoundingBoxMax(dim);

	std::cout << "dim is " << dim << ", dir is " << dir << std::endl;
	std::cout << "halo is " << moleculeContainer->get_halo_L(dim) << std::endl;
	Molecule* currentMolecule;


	double low_limit; // particles below this limit have to be copied or moved to the lower process
	double high_limit; // particles above(or equal) this limit have to be copied or moved to the higher process

	// set limits (outside "inner" region)
	low_limit = rmin;
	high_limit = rmax;
	currentMolecule = moleculeContainer->begin();

	std::cout << "low_limit: " << low_limit << " / high_limit: " << high_limit << std::endl;
	while(currentMolecule!=moleculeContainer->end()){

		const double& rd=currentMolecule->r(dim);
		if ((dir == 0) && (rd > high_limit)){
			//std::cout << "outmax++" << std::endl;
			//moleculeContainer->addParticle(m1);
			outmax++;
			deleteMap[currentMolecule->id()] = 1;
			//std::cout << currentMolecule->oldr(dim) << " -> " << currentMolecule->r(dim) << std::endl;
			currentMolecule = moleculeContainer->deleteCurrent ();
			//currentMolecule = moleculeContainer->next();
		}
		else if ((dir == 1) && (rd < low_limit)){

			deleteMap[currentMolecule->id()] = 1;
			//std::cout << currentMolecule->oldr(dim) << " -> " << currentMolecule->r(dim) << std::endl;
			currentMolecule = moleculeContainer->deleteCurrent();
			outmin++;
		}
		else {
			currentMolecule = moleculeContainer->next();
		}
	}
	std::cout << "outmin["<< dim << "] = " << outmin << std::endl;
	std::cout << "outmax["<< dim << "] = " << outmax << std::endl;

	currentMolecule = moleculeContainer->begin();
	int deleteCount = 0;
	while(currentMolecule!=moleculeContainer->end()){
		int actID = currentMolecule->id();
		int rOld = currentMolecule->oldr(dim);
		if (deleteMap.find(actID) != deleteMap.end())
		{
			//std::cout << currentMolecule->oldr(dim) - currentMolecule->r(dim) << std::endl;

			/* was the molecule inside before ? */
			if (((dir == 0) && (rOld <= high_limit)) ||
					((dir == 1) && (rOld >= low_limit))){
				currentMolecule = moleculeContainer->deleteCurrent();
				deleteCount++;
			}
			else {
				currentMolecule = moleculeContainer->next();
			}
		} else {
			currentMolecule = moleculeContainer->next();
		}
	}

	std::cout << "deleted " << deleteCount << " molecules for real" << std::endl;

	std::cout << "successfully executed SendCouplingCommand " << std::endl;
	if (stepInterval > 0)
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
		stepInterval = atoi (params.front().c_str ());
		params.pop_front();
	}
	if (listSize > 1)
	{
		borderToLook = atoi (params.front().c_str ());
		params.pop_front();
	}
}

SteereoCommand* SendCouplingMDCommand::generateNewInstance ()
{
	return new SendCouplingMDCommand;
}

bool SendCouplingMDCommand::condition ()
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

COMMAND_AUTOREG(SendCouplingMDCommand)

#endif
