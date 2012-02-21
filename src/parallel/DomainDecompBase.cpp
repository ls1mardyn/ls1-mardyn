/*
 * DomainDecompBase.cpp
 *
 *  Created on: Feb 21, 2012
 *      Author: eckhardw
 */
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"
#include <fstream>

void DomainDecompBase::writeMoleculesToFile(std::string filename, ParticleContainer* moleculeContainer) {
	for (int process = 0; process < getNumProcs(); process++) {
		if (getRank() == process) {
			std::ofstream checkpointfilestream(filename.c_str(), std::ios::app);
			checkpointfilestream.precision(20);
			Molecule* tempMolecule;
			for (tempMolecule = moleculeContainer->begin(); tempMolecule != moleculeContainer->end(); tempMolecule = moleculeContainer->next()) {
				tempMolecule->write(checkpointfilestream);
			}
			checkpointfilestream.close();
		}
		barrier();
	}
}
