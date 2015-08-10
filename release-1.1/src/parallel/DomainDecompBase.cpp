#include "parallel/DomainDecompBase.h"

#include <fstream>

#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"

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

void DomainDecompBase::getBoundingBoxMinMax(Domain *domain, double *min, double *max) {
	for(int d = 0; d < 3; d++) {
		min[d] = getBoundingBoxMin(d, domain);
		max[d] = getBoundingBoxMax(d, domain);
	}
}
