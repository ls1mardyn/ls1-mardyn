/*
 * RDFDummyDecomposition.cpp
 *
 * Domain decomposition to be used if RDF boundaries are used.
 * Makes sure no halo layers are created
 *
 *  Created on: Jul 5, 2012
 *      Author: tijana
 */

#include "RDFDummyDecomposition.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "utils/Logger.h"

using namespace std;
using Log::global_log;

RDFDummyDecomposition::RDFDummyDecomposition() {
	// TODO Auto-generated constructor stub

}

void RDFDummyDecomposition::balanceAndExchange(bool balance, ParticleContainer* moleculeContainer, const std::vector<Component>& components, Domain* domain) {
	exchangeMolecules(moleculeContainer, components, domain);

}

RDFDummyDecomposition::~RDFDummyDecomposition() {
	// TODO Auto-generated destructor stub
}

void RDFDummyDecomposition::exchangeMolecules(ParticleContainer* moleculeContainer, const std::vector<Component>& components, Domain* domain) {
	moleculeContainer->deleteOuterParticles();
	/*
	double rmin[3]; // lower corner of the process-specific domain //ENABLE_MPI
	double rmax[3];
	double halo_L[3]; // width of the halo strip //ENABLE_MPI
	for (int i = 0; i < 3; i++) {
		rmin[i] = moleculeContainer->getBoundingBoxMin(i);
		rmax[i] = moleculeContainer->getBoundingBoxMax(i);
		halo_L[i] = moleculeContainer->get_halo_L(i);
	}

	Molecule* currentMolecule;
	// molecules that have to be copied (because of halo), get a new position
	double new_position[3];

	double phaseSpaceSize[3];

	double low_limit; // particles below this limit have to be copied or moved to the lower process
	double high_limit; // particles above(or equal) this limit have to be copied or moved to the higher process

	for (unsigned short d = 0; d < 3; ++d) {


		//cout << "low_limit: " << low_limit << " / high_limit: " << high_limit << endl;
		//cout << "halo_L: " << halo_L[0] << " / " << halo_L[1] << " / " << halo_L[2] << endl;
		//cout << "proc_domain_L: " << proc_domain_L[0] << " / " << proc_domain_L[1] << " / " << proc_domain_L[2] << endl;
		while (currentMolecule != moleculeContainer->end()) {
			const double& rd = currentMolecule->r(d);
			if (rd < rmin[d]) {
				// determine the position for the copy of the molecule
				moleculeContainer->deleteCurrent();

			}
			else if (rd >= rmax[d]) {
				// determine the position for the copy of the molecule
				moleculeContainer->deleteCurrent();
			}
			moleculeContainer->update();
			currentMolecule = moleculeContainer->next();

		}
	}
	*/
}
