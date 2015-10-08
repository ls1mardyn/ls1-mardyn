#include "DomainDecompDummy.h"

#include "Domain.h"
#include "ensemble/EnsembleBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "utils/Logger.h"

using namespace std;
using Log::global_log;

DomainDecompDummy::DomainDecompDummy() {
}

DomainDecompDummy::~DomainDecompDummy() {
}

void DomainDecompDummy::readXML(XMLfileUnits& xmlconfig) {
	/* no parameters */
}

void DomainDecompDummy::exchangeMolecules(ParticleContainer* moleculeContainer, Domain* domain) {

	double rmin[3]; // lower corner of the process-specific domain
	double rmax[3]; // upper cornder of the process domain
	double halo_L[3]; // width of the halo stripe for boundary condition handling
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
		phaseSpaceSize[d] = rmax[d] - rmin[d];
		// set limits (outside "inner" region)
		low_limit = rmin[d] + halo_L[d];
		high_limit = rmax[d] - halo_L[d];
		

		/// @todo This loop depends on the fact, that inserted molecules will go before the current Molecule and will not be iterated in the loop! Shoule be made secure for general use.
		for (currentMolecule = moleculeContainer->begin(); currentMolecule != moleculeContainer->end(); currentMolecule = moleculeContainer->next()) {
			const double& rd = currentMolecule->r(d);
			if (rd < low_limit) {
				// determine the position for the copy of the molecule
				// when moving parallel to the coordinate d2 to another process, the
				// local coordinates in d2 change
				for (int d2 = 0; d2 < 3; d2++) {
					new_position[d2] = currentMolecule->r(d2);
				}
				new_position[d] = rd + phaseSpaceSize[d];
				Component* component = _simulation.getEnsemble()->component(currentMolecule->componentid());
				Molecule m1 = Molecule(currentMolecule->id(),component,
				                       new_position[0], new_position[1], new_position[2],
				                       currentMolecule->v(0),currentMolecule->v(1),currentMolecule->v(2),
				                       currentMolecule->q().qw(),currentMolecule->q().qx(),currentMolecule->q().qy(),currentMolecule->q().qz(),
				                       currentMolecule->D(0),currentMolecule->D(1),currentMolecule->D(2));
				moleculeContainer->addParticle(m1);
			}
			else if (rd >= high_limit) {
				// determine the position for the copy of the molecule
				// when moving parallel to the coordinate d2 to another process, the
				// local coordinates in d2 change
				for (unsigned short d2 = 0; d2 < 3; d2++) {
					new_position[d2] = currentMolecule->r(d2);
				}
				new_position[d] = rd - phaseSpaceSize[d];
				Component* component = _simulation.getEnsemble()->component(currentMolecule->componentid());
				Molecule m1 = Molecule(currentMolecule->id(),component,
				                       new_position[0], new_position[1], new_position[2],
				                       currentMolecule->v(0),currentMolecule->v(1),currentMolecule->v(2),
				                       currentMolecule->q().qw(),currentMolecule->q().qx(),currentMolecule->q().qy(),currentMolecule->q().qz(),
				                       currentMolecule->D(0),currentMolecule->D(1),currentMolecule->D(2));
				moleculeContainer->addParticle(m1);
			}
		}
	}
}

void DomainDecompDummy::balanceAndExchange(bool balance, ParticleContainer* moleculeContainer, Domain* domain) {
	exchangeMolecules(moleculeContainer, domain);
}


double DomainDecompDummy::getBoundingBoxMin(int dimension, Domain* domain) {
	return 0.0;
}
double DomainDecompDummy::getBoundingBoxMax(int dimension, Domain* domain) {
	return domain->getGlobalLength(dimension);
}


double DomainDecompDummy::getTime() {
	return double(clock()) / CLOCKS_PER_SEC;
}

unsigned DomainDecompDummy::Ndistribution(unsigned localN, float* minrnd, float* maxrnd) {
	*minrnd = 0.0;
	*maxrnd = 1.0;
	return localN;
}

void DomainDecompDummy::assertIntIdentity(int IX) {
}

void DomainDecompDummy::assertDisjunctivity(TMoleculeContainer* mm) {
}

void DomainDecompDummy::printDecomp(std::string filename, Domain* domain) {
	global_log->warning() << "printDecomp useless in serial mode" << endl;
}
