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
		phaseSpaceSize[d] = rmax[d] - rmin[d];

		// set limits (outside "inner" region)
		low_limit = rmin[d] + halo_L[d];
		high_limit = rmax[d] - halo_L[d];

		for (currentMolecule = moleculeContainer->begin();
			 currentMolecule != moleculeContainer->end();
			 currentMolecule = moleculeContainer->next()) {

			const double rd = currentMolecule->r(d);

			const bool copyToLow = rd < low_limit;
			const bool copyToHigh = rd >= high_limit;

			const bool copy = copyToLow or copyToHigh;

			if (copy) {
				// to copy the molecule across the lower boundary
				// we need to increment it's position,
				// otherwise - decrement it.
				const int sign = copyToLow ? +1 : -1;

				for (unsigned short d2 = 0; d2 < 3; d2++)
					new_position[d2] = currentMolecule->r(d2);
				new_position[d] += sign * phaseSpaceSize[d];

				Component* component = _simulation.getEnsemble()->component(currentMolecule->componentid());
				Molecule m1 = Molecule(currentMolecule->id(), component,
						new_position[0], new_position[1], new_position[2],
						currentMolecule->v(0), currentMolecule->v(1), currentMolecule->v(2),
						currentMolecule->q().qw(), currentMolecule->q().qx(), currentMolecule->q().qy(), currentMolecule->q().qz(),
						currentMolecule->D(0), currentMolecule->D(1), currentMolecule->D(2));
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
