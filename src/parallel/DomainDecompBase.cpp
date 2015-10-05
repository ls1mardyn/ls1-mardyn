#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "Domain.h"
#include "ensemble/EnsembleBase.h"


#include <fstream>

#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"

DomainDecompBase::DomainDecompBase() {
}

DomainDecompBase::~DomainDecompBase() {
}

void DomainDecompBase::readXML(XMLfileUnits& /* xmlconfig */) {
}

void DomainDecompBase::exchangeMolecules(ParticleContainer* moleculeContainer, Domain* domain) {

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

void DomainDecompBase::balanceAndExchange(bool /* balance */, ParticleContainer* moleculeContainer, Domain* domain) {
	exchangeMolecules(moleculeContainer, domain);
}

bool DomainDecompBase::procOwnsPos(double x, double y, double z, Domain* domain) {
	if (x < getBoundingBoxMin(0, domain)
			|| x >= getBoundingBoxMax(0, domain)
			|| y < getBoundingBoxMin(1, domain)
			|| y >= getBoundingBoxMax(1, domain)
			|| z < getBoundingBoxMin(2, domain)
			|| z >= getBoundingBoxMax(2, domain))
		return false;
	else
		return true;
}

double DomainDecompBase::getBoundingBoxMin(int dimension, Domain* domain) {
	return 0.0;
}
double DomainDecompBase::getBoundingBoxMax(int dimension, Domain* domain) {
	return domain->getGlobalLength(dimension);
}


double DomainDecompBase::getTime() {
	return double(clock()) / CLOCKS_PER_SEC;
}

unsigned DomainDecompBase::Ndistribution(unsigned localN, float* minrnd, float* maxrnd) {
	*minrnd = 0.0;
	*maxrnd = 1.0;
	return localN;
}

void DomainDecompBase::assertIntIdentity(int /* IX */) {
}

void DomainDecompBase::assertDisjunctivity(TMoleculeContainer* /* mm */) {
}

void DomainDecompBase::printDecomp(std::string /* filename */, Domain* /* domain */) {
	global_log->warning() << "printDecomp useless in serial mode" << std::endl;
}

int DomainDecompBase::getRank() {
	return 0;
}

int DomainDecompBase::getNumProcs() {
	return 1;
}

void DomainDecompBase::barrier() {
}

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

void DomainDecompBase::collCommInit(int numValues) {
	_collCommBase.init(numValues);
}

void DomainDecompBase::collCommFinalize() {
	_collCommBase.finalize();
}

void DomainDecompBase::collCommAppendInt(int intValue) {
	_collCommBase.appendInt(intValue);
}

void DomainDecompBase::collCommAppendUnsLong(unsigned long unsLongValue) {
	_collCommBase.appendUnsLong(unsLongValue);
}

void DomainDecompBase::collCommAppendFloat(float floatValue) {
	_collCommBase.appendFloat(floatValue);
}

void DomainDecompBase::collCommAppendDouble(double doubleValue) {
	_collCommBase.appendDouble(doubleValue);
}

void DomainDecompBase::collCommAppendLongDouble(long double longDoubleValue) {
	_collCommBase.appendLongDouble(longDoubleValue);
}

int DomainDecompBase::collCommGetInt() {
	return _collCommBase.getInt();
}

unsigned long DomainDecompBase::collCommGetUnsLong() {
	return _collCommBase.getUnsLong();
}

float DomainDecompBase::collCommGetFloat() {
	return _collCommBase.getFloat();
}

double DomainDecompBase::collCommGetDouble() {
	return _collCommBase.getDouble();
}

long double DomainDecompBase::collCommGetLongDouble() {
	return _collCommBase.getLongDouble();
}

void DomainDecompBase::collCommAllreduceSum() {
	_collCommBase.allreduceSum();
}

void DomainDecompBase::collCommBroadcast(int root) {
	_collCommBase.broadcast(root);
}
