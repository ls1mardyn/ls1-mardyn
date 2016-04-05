#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "Domain.h"
#include "ensemble/EnsembleBase.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"

#include <fstream>
#include <cmath>

DomainDecompBase::DomainDecompBase() : _rank(0), _numProcs(1) {
}

DomainDecompBase::~DomainDecompBase() {
}

void DomainDecompBase::readXML(XMLfileUnits& /* xmlconfig */) {
}

void DomainDecompBase::exchangeMolecules(ParticleContainer* moleculeContainer, Domain* /*domain*/) {

	for (unsigned d = 0; d < 3; ++d) {
		handleDomainLeavingParticles(d, moleculeContainer);
	}

	for (unsigned d = 0; d < 3; ++d) {
		populateHaloLayerWithCopies(d, moleculeContainer);
	}
}

void DomainDecompBase::handleDomainLeavingParticles(unsigned dim, ParticleContainer* moleculeContainer) const {

	std::vector<Molecule * > mols;
	std::vector<Molecule *>::iterator it;
	const double shiftMagnitude = moleculeContainer->getBoundingBoxMax(dim) - moleculeContainer->getBoundingBoxMin(dim);

	// molecules that have crossed the lower boundary need a positive shift
	// molecules that have crossed the higher boundary need a negative shift
	// loop over -+1 for dim=0, -+2 for dim=1, -+3 for dim=2
	const int sDim = dim+1;
	for(int direction = -sDim; direction < 2*sDim; direction += 2*sDim) {
		double shift = copysign(shiftMagnitude, static_cast<double>(-direction));
		const bool removeFromContainer = true;
		moleculeContainer->getHaloParticlesDirection(direction, mols, removeFromContainer);

		for (it = mols.begin(); it != mols.end(); ++it) {
			Molecule * m = *it;
			m->setr(dim, m->r(dim) + shift);
			moleculeContainer->addParticlePointer(m);
		}
		mols.clear();
	}
}

void DomainDecompBase::populateHaloLayerWithCopies(unsigned dim, ParticleContainer* moleculeContainer) const {

	std::vector<Molecule * > mols;
	std::vector<Molecule *>::iterator it;
	double shiftMagnitude = moleculeContainer->getBoundingBoxMax(dim) - moleculeContainer->getBoundingBoxMin(dim);

	// molecules that have crossed the lower boundary need a positive shift
	// molecules that have crossed the higher boundary need a negative shift
	// loop over -+1 for dim=0, -+2 for dim=1, -+3 for dim=2
	const int sDim = dim+1;
	for(int direction = -sDim; direction < 2*sDim; direction += 2*sDim) {
		moleculeContainer->getBoundaryParticlesDirection(direction, mols);

		double shift = copysign(shiftMagnitude, static_cast<double>(-direction));
		for (it = mols.begin(); it != mols.end(); ++it) {
			Molecule m = Molecule(**it);
			m.setr(dim, m.r(dim) + shift);
			moleculeContainer->addParticle(m);
		}
		mols.clear();
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
	return _rank;
}

int DomainDecompBase::getNumProcs() {
	return _numProcs;
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
