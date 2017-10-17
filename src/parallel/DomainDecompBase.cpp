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
	const double shiftMagnitude = moleculeContainer->getBoundingBoxMax(dim) - moleculeContainer->getBoundingBoxMin(dim);

	// molecules that have crossed the lower boundary need a positive shift
	// molecules that have crossed the higher boundary need a negative shift
	// loop over -+1 for dim=0, -+2 for dim=1, -+3 for dim=2
	const int sDim = dim+1;
	for(int direction = -sDim; direction < 2*sDim; direction += 2*sDim) {
		double shift = copysign(shiftMagnitude, static_cast<double>(-direction));

		double startRegion[3];
		double endRegion[3];

		moleculeContainer->getHaloRegionPerDirection(direction, &startRegion, &endRegion);

		#if defined (_OPENMP)
		#pragma omp parallel shared(startRegion, endRegion)
		#endif
		{
			RegionParticleIterator begin = moleculeContainer->iterateRegionBegin(startRegion, endRegion);
			RegionParticleIterator end = moleculeContainer->iterateRegionEnd();

			//traverse and gather all halo particles in the cells
			for(RegionParticleIterator i = begin; i != end; ++i){
				Molecule m = *i;
				m.setr(dim, m.r(dim) + shift);
				// some additional shifting to ensure that rounding errors do not hinder the correct placement
				if (shift < 0) {  // if the shift was negative, it is now in the lower part of the domain -> min
					if (m.r(dim) <= moleculeContainer->getBoundingBoxMin(dim)) { // in the lower part it was wrongly shifted if
						m.setr(dim, moleculeContainer->getBoundingBoxMin(dim));  // ensures that r is at least the boundingboxmin
					}
				} else {  // shift > 0
					if (m.r(dim) >= moleculeContainer->getBoundingBoxMax(dim)) { // in the lower part it was wrongly shifted if
						// std::nexttoward: returns the next bigger value of _boundingBoxMax
						vcp_real_calc r = moleculeContainer->getBoundingBoxMax(dim);
						m.setr(dim, std::nexttoward(r, r - 1.f));  // ensures that r is smaller than the boundingboxmax
					}
				}
				moleculeContainer->addParticle(m);
				i.deleteCurrentParticle(); //removeFromContainer = true;
			}
		}
	}
}

void DomainDecompBase::populateHaloLayerWithCopies(unsigned dim, ParticleContainer* moleculeContainer) const {
	double shiftMagnitude = moleculeContainer->getBoundingBoxMax(dim) - moleculeContainer->getBoundingBoxMin(dim);

	// molecules that have crossed the lower boundary need a positive shift
	// molecules that have crossed the higher boundary need a negative shift
	// loop over -+1 for dim=0, -+2 for dim=1, -+3 for dim=2
	const int sDim = dim+1;
	for(int direction = -sDim; direction < 2*sDim; direction += 2*sDim) {
		double shift = copysign(shiftMagnitude, static_cast<double>(-direction));

		double startRegion[3];
		double endRegion[3];

		moleculeContainer->getBoundaryRegionPerDirection(direction, &startRegion, &endRegion);

		#if defined (_OPENMP)
		#pragma omp parallel shared(startRegion, endRegion)
		#endif
		{
			RegionParticleIterator begin = moleculeContainer->iterateRegionBegin(startRegion, endRegion);
			RegionParticleIterator end = moleculeContainer->iterateRegionEnd();

			//traverse and gather all boundary particles in the cells
			for(RegionParticleIterator i = begin; i != end; ++i){
				Molecule m = *i;
				m.setr(dim, m.r(dim) + shift);
				// checks if the molecule has been shifted to inside the domain due to rounding errors.
				if (shift < 0) {  // if the shift was negative, it is now in the lower part of the domain -> min
					if (m.r(dim) >= moleculeContainer->getBoundingBoxMin(dim)) { // in the lower part it was wrongly shifted if
						vcp_real_calc r = moleculeContainer->getBoundingBoxMin(dim);
						m.setr(dim, std::nexttoward(r, r - 1.f));  // ensures that r is smaller than the boundingboxmin
					}
				} else {  // shift > 0
					if (m.r(dim) < moleculeContainer->getBoundingBoxMax(dim)) { // in the lower part it was wrongly shifted if
						// std::nextafter: returns the next bigger value of _boundingBoxMax
						vcp_real_calc r = moleculeContainer->getBoundingBoxMax(dim);
						m.setr(dim, std::nexttoward(r, r + 1.f));  // ensures that r is bigger than the boundingboxmax
					}
				}
				moleculeContainer->addHaloParticle(m);
			}
		}
	}
}

int DomainDecompBase::getNonBlockingStageCount(){
	return -1;
}

bool DomainDecompBase::queryBalanceAndExchangeNonBlocking(bool /*forceRebalancing*/, ParticleContainer* /*moleculeContainer*/, Domain* /*domain*/, double etime){
	return false;
}

void DomainDecompBase::balanceAndExchange(double /*lastTraversalTime*/, bool /* forceRebalancing */, ParticleContainer* moleculeContainer, Domain* domain) {
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

double DomainDecompBase::getBoundingBoxMin(int /*dimension*/, Domain* /*domain*/) {
	return 0.0;
}

double DomainDecompBase::getBoundingBoxMax(int dimension, Domain* domain) {
	return domain->getGlobalLength(dimension);
}

double DomainDecompBase::getTime() const {
	return double(clock()) / CLOCKS_PER_SEC;
}

unsigned DomainDecompBase::Ndistribution(unsigned localN, float* minrnd, float* maxrnd) {
	*minrnd = 0.0;
	*maxrnd = 1.0;
	return localN;
}

void DomainDecompBase::assertIntIdentity(int /* IX */) {
}

void DomainDecompBase::assertDisjunctivity(TMoleculeContainer* /* mm */) const {
}

void DomainDecompBase::printDecomp(std::string /* filename */, Domain* /* domain */) {
	global_log->warning() << "printDecomp useless in serial mode" << std::endl;
}

int DomainDecompBase::getRank() const {
	return _rank;
}

int DomainDecompBase::getNumProcs() const {
	return _numProcs;
}

void DomainDecompBase::barrier() const {
}

void DomainDecompBase::writeMoleculesToFile(std::string filename, ParticleContainer* moleculeContainer, bool binary) const{
	for (int process = 0; process < getNumProcs(); process++) {
		if (getRank() == process) {
			std::ofstream checkpointfilestream;
			if(binary == true){
//				checkpointfilestream.open((filename + ".xdr").c_str(), std::ios::binary | std::ios::out | std::ios::trunc);
				checkpointfilestream.open((filename + ".xdr").c_str(), std::ios::binary | std::ios::out | std::ios::app);
			}
			else {
				checkpointfilestream.open(filename.c_str(), std::ios::app);
				checkpointfilestream.precision(20);
			}

			ParticleIterator tempMolecule;
			for (tempMolecule = moleculeContainer->iteratorBegin(); tempMolecule != moleculeContainer->iteratorEnd(); ++tempMolecule) {
				if(binary == true){
					tempMolecule->writeBinary(checkpointfilestream);
				}
				else {
					tempMolecule->write(checkpointfilestream);
				}
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

void DomainDecompBase::collCommAllreduceSumAllowPrevious() {
	_collCommBase.allreduceSum();
}

void DomainDecompBase::collCommAllreduceCustom(ReduceType type) {
	_collCommBase.allreduceCustom(type);
}

void DomainDecompBase::collCommScanSum() {
	_collCommBase.scanSum();
}

void DomainDecompBase::collCommBroadcast(int root) {
	_collCommBase.broadcast(root);
}

double DomainDecompBase::getIOCutoffRadius(int dim, Domain* domain,
		ParticleContainer* moleculeContainer) {

	double length = domain->getGlobalLength(dim);
	double cutoff = moleculeContainer->getCutoff();
	mardyn_assert( ((int) length / cutoff ) == length / cutoff );
	return cutoff;
}



