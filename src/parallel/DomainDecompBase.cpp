#include <cmath>
#include <fstream>

#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "Domain.h"
#include "ensemble/EnsembleBase.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"
#include "utils/mardyn_assert.h"
#include "ZonalMethods/FullShell.h"
#include "ForceHelper.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#include "utils/MPI_Info_object.h"
#endif

DomainDecompBase::DomainDecompBase() : _rank(0), _numProcs(1) {
}

DomainDecompBase::~DomainDecompBase() {
}

void DomainDecompBase::readXML(XMLfileUnits& /* xmlconfig */) {
}

void DomainDecompBase::addLeavingMolecules(std::vector<Molecule>& invalidMolecules,
										   ParticleContainer* moleculeContainer) {
	for (auto& molecule : invalidMolecules) {
		for (auto dim : {0, 1, 2}) {
			auto shift = moleculeContainer->getBoundingBoxMax(dim) - moleculeContainer->getBoundingBoxMin(dim);
			auto r = molecule.r(dim);
			if (r < moleculeContainer->getBoundingBoxMin(dim)) {
				r = r + shift;
				if (r >= moleculeContainer->getBoundingBoxMax(dim)) {
					r = std::nextafter(moleculeContainer->getBoundingBoxMax(dim), -1);
				}
			} else if (r >= moleculeContainer->getBoundingBoxMax(dim)) {
				r = r - shift;
				if (r < moleculeContainer->getBoundingBoxMin(dim)) {
					r = moleculeContainer->getBoundingBoxMin(dim);
				}
			}
			molecule.setr(dim, r);
		}
	}
	moleculeContainer->addParticles(invalidMolecules);
}

void DomainDecompBase::exchangeMolecules(ParticleContainer* moleculeContainer, Domain* domain) {
	if (moleculeContainer->isInvalidParticleReturner()) {
		// autopas mode!
		global_log->debug() << "DDBase: Adding + shifting invalid particles." << std::endl;
		// in case the molecule container returns invalid particles using getInvalidParticlesRef(), we have to handle them directly.
		addLeavingMolecules(moleculeContainer->getInvalidParticlesRef(), moleculeContainer);
		// now use direct scheme to transfer the rest!
		FullShell fs;
		double rmin[3];  // lower corner
		double rmax[3];  // higher corner
		for (int d = 0; d < 3; d++) {
			rmin[d] = getBoundingBoxMin(d, domain);
			rmax[d] = getBoundingBoxMax(d, domain);
		}
		HaloRegion ownRegion = {rmin[0], rmin[1], rmin[2], rmax[0], rmax[1], rmax[2], 0, 0, 0, 0.};
		bool coversWholeDomain[3];
		double cellLengthDummy[3]{};
		global_log->debug() << "DDBase: Populating halo." << std::endl;
		auto haloExportRegions =
			fs.getHaloExportForceImportRegions(ownRegion, moleculeContainer->getCutoff(),
																	coversWholeDomain, cellLengthDummy);
		for (auto haloExportRegion : haloExportRegions) {
			populateHaloLayerWithCopiesDirect(haloExportRegion, moleculeContainer,
											  true /*positionCheck, same as doLeavingExchange*/);
		}
	} else {
	    // default ls1-mode (non-autopas, so linked-cells!)
		for (unsigned d = 0; d < 3; ++d) {
			handleDomainLeavingParticles(d, moleculeContainer);
		}
		for (unsigned d = 0; d < 3; ++d) {
			populateHaloLayerWithCopies(d, moleculeContainer);
		}
	}
}

void DomainDecompBase::exchangeForces(ParticleContainer* moleculeContainer, Domain* domain){
	for (unsigned d = 0; d < 3; ++d) {
		handleForceExchange(d,moleculeContainer);
	}
}

void DomainDecompBase::handleForceExchange(unsigned dim, ParticleContainer* moleculeContainer) const {
	const double shiftMagnitude = moleculeContainer->getBoundingBoxMax(dim) - moleculeContainer->getBoundingBoxMin(dim);

	// direction +1/-1 for dim 0, +2/-2 for dim 1, +3/-3 for dim 2
	//const int direction = dim+1;
	const int sDim = dim + 1;
	for (int direction = -sDim; direction < 2 * sDim; direction += 2 * sDim) {
		double shift = copysign(shiftMagnitude, static_cast<double>(-direction));

		// Loop over all halo particles in the positive direction
		double cutoff = moleculeContainer->getCutoff();
		double startRegion[3]{moleculeContainer->getBoundingBoxMin(0) - cutoff,
							  moleculeContainer->getBoundingBoxMin(1) - cutoff,
							  moleculeContainer->getBoundingBoxMin(2) - cutoff};
		double endRegion[3]{moleculeContainer->getBoundingBoxMax(0) + cutoff,
							moleculeContainer->getBoundingBoxMax(1) + cutoff,
							moleculeContainer->getBoundingBoxMax(2) + cutoff};

		if (direction < 0) {
			endRegion[dim] = moleculeContainer->getBoundingBoxMin(dim);
		} else {
			startRegion[dim] = moleculeContainer->getBoundingBoxMax(dim);
		}

#if defined (_OPENMP)
#pragma omp parallel shared(startRegion, endRegion)
#endif
		{
			auto begin = moleculeContainer->regionIterator(startRegion, endRegion, ParticleIterator::ALL_CELLS);

			double shiftedPosition[3];

			decltype(moleculeContainer->getMoleculeAtPosition(shiftedPosition)) originalPreviousIter{};

			for (auto haloIter = begin; haloIter.isValid(); ++haloIter) {

				// Add force of halo particle to original particle (or other duplicates)
				// that have a distance of -'shiftMagnitude' in the current direction
				shiftedPosition[0] = haloIter->r(0);
				shiftedPosition[1] = haloIter->r(1);
				shiftedPosition[2] = haloIter->r(2);
				shiftedPosition[dim] += shift;
				originalPreviousIter =
					addValuesAndGetIterator(moleculeContainer, shiftedPosition, originalPreviousIter, *haloIter);
			}
		}
	}
}

void DomainDecompBase::handleForceExchangeDirect(const HaloRegion& haloRegion, ParticleContainer* moleculeContainer) const {
	double shift[3];
	for (int dim=0;dim<3;dim++){
		shift[dim] = moleculeContainer->getBoundingBoxMax(dim) - moleculeContainer->getBoundingBoxMin(dim);
		shift[dim] *= haloRegion.offset[dim] * (-1);  // for halo regions the shift has to be multiplied with -1; as the offset is in negative direction of the shift
	}

#if defined (_OPENMP)
#pragma omp parallel
#endif
	{
		auto begin = moleculeContainer->regionIterator(haloRegion.rmin, haloRegion.rmax, ParticleIterator::ALL_CELLS);

		double shiftedPosition[3];

		decltype(moleculeContainer->getMoleculeAtPosition(shiftedPosition)) originalPreviousIter{};

		for (auto haloIter = begin; haloIter.isValid(); ++haloIter) {

			// Add force of halo particle to original particle (or other duplicates)
			// that have a distance of -'shiftMagnitude' in the current direction
			for (int dim = 0; dim < 3; dim++) {
				shiftedPosition[dim] = haloIter->r(dim) + shift[dim];
			}
			originalPreviousIter =
				addValuesAndGetIterator(moleculeContainer, shiftedPosition, originalPreviousIter, *haloIter);
		}
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

		double cutoff = moleculeContainer->getCutoff();
		double startRegion[3]{moleculeContainer->getBoundingBoxMin(0) - cutoff,
							  moleculeContainer->getBoundingBoxMin(1) - cutoff,
							  moleculeContainer->getBoundingBoxMin(2) - cutoff};
		double endRegion[3]{moleculeContainer->getBoundingBoxMax(0) + cutoff,
							moleculeContainer->getBoundingBoxMax(1) + cutoff,
							moleculeContainer->getBoundingBoxMax(2) + cutoff};

		if (direction < 0) {
			endRegion[dim] = moleculeContainer->getBoundingBoxMin(dim);
		} else {
			startRegion[dim] = moleculeContainer->getBoundingBoxMax(dim);
		}

		#if defined (_OPENMP)
		#pragma omp parallel shared(startRegion, endRegion)
		#endif
		{
			auto begin = moleculeContainer->regionIterator(startRegion, endRegion, ParticleIterator::ALL_CELLS);

			//traverse and gather all halo particles in the cells
			for(auto i = begin; i.isValid(); ++i){
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
				moleculeContainer->deleteMolecule(i, false); // removeFromContainer = true;
			}
		}
	}
}

void DomainDecompBase::handleDomainLeavingParticlesDirect(const HaloRegion& haloRegion,
														  ParticleContainer* moleculeContainer,
														  std::vector<Molecule>& invalidParticles) const {
	double shift[3];
	for (int dim = 0; dim < 3; dim++) {
		shift[dim] = moleculeContainer->getBoundingBoxMax(dim) - moleculeContainer->getBoundingBoxMin(dim);
		shift[dim] *= haloRegion.offset[dim] * (-1);  // for halo regions the shift has to be multiplied with -1; as the
													  // offset is in negative direction of the shift
	}

	auto shiftAndAdd = [&moleculeContainer, haloRegion, shift](Molecule& m) {
		if (not m.inBox(haloRegion.rmin, haloRegion.rmax)) {
			global_log->error() << "trying to remove a particle that is not in the halo region" << std::endl;
			Simulation::exit(456);
		}
		for (int dim = 0; dim < 3; dim++) {
			if (shift[dim] != 0) {
				m.setr(dim, m.r(dim) + shift[dim]);
				// some additional shifting to ensure that rounding errors do not hinder the correct
				// placement
				if (shift[dim] < 0) {  // if the shift was negative, it is now in the lower part of the domain -> min
					if (m.r(dim) <=
						moleculeContainer->getBoundingBoxMin(dim)) {  // in the lower part it was wrongly shifted if
						m.setr(dim, moleculeContainer->getBoundingBoxMin(
										dim));  // ensures that r is at least the boundingBoxMin
					}
				} else {  // shift > 0
					if (m.r(dim) >=
						moleculeContainer->getBoundingBoxMax(dim)) {  // in the lower part it was wrongly shifted if
						// std::nexttoward: returns the next bigger value of _boundingBoxMax
						vcp_real_calc r = moleculeContainer->getBoundingBoxMax(dim);
						m.setr(dim, std::nexttoward(r, r - 1.f));  // ensures that r is smaller than the boundingBoxMax
					}
				}
			}
		}
		moleculeContainer->addParticle(m);
	};

	if (moleculeContainer->isInvalidParticleReturner()) {
		// move all particles that will be inserted now to the end of the container
		auto removeBegin = std::partition(invalidParticles.begin(), invalidParticles.end(), [=](const Molecule& m) {
			// if this is true, it will be put in the first part of the partition, if it is false, in the second.
			return not m.inBox(haloRegion.rmin, haloRegion.rmax);
		});
		// now insert all particles that are in the second partition.
		std::for_each(removeBegin, invalidParticles.end(), shiftAndAdd);
		// remove them from the vector.
		invalidParticles.erase(removeBegin, invalidParticles.end());
	} else {
#if defined(_OPENMP)
#pragma omp parallel
#endif
		{
			auto begin =
				moleculeContainer->regionIterator(haloRegion.rmin, haloRegion.rmax, ParticleIterator::ALL_CELLS);

			// traverse and gather all halo particles in the cells
			for (auto i = begin; i.isValid(); ++i) {
				shiftAndAdd(*i);
				moleculeContainer->deleteMolecule(i, false);  // removeFromContainer = true;
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
	double interactionLength = moleculeContainer->getCutoff();

	for(int direction = -sDim; direction < 2*sDim; direction += 2*sDim) {
		double shift = copysign(shiftMagnitude, static_cast<double>(-direction));

		double startRegion[3]{moleculeContainer->getBoundingBoxMin(0) - interactionLength,
							  moleculeContainer->getBoundingBoxMin(1) - interactionLength,
							  moleculeContainer->getBoundingBoxMin(2) - interactionLength};
		double endRegion[3]{moleculeContainer->getBoundingBoxMax(0) + interactionLength,
							moleculeContainer->getBoundingBoxMax(1) + interactionLength,
							moleculeContainer->getBoundingBoxMax(2) + interactionLength};

		if (direction < 0) {
			startRegion[dim] = moleculeContainer->getBoundingBoxMin(dim);
			endRegion[dim] = moleculeContainer->getBoundingBoxMin(dim) + interactionLength;
		} else {
			startRegion[dim] = moleculeContainer->getBoundingBoxMax(dim) - interactionLength;
			endRegion[dim] = moleculeContainer->getBoundingBoxMax(dim);
		}


		#if defined (_OPENMP)
		#pragma omp parallel shared(startRegion, endRegion)
		#endif
		{
			auto begin = moleculeContainer->regionIterator(startRegion, endRegion, ParticleIterator::ALL_CELLS);

			//traverse and gather all boundary particles in the cells
			for(auto i = begin; i.isValid(); ++i){
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

void DomainDecompBase::populateHaloLayerWithCopiesDirect(const HaloRegion& haloRegion,
														 ParticleContainer* moleculeContainer, bool positionCheck) const {
	double shift[3];
	for (int dim = 0; dim < 3; dim++) {
		shift[dim] = moleculeContainer->getBoundingBoxMax(dim) - moleculeContainer->getBoundingBoxMin(dim);
		shift[dim] *= haloRegion.offset[dim] * (-1);  // for halo regions the shift has to be multiplied with -1; as the
													  // offset is in negative direction of the shift
	}

#if defined (_OPENMP)
#pragma omp parallel
#endif
	{
		auto begin = moleculeContainer->regionIterator(haloRegion.rmin, haloRegion.rmax, ParticleIterator::ONLY_INNER_AND_BOUNDARY);

		//traverse and gather all boundary particles in the cells
		for (auto i = begin; i.isValid(); ++i) {
			Molecule m = *i;
			for (int dim = 0; dim < 3; dim++) {
				if (shift[dim] != 0) {
					m.setr(dim, m.r(dim) + shift[dim]);
					if (positionCheck) {
						// checks if the molecule has been shifted to inside the domain due to rounding errors.
						if (shift[dim] < 0) {
							// if the shift was negative, it is now in the lower part of the domain -> min
							if (m.r(dim) >= moleculeContainer->getBoundingBoxMin(dim)) {
								// in the lower part it was wrongly shifted if it is at least boxMin
								vcp_real_calc r = moleculeContainer->getBoundingBoxMin(dim);
								m.setr(dim, std::nexttoward(
												r, r - 1.f));  // ensures that r is smaller than the boundingboxmin
							}
						} else {  // shift > 0
							if (m.r(dim) < moleculeContainer->getBoundingBoxMax(dim)) {
								// in the lower part it was wrongly shifted if
								// std::nextafter: returns the next bigger value of _boundingBoxMax
								vcp_real_calc r = moleculeContainer->getBoundingBoxMax(dim);
								m.setr(dim, r);  // ensures that r is at least boundingboxmax
							}
						}
					}
				}
			}
			moleculeContainer->addHaloParticle(m);
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
	return !(x < getBoundingBoxMin(0, domain) || x >= getBoundingBoxMax(0, domain) ||
			 y < getBoundingBoxMin(1, domain) || y >= getBoundingBoxMax(1, domain) ||
			 z < getBoundingBoxMin(2, domain) || z >= getBoundingBoxMax(2, domain));
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

void DomainDecompBase::assertDisjunctivity(ParticleContainer* /* moleculeContainer */) const {
}

void DomainDecompBase::printDecomp(const std::string &filename, Domain *domain, ParticleContainer *particleContainer) {
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
#ifdef ENABLE_MPI
void DomainDecompBase::writeMoleculesToMPIFileBinary(const std::string& filename, ParticleContainer* moleculeContainer) const {
	int rank = getRank();

	MPI_File mpifh;
	MPI_Info_object mpiinfo;
	auto extfilename = filename + ".dat";
	MPI_File_open(MPI_COMM_WORLD, extfilename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, mpiinfo,
	              &mpifh);

	uint64_t numParticles_local = 0;
	uint64_t numParticles_exscan = 0;
	auto begin = moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
	for (auto it = begin; it.isValid(); ++it) numParticles_local++;

	MPI_Exscan(&numParticles_local, &numParticles_exscan, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

	uint16_t particle_data_size = 0;
	// if no particle is found (begin is not valid) particle_data_size is zero and provides no problem.
	if(begin.isValid())
	{
		std::stringstream str;
		begin->writeBinary(str);
		particle_data_size = str.str().size();
	}
	uint64_t buffer_size = 32768;
	std::string __dummy(buffer_size, '\0');  // __dummy for preallocation of internal buffer with buffer_size.
	std::ostringstream write_buffer(__dummy, ios_base::binary);
	__dummy.clear();
	__dummy.shrink_to_fit();
	//char* write_buffer = new char[buffer_size];
	uint64_t offset = numParticles_exscan * particle_data_size;
	MPI_File_seek(mpifh, offset, MPI_SEEK_SET);
	uint64_t buffer_pos = 0;

	for (auto it = begin; it.isValid(); ++it) {
		it->writeBinary(write_buffer);
		buffer_pos += particle_data_size;

		if (buffer_pos > buffer_size - particle_data_size) {
			// we cannot add any more particles to this buffer, so we write the buffer.
			MPI_File_write(mpifh, write_buffer.str().c_str(), buffer_pos, MPI_BYTE, MPI_STATUS_IGNORE);
			// reset buffer position and clear stream.
			buffer_pos = 0;
			write_buffer.str("");
		}
	}
	MPI_File_write(mpifh, write_buffer.str().c_str(), buffer_pos, MPI_BYTE, MPI_STATUS_IGNORE);

	MPI_File_close(&mpifh);
}
#endif
void DomainDecompBase::writeMoleculesToFile(const std::string& filename, ParticleContainer* moleculeContainer,
                                            bool binary) const {
#ifdef ENABLE_MPI
	if (binary) {
		writeMoleculesToMPIFileBinary(filename, moleculeContainer);
	} else {
#else
		{
#endif
		for (int process = 0; process < getNumProcs(); process++) {
			if (getRank() == process) {
				std::ofstream checkpointfilestream;
				if (binary) {
					auto appendOrTruncate = getRank() != 0 ? std::ios::app : std::ios::trunc;
					// truncate if we are binary and in rank 0, otherwise append!
					checkpointfilestream.open((filename + ".dat").c_str(),
											  std::ios::binary | std::ios::out | appendOrTruncate);
				} else {
					// always append in ascii mode, as we have written the header already to the same file.
					checkpointfilestream.open(filename.c_str(), std::ios::app);
					checkpointfilestream.precision(20);
				}

				for (auto tempMolecule = moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
				     tempMolecule.isValid(); ++tempMolecule) {
					if (binary) {
						tempMolecule->writeBinary(checkpointfilestream);
					} else {
						tempMolecule->write(checkpointfilestream);
					}
				}
				checkpointfilestream.close();
			}
			barrier();
		}
	}
}

void DomainDecompBase::getBoundingBoxMinMax(Domain *domain, double *min, double *max) {
	for(int d = 0; d < 3; d++) {
		min[d] = getBoundingBoxMin(d, domain);
		max[d] = getBoundingBoxMax(d, domain);
	}
}

void DomainDecompBase::collCommInit(int numValues, int /*key*/) {
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



