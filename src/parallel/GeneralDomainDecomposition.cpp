/**
 * @file GeneralDomainDecomposition.cpp
 * @author seckler
 * @date 11.04.19
 */

#include "GeneralDomainDecomposition.h"
#ifdef ENABLE_ALLLBL
#include "ALLLoadBalancer.h"
#endif
#include "Domain.h"
#include "NeighborAcquirer.h"
#include "NeighbourCommunicationScheme.h"

GeneralDomainDecomposition::GeneralDomainDecomposition(double interactionLength, Domain* domain)
	: _boxMin{0.}, _boxMax{0.} {
	std::array<double, 3> domainLength = {domain->getGlobalLength(0), domain->getGlobalLength(1),
										  domain->getGlobalLength(2)};

	auto gridSize = getOptimalGrid(domainLength, this->getNumProcs());
	auto gridCoords = getCoordsFromRank(gridSize, _rank);
	_coversWholeDomain = {gridSize[0] == 1, gridSize[1] == 1, gridSize[2] == 1};
	global_log->info() << "gridSize:" << gridSize[0] << ", " << gridSize[1] << ", " << gridSize[2] << std::endl;
	global_log->info() << "gridCoords:" << gridCoords[0] << ", " << gridCoords[1] << ", " << gridCoords[2] << std::endl;
	std::tie(_boxMin, _boxMax) = initializeRegularGrid(domainLength, gridSize, gridCoords);
#ifdef ENABLE_ALLLBL
	// Increased slightly to prevent rounding errors.
	double minimalDomainSize{interactionLength * (1. + 1.e-10)};
	_loadBalancer = std::make_unique<ALLLoadBalancer>(_boxMin, _boxMax, 4 /*gamma*/, this->getCommunicator(), gridSize,
													  gridCoords, minimalDomainSize);
#else
	global_log->error() << "ALL load balancing library not enabled. Aborting." << std::endl;
	Simulation::exit(24235);
#endif

	global_log->info() << "GeneralDomainDecomposition initial box: [" << _boxMin[0] << ", " << _boxMax[0] << "] x ["
					   << _boxMin[1] << ", " << _boxMax[1] << "] x [" << _boxMin[2] << ", " << _boxMax[2] << "]"
					   << std::endl;
}

double GeneralDomainDecomposition::getBoundingBoxMin(int dimension, Domain* /*domain*/) { return _boxMin[dimension]; }

GeneralDomainDecomposition::~GeneralDomainDecomposition() = default;

double GeneralDomainDecomposition::getBoundingBoxMax(int dimension, Domain* /*domain*/) { return _boxMax[dimension]; }

bool GeneralDomainDecomposition::queryRebalancing(size_t step, size_t updateFrequency, size_t initPhase,
												  size_t initUpdateFrequency, double /*lastTraversalTime*/) {
	return step <= initPhase ? step % initUpdateFrequency == 0 : step % updateFrequency == 0;
}

void GeneralDomainDecomposition::balanceAndExchange(double lastTraversalTime, bool forceRebalancing,
													ParticleContainer* moleculeContainer, Domain* domain) {
	bool rebalance =
		queryRebalancing(_steps, _rebuildFrequency, _initPhase, _initFrequency, lastTraversalTime) or forceRebalancing;
	if (_steps == 0) {
		// ensure that there are no outer particles
		moleculeContainer->deleteOuterParticles();
		// init communication partners
		initCommPartners(moleculeContainer, domain);
		DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);
	} else {
		if (rebalance) {
			if (moleculeContainer->isInvalidParticleReturner() and not moleculeContainer->hasInvalidParticles()) {
				moleculeContainer->forcedUpdate();
			}
			// first transfer leaving particles
			DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_ONLY);

			// ensure that there are no outer particles
			moleculeContainer->deleteOuterParticles();

			// rebalance
			std::array<double, 3> newBoxMin{0.}, newBoxMax{0.};
			global_log->info() << "rebalancing..." << std::endl;

			global_log->set_mpi_output_all();
			global_log->debug() << "work:" << lastTraversalTime << std::endl;
			global_log->set_mpi_output_root(0);
			std::tie(newBoxMin, newBoxMax) = _loadBalancer->rebalance(lastTraversalTime);

			// migrate the particles, this will rebuild the moleculeContainer!
			global_log->info() << "migrating particles" << std::endl;
			migrateParticles(domain, moleculeContainer, newBoxMin, newBoxMax);

			// set new boxMin and boxMax
			_boxMin = newBoxMin;
			_boxMax = newBoxMax;

			// init communication partners
			global_log->info() << "updating communication partners" << std::endl;
			initCommPartners(moleculeContainer, domain);
			global_log->info() << "rebalancing finished" << std::endl;
			DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);
		} else {
			if (not moleculeContainer->isInvalidParticleReturner() or moleculeContainer->hasInvalidParticles()) {
				if (sendLeavingWithCopies()) {
					DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_AND_HALO_COPIES);
				} else {
					DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_ONLY);
					moleculeContainer->deleteOuterParticles();
					DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);
				}
			} else {
				DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES,
														  false /*dohaloPositionCheck*/);
			}
		}
	}
	++_steps;
}

void GeneralDomainDecomposition::migrateParticles(Domain* domain, ParticleContainer* particleContainer,
												  array<double, 3> newMin, array<double, 3> newMax) {
	std::array<double, 3> oldBoxMin{particleContainer->getBoundingBoxMin(0), particleContainer->getBoundingBoxMin(1),
									particleContainer->getBoundingBoxMin(2)};
	std::array<double, 3> oldBoxMax{particleContainer->getBoundingBoxMax(0), particleContainer->getBoundingBoxMax(1),
									particleContainer->getBoundingBoxMax(2)};

	HaloRegion ownDomain{}, newDomain{};
	for (size_t i = 0; i < 3; ++i) {
		ownDomain.rmin[i] = oldBoxMin[i];
		newDomain.rmin[i] = newMin[i];
		ownDomain.rmax[i] = oldBoxMax[i];
		newDomain.rmax[i] = newMax[i];
		ownDomain.offset[i] = 0;
		newDomain.offset[i] = 0;
	}
	global_log->set_mpi_output_all();
	global_log->debug() << "migrating from"
						<< " [" << oldBoxMin[0] << ", " << oldBoxMax[0] << "] x"
						<< " [" << oldBoxMin[1] << ", " << oldBoxMax[1] << "] x"
						<< " [" << oldBoxMin[2] << ", " << oldBoxMax[2] << "] " << std::endl;
	global_log->debug() << "to"
						<< " [" << newMin[0] << ", " << newMax[0] << "] x"
						<< " [" << newMin[1] << ", " << newMax[1] << "] x"
						<< " [" << newMin[2] << ", " << newMax[2] << "]." << std::endl;
	global_log->set_mpi_output_root(0);
	std::vector<HaloRegion> desiredDomain{newDomain};
	std::vector<CommunicationPartner> sendNeighbors{}, recvNeighbors{};

	std::array<double, 3> globalDomainLength{domain->getGlobalLength(0), domain->getGlobalLength(1),
											 domain->getGlobalLength(2)};
	// 0. skin, as it is not needed for the migration of particles!
	std::tie(recvNeighbors, sendNeighbors) =
		NeighborAcquirer::acquireNeighbors(globalDomainLength, &ownDomain, desiredDomain, 0. /*skin*/, _comm);

	std::vector<Molecule> dummy;
	for (auto& sender : sendNeighbors) {
		sender.initSend(particleContainer, _comm, _mpiParticleType, LEAVING_ONLY, dummy,
						false /*don't use invalid particles*/, true /*do halo position change*/,
						true /*removeFromContainer*/);
	}
	// TODO: copying own molecules out and reinserting them can be done within autopas more efficiently
	std::vector<Molecule> ownMolecules{};
	for (auto iter = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); iter.isValid(); ++iter) {
		ownMolecules.push_back(*iter);
		// TODO: This check should be in debug mode only
		if (not iter->inBox(newMin.data(), newMax.data())) {
			global_log->error_always_output()
				<< "particle still in domain that should have been migrated." << std::endl;
			Simulation::exit(2315);
		}
	}
	particleContainer->clear();
	particleContainer->rebuild(newMin.data(), newMax.data());
	particleContainer->addParticles(ownMolecules);
	bool allDone = false;
	double waitCounter = 30.0;
	double deadlockTimeOut = 360.0;
	double startTime = MPI_Wtime();
	while (not allDone) {
		allDone = true;

		// "kickstart" processing of all Isend requests
		for (auto& sender : sendNeighbors) {
			allDone &= sender.testSend();
		}

		// unpack molecules
		for (auto& recv : recvNeighbors) {
			allDone &= recv.iprobeCount(this->getCommunicator(), this->getMPIParticleType());
			allDone &= recv.testRecv(particleContainer, false);
		}

		// catch deadlocks
		double waitingTime = MPI_Wtime() - startTime;
		if (waitingTime > waitCounter) {
			global_log->warning() << "KDDecomposition::migrateParticles: Deadlock warning: Rank " << _rank
								  << " is waiting for more than " << waitCounter << " seconds" << std::endl;
			waitCounter += 1.0;
			for (auto& sender : sendNeighbors) {
				sender.deadlockDiagnosticSend();
			}
			for (auto& recv : recvNeighbors) {
				recv.deadlockDiagnosticRecv();
			}
		}

		if (waitingTime > deadlockTimeOut) {
			global_log->error() << "KDDecomposition::migrateParticles: Deadlock error: Rank " << _rank
								<< " is waiting for more than " << deadlockTimeOut << " seconds" << std::endl;
			for (auto& sender : sendNeighbors) {
				sender.deadlockDiagnosticSend();
			}
			for (auto& recv : recvNeighbors) {
				recv.deadlockDiagnosticRecv();
			}
			break;
		}
	}
}

void GeneralDomainDecomposition::initCommPartners(ParticleContainer* moleculeContainer,
												  Domain* domain) {  // init communication partners
	for (int d = 0; d < DIMgeom; ++d) {
		// this needs to be updated for proper initialization of the neighbours
		_neighbourCommunicationScheme->setCoverWholeDomain(d, _coversWholeDomain[d]);
	}
	_neighbourCommunicationScheme->initCommunicationPartners(moleculeContainer->getCutoff(), domain, this,
															 moleculeContainer);
}

void GeneralDomainDecomposition::readXML(XMLfileUnits& xmlconfig) {
	DomainDecompMPIBase::readXML(xmlconfig);
	global_log->info()
		<< "The GeneralDomainDecomposition is enforcing the direct-pp neighbor scheme using fs, so setting it."
		<< std::endl;
	setCommunicationScheme("direct-pp", "fs");

	xmlconfig.getNodeValue("updateFrequency", _rebuildFrequency);
	global_log->info() << "GeneralDomainDecomposition update frequency: " << _rebuildFrequency << endl;

	xmlconfig.getNodeValue("initialPhaseTime", _initPhase);
	global_log->info() << "GeneralDomainDecomposition time for initial rebalancing phase: " << _initPhase << endl;

	xmlconfig.getNodeValue("initialPhaseFrequency", _initFrequency);
	global_log->info() << "GeneralDomainDecomposition frequency for initial rebalancing phase: " << _initFrequency
					   << endl;
}

/**
 * Get the ordering of the input data.
 * The ordering will contain indices of elements of data, starting with the smallest going to the biggest.
 * e.g.:
 * returns {2, 0, 3, 1} for data={5, 16, 4, 7}
 * @tparam ArrayType should be some array type, e.g., vector of int
 * @param data the data to define the order
 * @return The ordering.
 */
template <typename ArrayType>
std::vector<size_t> getOrdering(const ArrayType& data) {
	std::vector<size_t> index(data.size(), 0);
	for (int i = 0; i != index.size(); i++) {
		index[i] = i;
	}
	std::sort(index.begin(), index.end(), [&](const size_t& a, const size_t& b) { return (data[a] < data[b]); });
	return index;
}

std::array<size_t, 3> GeneralDomainDecomposition::getOptimalGrid(const std::array<double, 3>& domainLength,
																 int numProcs) {
	// generate default grid
	std::array<int, 3> gridSize{0};
	MPI_CHECK(MPI_Dims_create(numProcs, 3, gridSize.data()));
	std::sort(gridSize.begin(), gridSize.end());

	// we want the grid to actually resemble the domain a bit (long sides -> many processes!)
	std::array<size_t, 3> grid{0};
	auto ordering = getOrdering(domainLength);
	for (size_t i = 0; i < 3; ++i) {
		grid[ordering[i]] = gridSize[i];
	}
	return grid;
}

std::array<size_t, 3> GeneralDomainDecomposition::getCoordsFromRank(const std::array<size_t, 3>& gridSize, int rank) {
	auto yzSize = gridSize[1] * gridSize[2];
	auto zSize = gridSize[2];
	auto x = rank / yzSize;
	auto y = (rank - x * yzSize) / zSize;
	auto z = (rank - x * yzSize - y * zSize);
	return {x, y, z};
}

std::tuple<std::array<double, 3>, std::array<double, 3>> GeneralDomainDecomposition::initializeRegularGrid(
	const std::array<double, 3>& domainLength, const std::array<size_t, 3>& gridSize,
	const std::array<size_t, 3>& gridCoords) {
	std::array<double, 3> boxMin{0.};
	std::array<double, 3> boxMax{0.};
	// initialize it as regular grid!
	for (size_t dim = 0; dim < 3; ++dim) {
		boxMin[dim] = gridCoords[dim] * domainLength[dim] / gridSize[dim];
		boxMax[dim] = (gridCoords[dim] + 1) * domainLength[dim] / gridSize[dim];
		if (gridCoords[dim] == gridSize[dim] - 1) {
			// ensure that the upper domain boundaries match.
			// lower domain boundaries always match, because they are 0.
			boxMax[dim] = domainLength[dim];
		}
	}
	return std::make_tuple(boxMin, boxMax);
}

void GeneralDomainDecomposition::printDecomp(const std::string &filename,
                                             Domain *domain,
                                             ParticleContainer *particleContainer) {
	if (_rank == 0) {
		ofstream povcfgstrm(filename.c_str());
		povcfgstrm << "size " << domain->getGlobalLength(0) << " " << domain->getGlobalLength(1) << " "
				   << domain->getGlobalLength(2) << endl;
		povcfgstrm << "decompData Regions Configuration" << endl;
		povcfgstrm.close();
	}

	stringstream output;
	output << getBoundingBoxMin(0, domain) << " " << getBoundingBoxMin(1, domain) << " " << getBoundingBoxMin(2, domain)
		   << " " << getBoundingBoxMax(0, domain) << " " << getBoundingBoxMax(1, domain) << " "
		   << getBoundingBoxMax(2, domain) << " " << particleContainer->getConfigurationAsString() << "\n";
	string output_str = output.str();

#ifdef ENABLE_MPI
	MPI_File fh;
	MPI_File_open(_comm, filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_APPEND | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	uint64_t write_size = output_str.size();
	uint64_t offset = 0;
	if (_rank == 0) {
		MPI_Offset file_end_pos;
		MPI_File_seek(fh, 0, MPI_SEEK_END);
		MPI_File_get_position(fh, &file_end_pos);
		write_size += file_end_pos;
		MPI_Exscan(&write_size, &offset, 1, MPI_UINT64_T, MPI_SUM, _comm);
		offset += file_end_pos;
	} else {
		MPI_Exscan(&write_size, &offset, 1, MPI_UINT64_T, MPI_SUM, _comm);
	}
	MPI_File_write_at(fh, offset, output_str.c_str(), output_str.size(), MPI_CHAR, MPI_STATUS_IGNORE);
	MPI_File_close(&fh);
#else
	ofstream povcfgstrm(filename.c_str(), ios::app);
	povcfgstrm << output_str;
	povcfgstrm.close();
#endif
}
