/**
 * @file GeneralDomainDecomposition.cpp
 * @author seckler, Georg von Bismarck
 * @date 23.09.2026
 */

#include "GeneralDomainDecomposition.h"
#include <csignal>
#include "utils/String_utils.h"
#ifdef ENABLE_ALLLBL
#include "ALLLoadBalancer.h"
#endif

#include "Domain.h"
#include "NeighborAcquirer.h"
#include "NeighbourCommunicationScheme.h"
#include "DomainBox.h"

#include "utils/mardyn_assert.h"
#include "utils/Math.h"

#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <tuple>
#include <sstream>
#include <array>


GeneralDomainDecomposition::GeneralDomainDecomposition(double cutoffRadius, double skin, Domain* domain) : GeneralDomainDecomposition(cutoffRadius, skin, domain, MPI_COMM_WORLD) {}

GeneralDomainDecomposition::GeneralDomainDecomposition(double cutoffRadius, double skin, Domain* domain, MPI_Comm comm) : 
DomainDecompMPIBase(comm),
_cutoffRadius{cutoffRadius},
_skin{skin},
_skinDims{skin,skin,skin},
_domainLength{domain->getGlobalLength(0), domain->getGlobalLength(1), domain->getGlobalLength(2)},
_gridSize({0,0,0}), 
_coords{0, 0,0} {
	initMPIGridDims();
}

void GeneralDomainDecomposition::initMPIGridDims() {
	mardyn_assert(DIMgeom == 3);
	int period[DIMgeom] = {1, 1, 1}; // 1(true) when using periodic boundary conditions in the corresponding dimension
	int reorder = 1; // 1(true) if the ranking may be reordered by MPI_Cart_create
	{
		auto numProcsGridSize = _gridSize[0] * _gridSize[1] * _gridSize[2];
		if (numProcsGridSize != _numProcs and numProcsGridSize != 0) {
			std::ostringstream error_message;
			error_message << "GeneralDomainDecomposition: Wrong grid size given!" << std::endl;
			error_message << "\tnumProcs is " << _numProcs << "," << std::endl;
			error_message << "\tbut grid is " << _gridSize[0] << " x " << _gridSize[1] << " x " << _gridSize[2] << std::endl;
			error_message << "\tresulting in " << numProcsGridSize << " subdomains!" << std::endl;
			error_message << "\tplease check your input file!" << std::endl;
			MARDYN_EXIT(error_message.str());
		}
	}
	MPI_Comm previousComm = _comm;

	MPI_CHECK(MPI_Dims_create( _numProcs, DIMgeom, _gridSize.data()));
	MPI_CHECK(MPI_Cart_create(_comm, DIMgeom, _gridSize.data(), period, reorder, &_comm));

	// If initMPIGridDims has already been executed, the previous communicator is deleted.
	if (_cartCommunicatorCreated) {
		MPI_CHECK(MPI_Comm_free(&previousComm));
	}
	_cartCommunicatorCreated = true;


	Log::global_log->info() << "MPI grid dimensions: " << _gridSize[0] << ", " << _gridSize[1] << ", " << _gridSize[2] << std::endl;
	MPI_CHECK(MPI_Comm_rank(_comm, &_rank));
	MPI_CHECK(MPI_Cart_coords(_comm, _rank, DIMgeom, _coords.data()));
	Log::global_log->info() << "MPI coordinate of current process: " << _coords[0] << ", " << _coords[1] << ", " << _coords[2] << std::endl;

	_localDomain = DomainBox::createAsRegularGrid(_domainLength, _gridSize, _coords);
	_proposedLocalDomain = _localDomain;

	Log::global_log->info() << "GeneralDomainDecomposition initial box: " << _localDomain << std::endl;
}

GeneralDomainDecomposition::~GeneralDomainDecomposition() {
	MPI_CHECK(MPI_Comm_free(&_comm));
}


void GeneralDomainDecomposition::initializeALLLoadBalancer() {
	Log::global_log->info() << "GeneralDomainDecomposition: initializing ALL load balancer..." << std::endl;
#ifdef ENABLE_ALLLBL
	_loadBalancer = std::make_unique<ALLLoadBalancer>(this->getCommunicator(), _gridSize,  _minimalDomainSize);
#else
	std::ostringstream error_message;
	error_message << "ALL load balancing library not enabled. Aborting." << std::endl;
	MARDYN_EXIT(error_message.str());
#endif
}

void GeneralDomainDecomposition::readXML(XMLfileUnits& xmlconfig) {
	// Ensures that the readXML() call to DomainDecompMPIBase forces the direct-pp communication scheme.
	_forceDirectPP = true;

	DomainDecompMPIBase::readXML(xmlconfig);

	#ifdef MARDYN_AUTOPAS
		Log::global_log->info() << "GeneralDomainDecomposition: AutoPas only supports FS, so setting it." << std::endl;
		setCommunicationScheme("direct-pp", "fs");
	#endif

	xmlconfig.getNodeValue("updateFrequency", _rebuildFrequency);
	Log::global_log->info() << "GeneralDomainDecomposition: update frequency: " << _rebuildFrequency << std::endl;

	xmlconfig.getNodeValue("initialPhaseTime", _initPhase);
	Log::global_log->info() << "GeneralDomainDecomposition: time for initial rebalancing phase: " << _initPhase << std::endl;

	xmlconfig.getNodeValue("initialPhaseFrequency", _initFrequency);
	Log::global_log->info() << "GeneralDomainDecomposition: frequency for initial rebalancing phase: " << _initFrequency
					   << std::endl;

	xmlconfig.getNodeValue("imbalanceThresholdCV", _imbalanceThresholdCV);
	xmlconfig.getNodeValue("imbalanceThresholdMinMax", _imbalanceThresholdMinMax);

	if (_imbalanceThresholdCV != 0 && _imbalanceThresholdMinMax != 0) {
		std::ostringstream error_message;
		error_message << "GeneralDomainDecomposition: Multiple imbalanceThresholds were detected, but only one can be supported simultaneously." << std::endl;
		MARDYN_EXIT(error_message.str());
	} else if (_imbalanceThresholdCV != 0) {
		if (_imbalanceThresholdCV < 0) {
			//Coefficient of variation cannot be less than 0
			std::ostringstream error_message;
			error_message << "GeneralDomainDecomposition: imbalanceThresholdCV settion of " << _imbalanceThresholdCV << " is illogical (cannot be less than 0). Aborting! Please select a valid option!";
			MARDYN_EXIT(error_message.str());
		}
		_imbalanceThresholdMode = 1;
		Log::global_log->info() << "GeneralDomainDecomposition: imbalance Threshold Coefficient of variation is active with the value: " << _imbalanceThresholdCV << std::endl;
	} else if (_imbalanceThresholdMinMax != 0) {
		if (_imbalanceThresholdMinMax < 1) {
			// max(data) / min(data) cannot be less than 1
			std::ostringstream error_message;
			error_message << "GeneralDomainDecomposition: imbalanceThresholdMinMax settion of " << _imbalanceThresholdMinMax << " is illogical (cannot be less than 1). Aborting! Please select a valid option!";
			MARDYN_EXIT(error_message.str());
		}
		_imbalanceThresholdMode = 2;
		Log::global_log->info() << "GeneralDomainDecomposition: imbalance Threshold MinMax is active with the value: " << _imbalanceThresholdMinMax << std::endl;
	}

	if (xmlconfig.getNodeValue("maximumRepeatedLoadChange", _maximumRepeatedLoadChange)) {
		if (_maximumRepeatedLoadChange <= 0 || _maximumRepeatedLoadChange > 1) {
			std::ostringstream error_message;
			error_message << "GeneralDomainDecomposition: maximumRepeatedLoadChange settion of " << _maximumRepeatedLoadChange << " is illogical (It is not a percentage). Aborting! Please select a valid option!";
			MARDYN_EXIT(error_message.str());
		}
		Log::global_log->info() << "GeneralDomainDecomposition: maximum Repeated LoadChange is active with the value: " << _maximumRepeatedLoadChange * 100 << "%" << std::endl;
	}

	if(xmlconfig.changecurrentnode("MPIGridDims")) {
		_gridSize[0] = xmlconfig.getNodeValue_int("x", 0);
		_gridSize[1] = xmlconfig.getNodeValue_int("y", 0);
		_gridSize[2] = xmlconfig.getNodeValue_int("z", 0);
		initMPIGridDims();
		xmlconfig.changecurrentnode("..");
	}

	if(xmlconfig.changecurrentnode("skinDims")) {
		_skinDims[0] = xmlconfig.getNodeValue_int("x", 0);
		_skinDims[1] = xmlconfig.getNodeValue_int("y", 0);
		_skinDims[2] = xmlconfig.getNodeValue_int("z", 0);
		if( _skinDims[0] < _skin or _skinDims[1] < _skin or _skinDims[2] < _skin){
				std::ostringstream error_message;
				error_message << "GeneralDomainDecomposition: Manual skin settings must not fall below the boundaries specified by the datastructure. " << std::endl;
				MARDYN_EXIT(error_message.str());
		}
		Log::global_log->info() << "GeneralDomainDecomposition: Using Manual skin values of (" << _skinDims[0] << ", " << _skinDims[1] << ", " << _skinDims[2] << ") for the minimal Domain Size." << std::endl;
		xmlconfig.changecurrentnode("..");
	}

	#ifdef MARDYN_AUTOPAS
		const double minimalDomainBoundary = _cutoffRadius;
	#else
		const double minimalDomainBoundary = 2 * _cutoffRadius;
	#endif

	_minimalDomainSize = {_skinDims[0] + minimalDomainBoundary, _skinDims[1] + minimalDomainBoundary, _skinDims[2] + minimalDomainBoundary};

	// NOTE Currently, forceLatchingToLinkedCellsGrid is not implemented.  
	std::string gridSizeString;
	if (xmlconfig.getNodeValue("gridSize", gridSizeString)) {
		Log::global_log->info() << "GeneralDomainDecomposition: grid size: " << gridSizeString << std::endl;

		if (gridSizeString.find(',') != std::string::npos) {
			auto strings = string_utils::split(gridSizeString, ',');
			if (strings.size() != 3) {
				std::ostringstream error_message;
				error_message
					<< "GeneralDomainDecomposition: gridSize should have three entries if a list is given, but has "
					<< strings.size() << "!" << std::endl;
				MARDYN_EXIT(error_message.str());
			}
			_latchedGridSize = {std::stod(strings[0]), std::stod(strings[1]), std::stod(strings[2])};
		} else {
			const double latchedGridMonoSize = std::stod(gridSizeString);
			_latchedGridSize = {latchedGridMonoSize, latchedGridMonoSize, latchedGridMonoSize};
		}
		_minimalDomainSize = {(*_latchedGridSize)[0], (*_latchedGridSize)[1], (*_latchedGridSize)[2]};
	}

	if(xmlconfig.changecurrentnode("minimalDomainSize")) {
		Log::global_log->info() << "GeneralDomainDecomposition: Manual minimal Domain Size setting found. This is overwriting (skin + minimalDomainBoundary) or latchedGridSize" << std::endl;
		_minimalDomainSize[0] = xmlconfig.getNodeValue_double("x", 0);
		_minimalDomainSize[1] = xmlconfig.getNodeValue_double("y", 0);
		_minimalDomainSize[2] = xmlconfig.getNodeValue_double("z", 0);
		xmlconfig.changecurrentnode("..");
	}
	Log::global_log->info() << "GeneralDomainDecomposition: Using minimal Domain Size of (" << _minimalDomainSize[0] << ", " << _minimalDomainSize[1] << ", " << _minimalDomainSize[2] << ") for the Load Balancer." << std::endl;
	checkMinimalDomainSize(minimalDomainBoundary);

	if (_latchedGridSize.has_value()) {
		// latch the initial RegularGrid
		// IMPORTANT: This must happen after the checkMinimalDomainSize check
		_localDomain.latchToGridSize((*_latchedGridSize), _domainLength);
	}

	if (xmlconfig.changecurrentnode("loadBalancer")) {
		std::string loadBalancerString = "None";
		xmlconfig.getNodeValue("@type", loadBalancerString);
		Log::global_log->info() << "GeneralDomainDecomposition: Chosen Load Balancer: " << loadBalancerString << std::endl;

		std::transform(loadBalancerString.begin(), loadBalancerString.end(), loadBalancerString.begin(), ::tolower);

		if (loadBalancerString.find("all") != std::string::npos) {
			initializeALLLoadBalancer();
		} else {
			std::ostringstream error_message;
			error_message << "GeneralDomainDecomposition: Unknown load balancer " << loadBalancerString
								<< ". Aborting! Please select a valid option! Valid options: ALL";
			MARDYN_EXIT(error_message.str());
		}
		_loadBalancer->readXML(xmlconfig);
		xmlconfig.changecurrentnode("..");
	} else {
		std::ostringstream error_message;
		error_message << "GeneralDomainDecomposition: loadBalancer section missing! Aborting!" << std::endl;
		MARDYN_EXIT(error_message.str());
	}
}

double GeneralDomainDecomposition::getBoundingBoxMin(int dimension, Domain* /*domain*/) { return _localDomain.min[dimension]; }

double GeneralDomainDecomposition::getBoundingBoxMax(int dimension, Domain* /*domain*/) { return _localDomain.max[dimension]; }

bool GeneralDomainDecomposition::checkNeedRebalance(double lastTraversalTime) {
	if (_imbalanceThresholdMode == 0){
		return true; // checkNeedRebalance is disabled
	}
	double globalTraversalTimes[_numProcs];
	MPI_CHECK(MPI_Allgather(&lastTraversalTime, 1, MPI_DOUBLE, globalTraversalTimes, 1, MPI_DOUBLE, _comm)); 
	
	if (_imbalanceThresholdMode == 1) {
		const double value = getCV(globalTraversalTimes, _numProcs);
		Log::global_log->debug() << "GeneralDomainDecomposition: Coefficient of variation: " <<  value << std::endl;
		return value > _imbalanceThresholdCV; 
		
	} else {
		const double value = getMaxdivMin(globalTraversalTimes, _numProcs);
		Log::global_log->debug() << "GeneralDomainDecomposition: Max div Min: " << value << std::endl;
		return value > _imbalanceThresholdMinMax;
	}
}

DomainBox GeneralDomainDecomposition::reviseNewRebalance(DomainBox proposedLocalDomain) noexcept {
	if (_latchedGridSize.has_value()) {
		proposedLocalDomain.latchToGridSize((*_latchedGridSize), _domainLength);
	}
	return proposedLocalDomain;
}

bool GeneralDomainDecomposition::checkForSensibleRebalance(const DomainBox& proposedLocalDomain) {
	if (_maximumRepeatedLoadChange == 1) {
		return true;
	}
	
	if (_previousDomainDecomposition.empty()) {
		const std::size_t numProcs = static_cast<std::size_t>(_numProcs);
		_previousDomainDecompositionChange.reserve(numProcs * (numProcs - 1));
    	_proposedDomainDecompositionChange.reserve(numProcs * (numProcs - 1));
		_mpiBufferDomainDecomposition.resize(6 * numProcs);

		_previousDomainDecomposition.resize(numProcs);
		_currentDomainDecomposition.resize(numProcs);
		_proposedDomainDecomposition.resize(numProcs);

		MPI_CHECK(MPI_Allgather(_localDomain.toArray().data(), 6, MPI_DOUBLE, _mpiBufferDomainDecomposition.data(), 6, MPI_DOUBLE, _comm));
		_previousDomainDecomposition = DomainBox::fromFlatArray(_mpiBufferDomainDecomposition);
		
		MPI_CHECK(MPI_Allgather(proposedLocalDomain.toArray().data(), 6, MPI_DOUBLE, _mpiBufferDomainDecomposition.data(), 6, MPI_DOUBLE, _comm));
		_currentDomainDecomposition = DomainBox::fromFlatArray(_mpiBufferDomainDecomposition);

		return true;
	}
	
	MPI_CHECK(MPI_Allgather(proposedLocalDomain.toArray().data(), 6, MPI_DOUBLE, _mpiBufferDomainDecomposition.data(), 6, MPI_DOUBLE, _comm));
	_proposedDomainDecomposition = DomainBox::fromFlatArray(_mpiBufferDomainDecomposition);

	const double repeatedChange = domainDecompositionPercentageOfRepeatedChanges();
	const bool issensibleRebalance = repeatedChange <= _maximumRepeatedLoadChange;

	Log::global_log->debug() << "GeneralDomainDecomposition: RepeatedLoadChange: " << repeatedChange << std::endl;

	if (issensibleRebalance) {
		//At this point, the load balancing procedure is expected to be carried out.
		_previousDomainDecomposition = _currentDomainDecomposition;
		_currentDomainDecomposition  = _proposedDomainDecomposition;
	}

	return issensibleRebalance;
}


bool GeneralDomainDecomposition::checkRebalancing(size_t step) const noexcept {
	return step <= _initPhase ? step % _initFrequency == 0 : step % _rebuildFrequency == 0;
}

void GeneralDomainDecomposition::balanceAndExchange(double lastTraversalTime, bool forceRebalancing,
													ParticleContainer* moleculeContainer, Domain* domain) {							
	if (_steps == 0) {
		// ensure that there are no outer particles
		moleculeContainer->deleteOuterParticles();
		initCommunicationPartners(domain, moleculeContainer);
		DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);
		++_steps;
		return;
	}

	const bool doRebalance = checkRebalancing(_steps) || forceRebalancing;
	if (doRebalance) {
		const bool needRebalance = checkNeedRebalance(lastTraversalTime);
		if (needRebalance || forceRebalancing) {
			rebalance(lastTraversalTime, moleculeContainer, domain);
		}
		else {
			Log::global_log->info() << "GeneralDomainDecomposition: Skiping rebalancing" << std::endl;
		}
		
	} else {
		if (sendLeavingWithCopies()) {
			Log::global_log->debug() << "GeneralDomainDecomposition: Sending Leaving and Halos." << std::endl;
			DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_AND_HALO_COPIES);
		} else {
			Log::global_log->debug() << "GeneralDomainDecomposition: Sending Leaving." << std::endl;
			DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_ONLY);
			#ifndef MARDYN_AUTOPAS
				moleculeContainer->deleteOuterParticles();
			#endif
			Log::global_log->debug() << "GeneralDomainDecomposition: Sending Halos." << std::endl;
			DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);
		}
	}
	++_steps;		
}

void GeneralDomainDecomposition::initCommunicationPartners(Domain* domain, ParticleContainer* moleculeContainer) { 
	auto coversWholeDomain = _loadBalancer->getCoversWholeDomain();
	for (int d = 0; d < DIMgeom; ++d) {
		_neighbourCommunicationScheme->setCoverWholeDomain(d, coversWholeDomain[d]);
	}
	_neighbourCommunicationScheme->initCommunicationPartners(moleculeContainer->getCutoff(), domain, this,
															 moleculeContainer);
}


void GeneralDomainDecomposition::rebalance(double lastTraversalTime, ParticleContainer* moleculeContainer, Domain* domain) {
	Log::global_log->info() << "GeneralDomainDecomposition: rebalancing..." << std::endl;
	Log::global_log->debug() << "GeneralDomainDecomposition: Sending Leaving." << std::endl;
	DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_ONLY);

	Log::global_log->set_mpi_output_all();
	Log::global_log->debug() << "GeneralDomainDecomposition: work:" << lastTraversalTime << std::endl;
	Log::global_log->set_mpi_output_root(0);
	
	DomainBox proposedLocalDomain = _loadBalancer->rebalance(_proposedLocalDomain, lastTraversalTime);
	DomainBox revisedProposedLocalDomain = reviseNewRebalance(proposedLocalDomain);

	if (!checkForSensibleRebalance(revisedProposedLocalDomain)) {
		Log::global_log->info() << "GeneralDomainDecomposition: rebalancing will be discontinued" << std::endl;

		// Conclude as a normal Particle exchange
		#ifndef MARDYN_AUTOPAS
			moleculeContainer->deleteOuterParticles();
		#endif
		Log::global_log->debug() << "GeneralDomainDecomposition: Sending Halos." << std::endl;
		DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);

		return;
	}

	moleculeContainer->deleteOuterParticles();
																	
	Log::global_log->debug() << "GeneralDomainDecomposition: migrating particles" << std::endl;
	migrateParticles(domain, moleculeContainer, revisedProposedLocalDomain);

	_localDomain = revisedProposedLocalDomain;
	_proposedLocalDomain = proposedLocalDomain;

	#ifndef MARDYN_AUTOPAS
			moleculeContainer->update();
	#endif

	Log::global_log->debug() << "GeneralDomainDecomposition: updating communication partners" << std::endl;
	initCommunicationPartners(domain, moleculeContainer);
	Log::global_log->debug() << "GeneralDomainDecomposition: rebalancing finished" << std::endl;

	Log::global_log->debug() << "GeneralDomainDecomposition: Sending Halos." << std::endl;
	DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);

	_boundaryHandler.setLocalRegion(_localDomain.min.data(),_localDomain.max.data());
	_boundaryHandler.updateGlobalWallLookupTable();
}

void GeneralDomainDecomposition::migrateParticles(Domain* domain, ParticleContainer* moleculeContainer, DomainBox newLocalDomain) {
	HaloRegion ownDomain{}, newDomain{};
	for (size_t i = 0; i < DIMgeom; ++i) {
		ownDomain.rmin[i] = _localDomain.min[i];
		newDomain.rmin[i] = newLocalDomain.min[i];
		ownDomain.rmax[i] = _localDomain.max[i];
		newDomain.rmax[i] = newLocalDomain.max[i];
		ownDomain.offset[i] = 0;
		newDomain.offset[i] = 0;
	}
	Log::global_log->set_mpi_output_all();
	Log::global_log->debug() << "GeneralDomainDecomposition: migrating from " << _localDomain << " to " << newLocalDomain << std::endl;
	Log::global_log->set_mpi_output_root(0);
	std::vector<HaloRegion> desiredDomain{newDomain};
	std::vector<CommunicationPartner> sendNeighbors{}, recvNeighbors{};
	std::vector<Molecule> emigrants;

	std::tie(recvNeighbors, sendNeighbors) =
		NeighborAcquirer::acquireNeighbors(_domainLength, &ownDomain, desiredDomain, _comm);
	if (moleculeContainer->isInvalidParticleReturner()) {
		//AutoPas
		#if false
		{
			//TODO: In rare cases, the code crashes when using Autopass. This can be reproduced by setting the load balancing input to the process rank.
			//Note: This code would be faster than the alternative if it worked without errors
			emigrants = moleculeContainer->rebuildFilter(newMin.data(), newMax.data());
			for (auto& sender : sendNeighbors) {
				sender.initSend(moleculeContainer, _comm, _mpiParticleType, LEAVING_ONLY, emigrants,
								true , false);
			}
		}
		#else
		{
			std::vector<Molecule> dummy;
			for (auto& sender : sendNeighbors) {
				sender.initSend(moleculeContainer, _comm, _mpiParticleType, LEAVING_ONLY, dummy,
								false /*don't use invalid particles*/, false /*do halo position change*/,
								true /*removeFromContainer*/);
			}
			emigrants = moleculeContainer->rebuildFilter(newLocalDomain.min.data(), newLocalDomain.max.data());
			}
		#endif
	} else {
		//LinkedCells
		std::vector<Molecule> dummy;
		for (auto& sender : sendNeighbors) {
			sender.initSend(moleculeContainer, _comm, _mpiParticleType, LEAVING_ONLY, dummy,
							false /*don't use invalid particles*/, false /*do halo position change*/,
							true /*removeFromContainer*/);
		}
		//Note: Changing the domain of the LinkedCells container that contains particles may results in the deletion of particles. 
		std::vector<Molecule> ownMolecules{};
		ownMolecules.reserve(moleculeContainer->getNumberOfParticles());
		for (auto iter = moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); iter.isValid(); ++iter) {
			ownMolecules.push_back(*iter);
		}
		moleculeContainer->clear();
		moleculeContainer->rebuild(newLocalDomain.min.data(), newLocalDomain.max.data());
		moleculeContainer->addParticles(ownMolecules);
	}	

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
			allDone &= recv.testRecv(moleculeContainer, false);
		}

		// catch deadlocks
		double waitingTime = MPI_Wtime() - startTime;
		if (waitingTime > waitCounter) {
			Log::global_log->warning() << "GeneralDomainDecomposition::migrateParticles: Deadlock warning: Rank " << _rank
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
			Log::global_log->error() << "GeneralDomainDecomposition::migrateParticles: Deadlock error: Rank " << _rank
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

	if(not emigrants.empty()){
		std::ostringstream error_message;
		error_message << "GeneralDomainDecomposition: Invalid particles that should have been sent, are still existent. They would be lost. Aborting...\n";						  
		MARDYN_EXIT(error_message.str());
	}

	if (not allDone) {
		std::ostringstream error_message;
		error_message << "A problem occurred during particle migration between old decomposition and new decomposition of the GeneralDomainDecomposition. Aborting." << std::endl;
		MARDYN_EXIT(error_message.str());
	}
	
}

void GeneralDomainDecomposition::checkMinimalDomainSize(double minimalDomainBoundary) {
	for (int i = 0; i < DIMgeom; ++i) {
		constexpr std::array<const char*, 3> axisNames{"X", "Y", "Z"};
		const char* axisName = axisNames[i];
		if (_minimalDomainSize[i] < (minimalDomainBoundary + _skin) or _minimalDomainSize[i] > _localDomain.max[i] - _localDomain.min[i]) {
			std::ostringstream error_message;
			error_message << "GeneralDomainDecomposition: The specified minimal Domain Size is invalid. on the " << axisName << "-axis.Aborting." << std::endl;
			MARDYN_EXIT(error_message.str());
		}
		if (_minimalDomainSize[i] < minimalDomainBoundary + _skinDims[i]) {
			Log::global_log->warning() << "GeneralDomainDecomposition: For the current minimal Domain Size, violates the specified skinDims on the " << axisName << "-axis. " << std::endl;
		}
		if (isNearRel(_minimalDomainSize[i], minimalDomainBoundary)) {
			Log::global_log->warning() << "GeneralDomainDecomposition: For the current minimal Domain Size, the load balancer can reduce the domain to the absolute minimum on the " << axisName << "-axis. " << std::endl;
		}
	}
}

double GeneralDomainDecomposition::getCV(double* data, const int size) const noexcept {
	double sum = 0;
	for( size_t i = 0; i < size; i++ ) {
		sum += data[i];
	}
	const double mean = sum / size;

	double stddev = 0;
	for( size_t i = 0; i < size; i++ ) {
		double diff = data[i] - mean;
		stddev += diff*diff;
	}
	stddev /= size;
	stddev = sqrt(stddev);
	return stddev / mean;
}

double GeneralDomainDecomposition::getMaxdivMin(double* data, const int size) const noexcept {
	double min = data[0];
	double max = data[0];

	for( size_t i = 1; i < size; i++ ) {
		min = std::min(min, data[i]);
		max = std::max(max, data[i]);
	}

	return max / min;
}

double GeneralDomainDecomposition::domainDecompositionPercentageOfRepeatedChanges() noexcept {
    for (int i = 0; i < _numProcs; ++i) {
		for (int j = 0; j < _numProcs; ++j) {
			if (i == j) continue;
			
			std::optional<DomainBox> maybeIntersectionPrevious =  _previousDomainDecomposition[i].intersection(_currentDomainDecomposition[j]);
			if (maybeIntersectionPrevious.has_value()) {
				_previousDomainDecompositionChange.push_back(*maybeIntersectionPrevious);
			}

			std::optional<DomainBox> maybeIntersectionProposed =  _currentDomainDecomposition[i].intersection(_proposedDomainDecomposition[j]);
			if (maybeIntersectionProposed.has_value()) {
				_proposedDomainDecompositionChange.push_back(*maybeIntersectionProposed);
			}
		}
	}

    double total_volume = 0.0;
    double total_intersection = 0.0;

    for (const DomainBox& box : _proposedDomainDecompositionChange) {
        total_volume += box.volume();
    }

    for (const DomainBox& left_box : _previousDomainDecompositionChange) {
        total_volume += left_box.volume();

        for (const DomainBox& right_box : _proposedDomainDecompositionChange) {
			const double inter = left_box.intersectionVolume(right_box);
            total_volume -= inter;
            total_intersection += inter;
        }
    }
	
	_previousDomainDecompositionChange.clear();
	_proposedDomainDecompositionChange.clear();

	if (total_volume == 0.0) {
		Log::global_log->error() << "GeneralDomainDecomposition: maximumRepeatedLoadChange calculated a value that should be impossible. Featcher is disabled for the currently running process." << std::endl;
		return 0.0;
	}
    return total_intersection / total_volume;
}
