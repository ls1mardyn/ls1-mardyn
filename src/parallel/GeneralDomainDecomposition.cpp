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

#include "utils/mardyn_assert.h"

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
		MPI_Comm_free(&previousComm);
	}
	_cartCommunicatorCreated = true;


	Log::global_log->info() << "MPI grid dimensions: " << _gridSize[0] << ", " << _gridSize[1] << ", " << _gridSize[2] << std::endl;
	MPI_CHECK(MPI_Comm_rank(_comm, &_rank));
	MPI_CHECK(MPI_Cart_coords(_comm, _rank, DIMgeom, _coords.data()));
	Log::global_log->info() << "MPI coordinate of current process: " << _coords[0] << ", " << _coords[1] << ", " << _coords[2] << std::endl;

	initializeRegularGrid(_domainLength, _gridSize, _coords);
}

GeneralDomainDecomposition::~GeneralDomainDecomposition() {
	MPI_Comm_free(&_comm);
}


void GeneralDomainDecomposition::initializeALLLoadBalancer() {
	Log::global_log->info() << "GeneralDomainDecomposition: initializing ALL load balancer..." << std::endl;
#ifdef ENABLE_ALLLBL
	_loadBalancer = std::make_unique<ALLLoadBalancer>(4 /*gamma*/, this->getCommunicator(), _gridSize,  _minimalDomainSize);
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

	#ifdef MARDYN_AUTOPAS
		const double minimalDomainBoundary = _cutoffRadius;
	#else
		const double minimalDomainBoundary = 2 * _cutoffRadius;
	#endif

	_minimalDomainSize = {_skin + minimalDomainBoundary, _skin + minimalDomainBoundary, _skin + minimalDomainBoundary};

	// NOTE Currently, forceLatchingToLinkedCellsGrid is not implemented.  
	std::string gridSizeString;
	if (xmlconfig.getNodeValue("gridSize", gridSizeString)) {
		Log::global_log->info() << "GeneralDomainDecomposition grid size: " << gridSizeString << std::endl;

		if (gridSizeString.find(',') != std::string::npos) {
			auto strings = string_utils::split(gridSizeString, ',');
			if (strings.size() != 3) {
				std::ostringstream error_message;
				error_message
					<< "GeneralDomainDecomposition's gridSize should have three entries if a list is given, but has "
					<< strings.size() << "!" << std::endl;
				MARDYN_EXIT(error_message.str());
			}
			_latchedGridSize = {std::stod(strings[0]), std::stod(strings[1]), std::stod(strings[2])};
		} else {
			const double latchedGridMonoSize = std::stod(gridSizeString);
			_latchedGridSize = {latchedGridMonoSize, latchedGridMonoSize, latchedGridMonoSize};
		}
		// latch the initial RegularGrid
		_localDomain = latchToGridSize(_localDomain);

		// TODO Fix Check
		// for (auto gridSize : *_latchedGridSize) {
		// 	if (gridSize < minimalDomainBoundary) {
		// 		std::ostringstream error_message;
		// 		error_message << "GeneralDomainDecomposition's gridSize (" << gridSize
		// 							<< ") is smaller than the interactionLength (" << minimalDomainBoundary
		// 							<< "). This is forbidden, as it leads to errors! " << std::endl;
		// 		MARDYN_EXIT(error_message.str());
		// 	}
		// }
		// _minimalDomainSize = {(*_latchedGridSize)[0], (*_latchedGridSize)[1], (*_latchedGridSize)[2]}; // TODO This line causes problems in the current code. 
	}

	if(xmlconfig.changecurrentnode("minimalDomainSize")) {
		Log::global_log->info() << "GeneralDomainDecomposition: minimalDomainSize setting is overwriting (skin + minimalDomainBoundary) or latchedGridSize" << std::endl;
		_minimalDomainSize[0] = xmlconfig.getNodeValue_double("x", 0);
		_minimalDomainSize[1] = xmlconfig.getNodeValue_double("y", 0);
		_minimalDomainSize[2] = xmlconfig.getNodeValue_double("z", 0);
		xmlconfig.changecurrentnode("..");
	}
	Log::global_log->info() << "GeneralDomainDecomposition: Using minimal Domain Size of (" << _minimalDomainSize[0] << ", " << _minimalDomainSize[1] << ", " << _minimalDomainSize[2] << ") for the Load Balancer." << std::endl;
	checkMinimalDomainSize(minimalDomainBoundary);

	if (xmlconfig.changecurrentnode("loadBalancer")) {
		std::string loadBalancerString = "None";
		xmlconfig.getNodeValue("@type", loadBalancerString);
		Log::global_log->info() << "Chosen Load Balancer: " << loadBalancerString << std::endl;

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
		error_message << "loadBalancer section missing! Aborting!" << std::endl;
		MARDYN_EXIT(error_message.str());
	}
}

double GeneralDomainDecomposition::getBoundingBoxMin(int dimension, Domain* /*domain*/) { return _localDomain[0][dimension]; }

double GeneralDomainDecomposition::getBoundingBoxMax(int dimension, Domain* /*domain*/) { return _localDomain[1][dimension]; }

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

DomainBox GeneralDomainDecomposition::reviseNewRebalance(DomainBox proposedLocalDomain) {
	if (_latchedGridSize.has_value()) {
		proposedLocalDomain = latchToGridSize(proposedLocalDomain);
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
    	_futureDomainDecompositionChange.reserve(numProcs * (numProcs - 1));

		_previousDomainDecomposition.resize(6 * numProcs);
		_currentDomainDecomposition.resize(6 * numProcs);
		_futureDomainDecomposition.resize(6 * numProcs);

		std::array<double, 6> oldDomainBox = {_localDomain[0][0], _localDomain[0][1], _localDomain[0][2], _localDomain[1][0], _localDomain[1][1], _localDomain[1][2]};
		MPI_CHECK(MPI_Allgather(oldDomainBox.data(), 6, MPI_DOUBLE, _previousDomainDecomposition.data(), 6, MPI_DOUBLE, _comm));

		std::array<double, 6> newDomainBox = {proposedLocalDomain[0][0], proposedLocalDomain[0][1], proposedLocalDomain[0][2], proposedLocalDomain[1][0], proposedLocalDomain[1][1], proposedLocalDomain[1][2]};
		MPI_CHECK(MPI_Allgather(newDomainBox.data(), 6, MPI_DOUBLE, _currentDomainDecomposition.data(), 6, MPI_DOUBLE, _comm));

		return true;
	}
	
	std::array<double, 6> newDomainBox = {proposedLocalDomain[0][0], proposedLocalDomain[0][1], proposedLocalDomain[0][2], proposedLocalDomain[1][0], proposedLocalDomain[1][1], proposedLocalDomain[1][2]};
	MPI_Allgather(newDomainBox.data(), 6, MPI_DOUBLE, _futureDomainDecomposition.data(), 6, MPI_DOUBLE, _comm);

	const double repeatedChange = domainDecompositionPercentageOfRepeatedChanges();
	const bool issensibleRebalance = repeatedChange <= _maximumRepeatedLoadChange;

	Log::global_log->debug() << "GeneralDomainDecomposition: RepeatedLoadChange: " << repeatedChange << std::endl;

	if (issensibleRebalance) {
		//At this point, the load balancing procedure is expected to be carried out.
		_previousDomainDecomposition = _currentDomainDecomposition;
		_currentDomainDecomposition  = _futureDomainDecomposition;
	}

	return issensibleRebalance;
}


bool GeneralDomainDecomposition::checkRebalancing(size_t step) {
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

	_boundaryHandler.setLocalRegion(_localDomain[0].data(),_localDomain[1].data());
	_boundaryHandler.updateGlobalWallLookupTable();
}

void GeneralDomainDecomposition::migrateParticles(Domain* domain, ParticleContainer* particleContainer, DomainBox newLocalDomain) {
	HaloRegion ownDomain{}, newDomain{};
	for (size_t i = 0; i < DIMgeom; ++i) {
		ownDomain.rmin[i] = _localDomain[0][i];
		newDomain.rmin[i] = newLocalDomain[0][i];
		ownDomain.rmax[i] = _localDomain[1][i];
		newDomain.rmax[i] = newLocalDomain[1][i];
		ownDomain.offset[i] = 0;
		newDomain.offset[i] = 0;
	}
	Log::global_log->set_mpi_output_all();
	Log::global_log->debug() << "GeneralDomainDecomposition: migrating from"
						<< " [" << _localDomain[0][0] << ", " << _localDomain[1][0] << "] x"
						<< " [" << _localDomain[0][1] << ", " << _localDomain[1][1] << "] x"
						<< " [" << _localDomain[0][2] << ", " << _localDomain[1][2] << "] " << std::endl;
	Log::global_log->debug() << "GeneralDomainDecomposition: to"
						<< " [" << newLocalDomain[0][0] << ", " << newLocalDomain[1][0] << "] x"
						<< " [" << newLocalDomain[0][1] << ", " << newLocalDomain[1][1] << "] x"
						<< " [" << newLocalDomain[0][2] << ", " << newLocalDomain[1][2] << "]." << std::endl;
	Log::global_log->set_mpi_output_root(0);
	std::vector<HaloRegion> desiredDomain{newDomain};
	std::vector<CommunicationPartner> sendNeighbors{}, recvNeighbors{};
	std::vector<Molecule> emigrants;

	std::tie(recvNeighbors, sendNeighbors) =
		NeighborAcquirer::acquireNeighbors(_domainLength, &ownDomain, desiredDomain, _comm);
	if (particleContainer->isInvalidParticleReturner()) {
		//AutoPas
		#if false
		{
			//TODO: In rare cases, the code crashes when using Autopass. This can be reproduced by setting the load balancing input to the process rank. 
			emigrants = particleContainer->rebuildFilter(newMin.data(), newMax.data());
			for (auto& sender : sendNeighbors) {
				sender.initSend(particleContainer, _comm, _mpiParticleType, LEAVING_ONLY, emigrants,
								true , false);
			}
		}
		#else
		{
			std::vector<Molecule> dummy;
			for (auto& sender : sendNeighbors) {
				sender.initSend(particleContainer, _comm, _mpiParticleType, LEAVING_ONLY, dummy,
								false /*don't use invalid particles*/, false /*do halo position change*/,
								true /*removeFromContainer*/);
			}
			emigrants = particleContainer->rebuildFilter(newLocalDomain[0].data(), newLocalDomain[1].data());
			}
		#endif
	} else {
		//LinkedCells
		std::vector<Molecule> dummy;
		for (auto& sender : sendNeighbors) {
			sender.initSend(particleContainer, _comm, _mpiParticleType, LEAVING_ONLY, dummy,
							false /*don't use invalid particles*/, false /*do halo position change*/,
							true /*removeFromContainer*/);
		}
		//Note: Changing the domain of the LinkedCells container that contains particles may results in the deletion of particles. 
		std::vector<Molecule> ownMolecules{};
		ownMolecules.reserve(particleContainer->getNumberOfParticles());
		for (auto iter = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); iter.isValid(); ++iter) {
			ownMolecules.push_back(*iter);
		}
		particleContainer->clear();
		particleContainer->rebuild(newLocalDomain[0].data(), newLocalDomain[1].data());
		particleContainer->addParticles(ownMolecules);
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
			allDone &= recv.testRecv(particleContainer, false);
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

void GeneralDomainDecomposition::initializeRegularGrid(const DomainPoint& domainLength, const DomainGridPoint& gridSize, const DomainGridPoint& gridCoords) {
	_localDomain = {{{0., 0., 0.}, {0., 0., 0.}}}; 

	// initialize it as regular grid!
	for (int dim = 0; dim < DIMgeom; ++dim) {
		_localDomain[0][dim] = gridCoords[dim] * domainLength[dim] / gridSize[dim];
		_localDomain[1][dim] = (gridCoords[dim] + 1) * domainLength[dim] / gridSize[dim];
		if (gridCoords[dim] == gridSize[dim] - 1) {
			// ensure that the upper domain boundaries match.
			// lower domain boundaries always match, because they are 0.
			_localDomain[1][dim] = domainLength[dim];
		}
	}
	_proposedLocalDomain = _localDomain;

	Log::global_log->info() << "GeneralDomainDecomposition initial box: [" << _localDomain[0][0] << ", " << _localDomain[1][0] << "] x ["
			<< _localDomain[0][1] << ", " << _localDomain[1][1] << "] x [" << _localDomain[0][2] << ", " << _localDomain[1][2] << "]"
			<< std::endl;
}

void GeneralDomainDecomposition::checkMinimalDomainSize(double minimalDomainBoundary) {
	for (int i = 0; i < DIMgeom; ++i) {
		if (_minimalDomainSize[i] < minimalDomainBoundary or _minimalDomainSize[i] > _localDomain[1][i] - _localDomain[0][i]) {
			std::ostringstream error_message;
			error_message << "GeneralDomainDecomposition: The specified minimal DomainSize is invalid. Aborting." << std::endl;
			MARDYN_EXIT(error_message.str());
		}
	}
}

double GeneralDomainDecomposition::getCV(double* data, const int size) {
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

double GeneralDomainDecomposition::getMaxdivMin(double* data, const int size) {
	double min = data[0];
	double max = data[0];

	for( size_t i = 1; i < size; i++ ) {
		min = std::min(min, data[i]);
		max = std::max(max, data[i]);
	}

	return max / min;
}

DomainBox GeneralDomainDecomposition::latchToGridSize(DomainBox proposedLocalDomain) const {
	for (size_t ind = 0; ind < DIMgeom; ++ind) {
		double currentGridSize = (*_latchedGridSize)[ind];
		// For boxmin, the lower domain boundary is 0, so that's always fine!
		proposedLocalDomain[0][ind] = std::round(proposedLocalDomain[0][ind] / currentGridSize) * currentGridSize;
		// update boxmax only if it isn't at the very top of the domain!
		if (proposedLocalDomain[1][ind] != _domainLength[ind]) {
			proposedLocalDomain[1][ind] = std::round(proposedLocalDomain[1][ind] / currentGridSize) * currentGridSize;
		}
	}
	return proposedLocalDomain;
}

inline double GeneralDomainDecomposition::bboxVolume(const DomainBox& bbox) const {
    return std::max(bbox[1][0] - bbox[0][0], 0.0) *
           std::max(bbox[1][1] - bbox[0][1], 0.0) *
           std::max(bbox[1][2] - bbox[0][2], 0.0);
}

inline double GeneralDomainDecomposition::bboxIntersectionVolume(const DomainBox& bbox1, const DomainBox& bbox2) const {
    const DomainPoint lower{
        std::max(bbox1[0][0], bbox2[0][0]),
        std::max(bbox1[0][1], bbox2[0][1]),
        std::max(bbox1[0][2], bbox2[0][2])
    };

    const DomainPoint upper{
        std::min(bbox1[1][0], bbox2[1][0]),
        std::min(bbox1[1][1], bbox2[1][1]),
        std::min(bbox1[1][2], bbox2[1][2])
    };

    return bboxVolume(DomainBox{lower, upper});
}

double GeneralDomainDecomposition::domainDecompositionPercentageOfRepeatedChanges() {
    for (int i = 0; i < _numProcs; ++i) {
		for (int j = 0; j < _numProcs; ++j) {
			if (i == j) continue;

			const int ci = 6 * i;
			const int lj = 6 * j;

			const DomainPoint previousLower{
				std::max(_currentDomainDecomposition[ci + 0], _previousDomainDecomposition[lj + 0]),
				std::max(_currentDomainDecomposition[ci + 1], _previousDomainDecomposition[lj + 1]),
				std::max(_currentDomainDecomposition[ci + 2], _previousDomainDecomposition[lj + 2])
			};

			const DomainPoint previousUpper{
				std::min(_currentDomainDecomposition[ci + 3], _previousDomainDecomposition[lj + 3]),
				std::min(_currentDomainDecomposition[ci + 4], _previousDomainDecomposition[lj + 4]),
				std::min(_currentDomainDecomposition[ci + 5], _previousDomainDecomposition[lj + 5])
			};

			if (previousLower[0] < previousUpper[0] &&
				previousLower[1] < previousUpper[1] &&
				previousLower[2] < previousUpper[2]) {
				_previousDomainDecompositionChange.push_back(DomainBox{previousLower, previousUpper});
			}

			const DomainPoint futureLower{
				std::max(_currentDomainDecomposition[ci + 0], _futureDomainDecomposition[lj + 0]),
				std::max(_currentDomainDecomposition[ci + 1], _futureDomainDecomposition[lj + 1]),
				std::max(_currentDomainDecomposition[ci + 2], _futureDomainDecomposition[lj + 2])
			};

			const DomainPoint futureUpper{
				std::min(_currentDomainDecomposition[ci + 3], _futureDomainDecomposition[lj + 3]),
				std::min(_currentDomainDecomposition[ci + 4], _futureDomainDecomposition[lj + 4]),
				std::min(_currentDomainDecomposition[ci + 5], _futureDomainDecomposition[lj + 5])
			};

			if (futureLower[0] < futureUpper[0] &&
				futureLower[1] < futureUpper[1] &&
				futureLower[2] < futureUpper[2]) {
				_futureDomainDecompositionChange.push_back(DomainBox{futureLower, futureUpper});
			}
		}
	}

    double total_volume = 0.0;
    double total_intersection = 0.0;

    for (const DomainBox& box : _futureDomainDecompositionChange) {
        total_volume += bboxVolume(box);
    }

    for (const DomainBox& left_box : _previousDomainDecompositionChange) {
        total_volume += bboxVolume(left_box);

        for (const DomainBox& right_box : _futureDomainDecompositionChange) {
            const double inter = bboxIntersectionVolume(left_box, right_box);
			mardyn_assert(inter < 0);
            total_volume -= inter;
            total_intersection += inter;
        }
    }
	_previousDomainDecompositionChange.clear();
	_futureDomainDecompositionChange.clear();

    return total_intersection / total_volume;
}
