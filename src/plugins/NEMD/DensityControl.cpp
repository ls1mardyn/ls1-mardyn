/*
 * DensityControl.cpp
 *
 *  Created on: 16.11.2021
 *  Author: mheinen
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>
#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "DensityControl.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/CommVar.h"
#include "utils/Logger.h"

DensityControl::DensityControl() = default;
DensityControl::~DensityControl() = default;

void DensityControl::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
	Log::global_log->debug() << "DensityControl enabled" << std::endl;

#ifdef ENABLE_MPI
	// Create new MPI data type to transport particle ID and component ID together in one object
	const int nitems = 2;
	const int blocklengths[2] = {1, 1};
	MPI_Datatype types[2] = {MPI_UINT64_T, MPI_UINT32_T};
	MPI_Aint offsets[2];
	offsets[0] = offsetof(pacIDtype, pid);
	offsets[1] = offsetof(pacIDtype, cid);

	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &pacID_mpi_type);
	MPI_Type_commit(&pacID_mpi_type);
#endif
}

void DensityControl::readXML(XMLfileUnits& xmlconfig) {
	// Timestep control
	_control.start = 0;
	_control.freq = 1;
	_control.stop = 10000;
	xmlconfig.getNodeValue("control/start", _control.start);
	xmlconfig.getNodeValue("control/frequency", _control.freq);
	xmlconfig.getNodeValue("control/stop", _control.stop);
	Log::global_log->info() << "[DensityControl] is acting start:freq:stop = " << _control.start << ":" << _control.freq
					   << ":" << _control.stop << std::endl;

	// range
	Domain* domain = global_simulation->getDomain();
	_range.xmin = 0.;
	_range.xmax = domain->getGlobalLength(0);
	_range.ymin = 0.;
	_range.ymax = domain->getGlobalLength(1);
	_range.zmin = 0.;
	_range.zmax = domain->getGlobalLength(2);
	xmlconfig.getNodeValue("range/xmin", _range.xmin);
	xmlconfig.getNodeValue("range/xmax", _range.xmax);
	xmlconfig.getNodeValue("range/ymin", _range.ymin);
	xmlconfig.getNodeValue("range/ymax", _range.ymax);
	xmlconfig.getNodeValue("range/zmin", _range.zmin);
	xmlconfig.getNodeValue("range/zmax", _range.zmax);
	_range.volume = (_range.xmax - _range.xmin) * (_range.ymax - _range.ymin) * (_range.zmax - _range.zmin);

	// priority (Usually: Molecule size sorted in descending order)
	const uint32_t numComponents = global_simulation->getEnsemble()->getComponents()->size();
	std::string strPrio = "";
	xmlconfig.getNodeValue("priority", strPrio);
	_vecPriority.push_back(0);
	const uint32_t nRet = this->tokenize_int_list(_vecPriority, strPrio);
	if (nRet != numComponents) {
		Log::global_log->error() << "[DensityControl] Number of component IDs specified in element <priority>...</priority>"
							<< " does not match the number of components in the simulation. Programm exit ..."
							<< std::endl;
		Simulation::exit(-1);
	}

	// targets
	{  // init target values
		_densityTarget.density.resize(numComponents + 1);
		_densityTarget.count.resize(numComponents + 1);
		_densityTarget.specified.resize(numComponents + 1);
		for (auto it : _densityTarget.density) {
			it = 0.;
		}
		for (auto it : _densityTarget.count) {
			it = 0.;
		}
		for (auto it : _densityTarget.specified) {
			it = false;
		}
	}

	uint32_t numTargets = 0;
	const XMLfile::Query query = xmlconfig.query("targets/target");
	numTargets = query.card();
	Log::global_log->info() << "[DensityControl] Number of component targets: " << numTargets << std::endl;
	if (numTargets < 1) {
		Log::global_log->error() << "[DensityControl] No target parameters specified. Program exit ..." << std::endl;
		Simulation::exit(-1);
	}
	const std::string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query::const_iterator nodeIter;
	for (nodeIter = query.begin(); nodeIter; nodeIter++) {
		xmlconfig.changecurrentnode(nodeIter);
		uint32_t cid_ub;
		xmlconfig.getNodeValue("@cid", cid_ub);
		double density = 0.;
		xmlconfig.getNodeValue("density", density);
		_densityTarget.density[cid_ub] = density;
		const uint64_t count = round(_range.volume * density);
		_densityTarget.count[cid_ub] = count;
		_densityTarget.specified.at(cid_ub) = true;
	}  // loop over 'targets/target' nodes
	xmlconfig.changecurrentnode(oldpath);
}

void DensityControl::beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
								  unsigned long simstep) {
	if (simstep == _control.start) {
		this->initTargetValues(particleContainer, domainDecomp);
	}

	if (simstep < _control.start || simstep > _control.stop || simstep % _control.freq != 0) return;
	this->controlDensity(particleContainer, domainDecomp, simstep);
}

void DensityControl::controlDensity(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
									unsigned long simstep) {
	Domain* domain = global_simulation->getDomain();
	domain->updateMaxMoleculeID(particleContainer, domainDecomp);
	CommVar<uint64_t> maxMoleculeID = domain->getMaxMoleculeID();
	const int numProcs = domainDecomp->getNumProcs();
	uint32_t numComponents = global_simulation->getEnsemble()->getComponents()->size();
	double regionLowCorner[3], regionHighCorner[3];

	// if linked cell in the region of interest
	for (unsigned d = 0; d < 3; d++) {
		regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
		regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
	}

	// ensure that we do not iterate over things outside of the container.
	regionLowCorner[0] = std::max(_range.xmin, regionLowCorner[0]);
	regionLowCorner[1] = std::max(_range.ymin, regionLowCorner[1]);
	regionLowCorner[2] = std::max(_range.zmin, regionLowCorner[2]);

	regionHighCorner[0] = std::min(_range.xmax, regionHighCorner[0]);
	regionHighCorner[1] = std::min(_range.ymax, regionHighCorner[1]);
	regionHighCorner[2] = std::min(_range.zmax, regionHighCorner[2]);

	CommVar<uint64_t> numMolecules;
	numMolecules.local = numMolecules.global = 0;
	std::array<double, 3> pos;

	CommVar<std::vector<pacIDtype> > vec_pacID;

	auto begin = particleContainer->regionIterator(regionLowCorner, regionHighCorner,
												   ParticleIterator::ONLY_INNER_AND_BOUNDARY);  // over all cell types
	for (auto it = begin; it.isValid(); ++it) {
		// check position
		pos[0] = it->r(0);
		pos[1] = it->r(1);
		pos[2] = it->r(2);
		if (not this->moleculeInsideRange(pos)) continue;

		const uint64_t pid = it->getID();
		const uint32_t cid_ub = it->componentid() + 1;
		pacIDtype pacID;
		pacID.pid = pid;
		pacID.cid = cid_ub;

		vec_pacID.local.push_back(pacID);
	}
	numMolecules.local = vec_pacID.local.size();

#ifdef ENABLE_MPI
	domainDecomp->collCommInit(1);
	domainDecomp->collCommAppendUnsLong(numMolecules.local);
	domainDecomp->collCommAllreduceSum();
	numMolecules.global = domainDecomp->collCommGetUnsLong();
	domainDecomp->collCommFinalize();
#else
	numMolecules.global = numMolecules.local;
#endif

#ifdef ENABLE_MPI
	uint64_t displs = 0;
	MPI_Exscan(&numMolecules.local, &displs, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

	int displs_int = static_cast<int>(displs);
	std::vector<int> displs_vec;
	displs_vec.resize(numProcs);
	MPI_Allgather(&displs_int, 1, MPI_INT, displs_vec.data(), 1, MPI_INT, MPI_COMM_WORLD);

	std::vector<int> recvcounts;
	recvcounts.resize(numProcs);
	MPI_Allgather(&numMolecules.local, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

	// gather particle and component IDs
	{
		CommVar<std::vector<pacIDtype> >& vec = vec_pacID;
		vec.global.resize(numMolecules.global);
		MPI_Allgatherv(vec.local.data(), numMolecules.local, pacID_mpi_type, vec.global.data(), recvcounts.data(),
					   displs_vec.data(), pacID_mpi_type, MPI_COMM_WORLD);
	}
#else
	vec_pacID.global = vec_pacID.local;
#endif

	// shuffle particle and component ID vector (global)
	auto rng = std::default_random_engine{};
	std::shuffle(std::begin(vec_pacID.global), std::end(vec_pacID.global), rng);

	// sort particle IDs to map with respect to component ID
	std::map<int, std::vector<uint64_t> > pidMap;
	for (int32_t i = 0; i < numComponents + 1; i++) {
		const std::vector<uint64_t> vec;
		pidMap.insert(std::pair<int, std::vector<uint64_t> >(i, vec));
	}
	for (auto it : vec_pacID.global) {
		pidMap[it.cid].push_back(it.pid);
		pidMap[0].push_back(it.pid);  // all particles
	}

#ifndef NDEBUG
	std::cout << "[" << domainDecomp->getRank() << "] pidMap: ";
	for (int32_t i = 0; i < numComponents + 1; i++) {
		std::cout << pidMap[i].size() << ' ';
	}
	std::cout << std::endl;
#endif

	// determine particles to swap, delete and add
	std::map<int, std::vector<uint64_t> > swpMap;
	std::map<int, std::vector<uint64_t> > delMap;
	std::map<int, std::vector<uint64_t> > addMap;
	std::vector<int64_t> vecBalance;
	vecBalance.resize(numComponents + 1);

	// update particle count balance
	this->updateBalanceVector(vecBalance, pidMap, numComponents);

	// Prepare: swap components
	{
		uint32_t cid_ub_send;
		uint32_t cid_ub_recv;
		uint32_t prio_send;
		uint32_t prio_recv;
		bool bQuit = false;

		while (not bQuit) {
			cid_ub_send = 0;
			cid_ub_recv = 0;
			prio_send = 0;
			prio_recv = 0;

			for (int32_t i = 1; i < numComponents + 1; i++) {
				const int64_t balance = vecBalance[i];
				const uint32_t prio = _vecPriority[i];
				if (balance > 0 && prio > prio_send) cid_ub_send = i;
				if (balance < 0 && prio > prio_recv) cid_ub_recv = i;
			}

			if (cid_ub_send > 0 && cid_ub_recv > 0) {
				const uint64_t numSwap = std::min(vecBalance[cid_ub_send], abs(vecBalance[cid_ub_recv]));
				Log::global_log->debug() << "[DensityControl] numSwap = " << numSwap << std::endl;
				for (uint64_t i = 0; i < numSwap; i++) {
					const uint64_t pid = pidMap[cid_ub_send].back();
					pidMap[cid_ub_send].pop_back();
					pidMap[cid_ub_recv].push_back(pid);
					swpMap[cid_ub_recv].push_back(pid);
				}
			} else {
				bQuit = true;
			}
			// update particle count balance
			this->updateBalanceVector(vecBalance, pidMap, numComponents);
		}
	}

	// Prepare: delete particles
	{
		for (auto i = 1; i < numComponents + 1; i++) {
			const uint32_t cid_ub = i;
			const int64_t balance = vecBalance[i];
			for (int64_t j = 0; j < balance; j++) {
				const uint64_t pid = pidMap[cid_ub].back();
				delMap[cid_ub].push_back(pid);
				pidMap[cid_ub].pop_back();
			}
		}
		// update particle count balance
		this->updateBalanceVector(vecBalance, pidMap, numComponents);
	}

	// swap components, delete particles
	begin =
		particleContainer->regionIterator(regionLowCorner, regionHighCorner, ParticleIterator::ONLY_INNER_AND_BOUNDARY);
	for (auto it = begin; it.isValid(); ++it) {
		const uint64_t pid = it->getID();

		// check position
		pos[0] = it->r(0);
		pos[1] = it->r(1);
		pos[2] = it->r(2);
		if (not this->moleculeInsideRange(pos)) continue;

		// swap components
		std::vector<Component>* ptrComps = global_simulation->getEnsemble()->getComponents();
		Component* compNew;
		for (uint32_t i = 1; i < numComponents + 1; i++) {
			compNew = &(ptrComps->at(i - 1));

			// search for pid
			for (auto sid : swpMap[i]) {
				if (pid == sid) {
					Log::global_log->debug() << "[DensityControl] Swap pid=" << pid << " with cid=" << it->componentid() + 1
										<< " to cid=" << i << std::endl;
					it->setComponent(compNew);
				}
			}
		}

		// delete particles
		for (uint32_t i = 1; i < numComponents + 1; i++) {
			std::vector<uint64_t>& vec = delMap[i];
			if (std::find(vec.begin(), vec.end(), pid) != vec.end()) {  // found ID of to be deleted particle in vector
				Log::global_log->debug() << "[DensityControl] Delete particle with pid= " << pid
									<< " ; cid_ub= " << it->componentid() + 1 << std::endl;
				particleContainer->deleteMolecule(it, false);
			}
		}
	}  // loop over particles in region

	// add particles
	{
		for (auto cid_ub = 1; cid_ub < numComponents + 1; cid_ub++) {
			int64_t balance = vecBalance[cid_ub];
			if (balance < 0) balance = abs(balance);
			for (int64_t j = 0; j < balance; j++) {
				const uint64_t pid = ++maxMoleculeID.global;
				addMap[cid_ub].push_back(pid);
				pidMap[cid_ub].push_back(pid);
			}
		}
		// update particle count balance
		this->updateBalanceVector(vecBalance, pidMap, numComponents);

		std::vector<Molecule> particles;
		std::vector<Component>* ptrComps = global_simulation->getEnsemble()->getComponents();
		Component* compNew;
		for (int32_t i = 1; i < numComponents + 1; i++) {
			compNew = &(ptrComps->at(i - 1));
			const std::vector<uint64_t>& vec = addMap[i];
			for (auto pid : vec) {
				double rx, ry, rz;
				rx = _range.xmin + _rnd.rnd() * _range.xmax;
				ry = _range.ymin + _rnd.rnd() * _range.ymax;
				rz = _range.zmin + _rnd.rnd() * _range.zmax;
				const Molecule mol(pid, compNew, rx, ry, rz);
				bool bAdd = true;
				bAdd = bAdd && rx > particleContainer->getBoundingBoxMin(0) &&
					   rx < particleContainer->getBoundingBoxMax(0);
				bAdd = bAdd && ry > particleContainer->getBoundingBoxMin(1) &&
					   ry < particleContainer->getBoundingBoxMax(1);
				bAdd = bAdd && rz > particleContainer->getBoundingBoxMin(2) &&
					   rz < particleContainer->getBoundingBoxMax(2);
				if (bAdd) {
					particles.push_back(mol);
				}
			}
		}
		particleContainer->addParticles(particles, true);
	}
}

bool DensityControl::moleculeInsideRange(const std::array<double, 3>& r) {
	return r[0] >= _range.xmin && r[0] < _range.xmax && r[1] >= _range.ymin && r[1] < _range.ymax &&
		   r[2] >= _range.zmin && r[2] < _range.zmax;
}

void DensityControl::updateBalanceVector(std::vector<int64_t>& vecBalance,
										 std::map<int, std::vector<uint64_t> >& pidMap, uint32_t& numComponents) {
	vecBalance[0] = 0;
	for (int32_t i = 1; i < numComponents + 1; i++) {
		const int64_t diff = pidMap[i].size() - _densityTarget.count[i];
		vecBalance[i] = diff;
		vecBalance[0] += diff;
	}

#ifndef NDEBUG
	std::cout << "_densityTarget.count: ";
	for (auto i : _densityTarget.count) std::cout << i << ' ';
	std::cout << std::endl;

	std::cout << "vecBalance: ";
	for (auto i : vecBalance) std::cout << i << ' ';
	std::cout << std::endl;
#endif
}

uint32_t DensityControl::tokenize_int_list(std::vector<uint32_t>& vec, std::string str, std::string del) {
	uint32_t numTokens = vec.size();
	int start = 0;
	int end = str.find(del);
	while (end != -1) {
		const std::string ss = str.substr(start, end - start);
		vec.push_back(atoi(ss.c_str()));
		start = end + del.size();
		end = str.find(del, start);
	}
	const std::string ss = str.substr(start, end - start);
	vec.push_back(atoi(ss.c_str()));
	numTokens = vec.size() - numTokens;
	return numTokens;
}

void DensityControl::initTargetValues(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp) {
	const uint32_t numComponents = global_simulation->getEnsemble()->getComponents()->size();
	double regionLowCorner[3], regionHighCorner[3];

	// if linked cell in the region of interest
	for (unsigned d = 0; d < 3; d++) {
		regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
		regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
	}
	// ensure that we do not iterate over things outside of the container.
	regionLowCorner[0] = std::max(_range.xmin, regionLowCorner[0]);
	regionLowCorner[1] = std::max(_range.ymin, regionLowCorner[1]);
	regionLowCorner[2] = std::max(_range.zmin, regionLowCorner[2]);

	regionHighCorner[0] = std::min(_range.xmax, regionHighCorner[0]);
	regionHighCorner[1] = std::min(_range.ymax, regionHighCorner[1]);
	regionHighCorner[2] = std::min(_range.zmax, regionHighCorner[2]);

	CommVar<std::vector<uint64_t> > numMolecules;
	numMolecules.local.resize(numComponents + 1);
	numMolecules.global.resize(numComponents + 1);
	for (uint32_t i = 0; i < numComponents + 1; i++) {
		numMolecules.local.at(i) = numMolecules.global.at(i) = 0;
	}
	std::array<double, 3> pos;

	auto begin = particleContainer->regionIterator(regionLowCorner, regionHighCorner,
												   ParticleIterator::ONLY_INNER_AND_BOUNDARY);  // over all cell types
	for (auto it = begin; it.isValid(); ++it) {
		// check position
		pos[0] = it->r(0);
		pos[1] = it->r(1);
		pos[2] = it->r(2);

		if (not this->moleculeInsideRange(pos)) continue;

		const uint32_t cid_ub = it->componentid() + 1;
		numMolecules.local.at(cid_ub)++;
		numMolecules.local.at(0)++;
	}

#ifdef ENABLE_MPI
	domainDecomp->collCommInit(numComponents + 1);
	for (uint32_t i = 0; i < numComponents + 1; i++) {
		domainDecomp->collCommAppendUnsLong(numMolecules.local.at(i));
	}
	domainDecomp->collCommAllreduceSum();
	for (uint32_t i = 0; i < numComponents + 1; i++) {
		numMolecules.global.at(i) = domainDecomp->collCommGetUnsLong();
	}
	domainDecomp->collCommFinalize();
#else
	for (uint32_t i = 0; i < numComponents + 1; i++) {
		numMolecules.global[i] = numMolecules.local[i];
	}
#endif

	// Set target values that were omitted in XML config
	const double invVolume = 1. / _range.volume;

	for (uint32_t i = 1; i < numComponents + 1; i++) {
		if (_densityTarget.specified.at(i)) {
			continue;
		}
		_densityTarget.count.at(i) = numMolecules.global.at(i);
		_densityTarget.density.at(i) = _densityTarget.count.at(i) * invVolume;
	}

	// calc total target density and particle count
	uint64_t countSum = 0;
	for (int32_t i = 1; i < numComponents + 1; i++) {
		countSum += _densityTarget.count[i];
	}
	_densityTarget.count[0] = countSum;
	_densityTarget.density[0] = countSum * invVolume;

#ifndef NDEBUG
	// show target values
	for (uint32_t i = 0; i < numComponents + 1; i++) {
		Log::global_log->info() << "cID=" << i << ": target count=" << _densityTarget.count[i]
						   << ", target density=" << _densityTarget.density[i] << std::endl;
	}
#endif
}
