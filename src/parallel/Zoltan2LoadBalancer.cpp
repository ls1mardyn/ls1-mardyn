/**
 * @file Zoltan2LoadBalancer.cpp
 * @author seckler
 * @date 14.11.19
 */

#include "Zoltan2LoadBalancer.h"

#include <Zoltan2_Adapter.hpp>
#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_MultiJagged_ReductionOps.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include "utils/Logger.h"

Zoltan2LoadBalancer::Zoltan2LoadBalancer(std::array<double, 3> boxMin, std::array<double, 3> boxMax, MPI_Comm comm,
										 double /*minimalDomainSize*/, std::array<double, 3> domainLength)
	: _comm{comm}, _boxMin{boxMin}, _boxMax{boxMax}, _domainLength{domainLength} {
	MPI_Comm_size(comm, &_numRanks);
	MPI_Comm_rank(comm, &_rank);
	Log::global_log->warning() << "Zoltan2LoadBalancer: minimalDomainSize currently ignored!" << std::endl;
	// TODO: don't ignore minimal domain size
}

typedef Tpetra::Map<> Map_t;
typedef Map_t::local_ordinal_type localId_t;
typedef Map_t::global_ordinal_type globalId_t;
typedef double scalar_t;
typedef Zoltan2::BasicUserTypes<scalar_t, localId_t, globalId_t> myTypes;
typedef Zoltan2::BasicVectorAdapter<myTypes> inputAdapter_t;

std::tuple<std::array<double, 3>, std::array<double, 3>> Zoltan2LoadBalancer::rebalance(double work) {
	// Number of local samples.
	size_t localCount = 10;
	if (_rank == 0) {
		localCount = localCount + 2;  // min and max of domain, s.t., it will remain
	}

	std::vector<scalar_t> x(localCount);
	std::vector<scalar_t> y(localCount);
	std::vector<scalar_t> z(localCount);

	// Create coordinates that range from 0 to 10.0

	srand(_rank);
	scalar_t scalingFactorX = (_boxMax[0] - _boxMin[0]) / RAND_MAX;
	scalar_t scalingFactorY = (_boxMax[1] - _boxMin[1]) / RAND_MAX;
	scalar_t scalingFactorZ = (_boxMax[2] - _boxMin[2]) / RAND_MAX;

	for (size_t i = 0; i < localCount; i++) {
		x[i] = scalar_t(rand()) * scalingFactorX + _boxMin[0];
		y[i] = scalar_t(rand()) * scalingFactorY + _boxMin[1];
		z[i] = scalar_t(rand()) * scalingFactorZ + _boxMin[2];
	}
	std::vector<scalar_t> weights(localCount);
	for (size_t i = 0; i < localCount; i++) {
		weights[i] = work;
	}

	if (_rank == 0) {
		x[0] = y[0] = z[0] = 0;
		x[1] = _domainLength[0];
		y[1] = _domainLength[1];
		z[1] = _domainLength[2];
		weights[0] = 0.;
		weights[1] = 0.;
	}

	// Create global ids for the coordinates.

	std::vector<globalId_t> globalIds(localCount);
	globalId_t offset = _rank * localCount + (_rank ? 2 : 0);

	for (size_t i = 0; i < localCount; i++) {
		globalIds[i] = offset++;
	}

	// Create a Zoltan2 input adapter for this geometry. TODO explain

	std::vector<const scalar_t *> coordVec(3);
	std::vector<int> coordStrides(3);

	coordVec[0] = x.data();
	coordStrides[0] = 1;
	coordVec[1] = y.data();
	coordStrides[1] = 1;
	coordVec[2] = z.data();
	coordStrides[2] = 1;

	// 1-d array as multiple weights are supported
	std::vector<const scalar_t *> weightVec(1);
	// 1-d array as multiple weights are supported
	std::vector<int> weightStrides(1);

	// point to weights vector
	weightVec[0] = weights.data();
	weightStrides[0] = 1;

	inputAdapter_t ia(localCount, globalIds.data(), coordVec, coordStrides, weightVec, weightStrides);

	Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, &_params, _comm);

	problem.solve();

	auto view = problem.getSolution().getPartBoxesView();

	auto lmins = view[_rank].getlmins();
	auto lmaxs = view[_rank].getlmaxs();

	_boxMin = {lmins[0], lmins[1], lmins[2]};
	_boxMax = {lmaxs[0], lmaxs[1], lmaxs[2]};
	return {_boxMin, _boxMax};
}

void Zoltan2LoadBalancer::readXML(XMLfileUnits &xmlconfig) {
	// set the algorithm
	_params.set("algorithm", "multijagged");

	// don't allow particles on one line to belong to different partitions
	_params.set("rectilinear", true);

	// needed to be able to get the boundary boxes.
	_params.set("mj_keep_part_boxes", true);

	// set the tolerance of imbalance, should be > 1.
	_params.set("imbalance_tolerance", 1.1);

	// option does not exist in trilinos 12.12.1
	// params->set("mj_premigration_option", mj_premigration_option);

	// predefined partitioning grid, e.g., 2x2x4 or so.
	// if (pqParts != "") params->set("mj_parts", pqParts);

	// set the number of ranks as the number of partitions.
	_params.set("num_global_parts", _numRanks);

	// if (migration_check_option >= 0) params->set("mj_migration_option", migration_check_option);

	// if (migration_imbalance_cut_off >= 0)
	//	params->set("mj_minimum_migration_imbalance", double(migration_imbalance_cut_off));
}
