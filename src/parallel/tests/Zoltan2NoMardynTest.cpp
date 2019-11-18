/**
 * @file Zoltan2NoMardynTest.cpp
 * @author seckler
 * @date 18.11.19
 */

#include "Zoltan2NoMardynTest.h"

#include <Tpetra_Map.hpp>
#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <cstdlib>
#include <vector>

TEST_SUITE_REGISTRATION(Zoltan2NoMardynTest);

void Zoltan2NoMardynTest::multiJaggedTest() {
#ifdef HAVE_ZOLTAN2_MPI
	int rank, nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	int rank = 0, nprocs = 1;
#endif

	// For convenience, we'll use the Tpetra defaults for local/global ID types
	// Users can substitute their preferred local/global ID types
	typedef Tpetra::Map<> Map_t;
	typedef Map_t::local_ordinal_type localId_t;
	typedef Map_t::global_ordinal_type globalId_t;

	typedef double scalar_t;
	typedef Zoltan2::BasicUserTypes<scalar_t, localId_t, globalId_t> myTypes;

	// TODO explain
	typedef Zoltan2::BasicVectorAdapter<myTypes> inputAdapter_t;
	typedef Zoltan2::EvaluatePartition<inputAdapter_t> quality_t;
	typedef inputAdapter_t::part_t part_t;

	///////////////////////////////////////////////////////////////////////
	// Create input data.

	size_t localCount = 40;
	int dim = 3;

	std::vector<scalar_t> coords(dim * localCount);

	scalar_t *x = coords.data();
	scalar_t *y = x + localCount;
	scalar_t *z = y + localCount;

	// Create coordinates that range from 0 to 10.0

	srand(rank);
	scalar_t scalingFactor = 10.0 / RAND_MAX;

	for (size_t i = 0; i < localCount * dim; i++) {
		coords[i] = scalar_t(rand()) * scalingFactor;
	}
	if (rank == 0) {
		x[0] = y[0] = z[0] = 0;
		x[1] = y[1] = z[1] = 10.;
	}

	printArray("x coords:", x, localCount, rank);

	// Create global ids for the coordinates.

	std::vector<globalId_t> globalIds(localCount);
	globalId_t offset = rank * localCount;

	for (size_t i = 0; i < localCount; i++) globalIds[i] = offset++;

	///////////////////////////////////////////////////////////////////////
	// Create parameters for an MJ problem

	double tolerance = 1.1;

	if (rank == 0) std::cout << "Imbalance tolerance is " << tolerance << std::endl;

	Teuchos::ParameterList params("test params");
	params.set("debug_level", "basic_status");
	params.set("debug_procs", "0");
	params.set("error_check_level", "debug_mode_assertions");

	params.set("algorithm", "multijagged");
	params.set("imbalance_tolerance", tolerance);
	params.set("num_global_parts", nprocs);
	params.set("mj_keep_part_boxes", true);

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	// A simple problem with no weights.
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////

	// Create a Zoltan2 input adapter for this geometry. TODO explain

	inputAdapter_t ia1(localCount, globalIds.data(), x, y, z, 1, 1, 1);

	// Create a Zoltan2 partitioning problem

	Zoltan2::PartitioningProblem<inputAdapter_t> problem1(&ia1, &params);

	// Solve the problem

	problem1.solve();

	printDecomposition(rank, nprocs, dim, &problem1);

	checkMetric(rank, tolerance, params, ia1, problem1);

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	// Try a problem with weights -- this only uses 2 dimensions (x,y) of the previous input.
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////

	std::vector<scalar_t> weights(localCount);
	for (size_t i = 0; i < localCount; i++) {
		weights[i] = 1.0 + scalar_t(rank) / scalar_t(nprocs);
	}

	printArray("weights:", weights.data(), localCount, rank);

	// Create a Zoltan2 input adapter that includes weights.

	std::vector<const scalar_t *> coordVec(2);
	std::vector<int> coordStrides(2);

	coordVec[0] = x;
	coordStrides[0] = 1;
	coordVec[1] = y;
	coordStrides[1] = 1;

	std::vector<const scalar_t *> weightVec(1);
	std::vector<int> weightStrides(1);

	weightVec[0] = weights.data();
	weightStrides[0] = 1;

	inputAdapter_t ia2(localCount, globalIds.data(), coordVec, coordStrides, weightVec, weightStrides);

	// Create a Zoltan2 partitioning problem

	Zoltan2::PartitioningProblem<inputAdapter_t> problem2(&ia2, &params);

	// Solve the problem

	problem2.solve();

	printDecomposition(rank, nprocs, 2 /*only 2 dimensional*/, &problem2);

	// create metric object for MPI builds
	checkMetric(rank, tolerance, params, ia2, problem2);

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	// Try a problem with multiple weights.
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////

	// Add to the parameters the multicriteria objective.

	params.set("partitioning_objective", "multicriteria_minimize_total_weight");

	// Create the new weights.

	std::vector<scalar_t> weights3(localCount * 3);
	srand(rank);

	for (size_t i = 0; i < localCount * 3; i += 3) {
		weights3[i] = 1.0 + scalar_t(rank) / nprocs;         // weight idx 1
		weights3[i + 1] = rank < nprocs / 2 ? 1 : 2;         // weight idx 2
		weights3[i + 2] = scalar_t(rand()) / RAND_MAX + .5;  // weight idx 3
	}

	// Create a Zoltan2 input adapter with these weights.

	weightVec.resize(3);
	weightStrides.resize(3);

	weightVec[0] = weights3.data();
	weightStrides[0] = 3;
	weightVec[1] = weights3.data() + 1;
	weightStrides[1] = 3;
	weightVec[2] = weights3.data() + 2;
	weightStrides[2] = 3;

	inputAdapter_t ia3(localCount, globalIds.data(), coordVec, coordStrides, weightVec, weightStrides);

	// Create a Zoltan2 partitioning problem.

	Zoltan2::PartitioningProblem<inputAdapter_t> problem3(&ia3, &params);

	// Solve the problem

	problem3.solve();

	// check solution

	checkMetric(rank, tolerance, params, ia3, problem3);

	///////////////////////////////////////////////////////////////////////
	// Try the other multicriteria objectives.

	bool dataHasChanged = false;  // default is true

	params.set("partitioning_objective", "multicriteria_minimize_maximum_weight");
	problem3.resetParameters(&params);
	problem3.solve(dataHasChanged);

	// Solution changed!

	checkMetric(rank, tolerance, params, ia3, problem3);

	params.set("partitioning_objective", "multicriteria_balance_total_maximum");
	problem3.resetParameters(&params);
	problem3.solve(dataHasChanged);

	// Solution changed!
	checkMetric(rank, tolerance, params, ia3, problem3);

	if (rank == 0) std::cout << "PASS" << std::endl;
}

void Zoltan2NoMardynTest::multiJaggedTest2() {
#ifdef HAVE_ZOLTAN2_MPI
	int rank, nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	int rank = 0, nprocs = 1;
#endif

	// For convenience, we'll use the Tpetra defaults for local/global ID types
	// Users can substitute their preferred local/global ID types
	typedef Tpetra::Map<> Map_t;
	typedef Map_t::local_ordinal_type localId_t;
	typedef Map_t::global_ordinal_type globalId_t;

	typedef double scalar_t;
	typedef Zoltan2::BasicUserTypes<scalar_t, localId_t, globalId_t> myTypes;

	// TODO explain
	typedef Zoltan2::BasicVectorAdapter<myTypes> inputAdapter_t;
	typedef Zoltan2::EvaluatePartition<inputAdapter_t> quality_t;
	typedef inputAdapter_t::part_t part_t;

	///////////////////////////////////////////////////////////////////////
	// Create input data.

	size_t localCount = 3;
	int dim = 3;

	std::vector<scalar_t> coords(dim * localCount);

	scalar_t *x = coords.data();
	scalar_t *y = x + localCount;
	scalar_t *z = y + localCount;

	// Create coordinates that range from 0 to 10.0

	srand(rank);
	scalar_t scalingFactor = 10.0 / RAND_MAX;

	for (size_t i = 0; i < localCount * dim; i++) {
		coords[i] = scalar_t(rand()) * scalingFactor;
	}
	std::vector<scalar_t> weights(localCount);
	for (size_t i = 0; i < localCount; i++) {
		weights[i] = 1.0 + scalar_t(rank) / scalar_t(nprocs);
	}

	if (rank == 0) {
		x[0] = y[0] = z[0] = 0;
		x[1] = y[1] = z[1] = 10.;
		// weights[0] = 0.;
		// weights[1] = 0.;
	}

	printArray("x coords: ", x, localCount, rank);
	printArray("weights: ", weights.data(), localCount, rank);

	// Create global ids for the coordinates.

	std::vector<globalId_t> globalIds(localCount);
	globalId_t offset = rank * localCount;

	for (size_t i = 0; i < localCount; i++) globalIds[i] = offset++;

	///////////////////////////////////////////////////////////////////////
	// Create parameters for an MJ problem

	double tolerance = 1.1;

	if (rank == 0) {
		std::cout << "Imbalance tolerance is " << tolerance << std::endl;
	}

	Teuchos::ParameterList params("test params");
	params.set("debug_level", "basic_status");
	params.set("debug_procs", "0");
	params.set("error_check_level", "debug_mode_assertions");

	params.set("algorithm", "multijagged");
	params.set("imbalance_tolerance", tolerance);
	params.set("num_global_parts", nprocs);
	params.set("mj_keep_part_boxes", true);
	params.set("rectilinear", true);

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	// A simple problem with one weight.
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////

	// Create a Zoltan2 input adapter for this geometry. TODO explain

	std::vector<const scalar_t *> coordVec(3);
	std::vector<int> coordStrides(3);

	coordVec[0] = x;
	coordStrides[0] = 1;
	coordVec[1] = y;
	coordStrides[1] = 1;
	coordVec[2] = z;
	coordStrides[2] = 1;

	// 1-d array as multiple weights are supported
	std::vector<const scalar_t *> weightVec(1);
	// 1-d array as multiple weights are supported
	std::vector<int> weightStrides(1);

	// point to weights vector
	weightVec[0] = weights.data();
	weightStrides[0] = 1;

	inputAdapter_t ia1(localCount, globalIds.data(), coordVec, coordStrides, weightVec, weightStrides);

	// Create a Zoltan2 partitioning problem

	Zoltan2::PartitioningProblem<inputAdapter_t> problem1(&ia1, &params);

	// Solve the problem

	problem1.solve();

	printDecomposition(rank, nprocs, 3, &problem1);
	// create metric object where communicator is Teuchos default

	checkMetric(rank, tolerance, params, ia1, problem1);

	if (rank == 0) {
		std::cout << "PASS" << std::endl;
	}
}

void Zoltan2NoMardynTest::printArray(const std::string &info, double *array, size_t count, int rank) {
	if (rank == 0) {
		std::cout << info << std::flush;
	}
#ifdef HAVE_ZOLTAN2_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	for (size_t i = 0; i < count; i++) {
		std::cout << array[i] << ", ";
	}
	std::cout << std::flush;
#ifdef HAVE_ZOLTAN2_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	if (rank == 0) {
		std::cout << std::endl;
	}
}
