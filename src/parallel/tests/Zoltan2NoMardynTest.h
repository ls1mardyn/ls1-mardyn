/**
 * @file Zoltan2NoMardynTest.h
 * @author seckler
 * @date 18.11.19
 */

#pragma once

#include <utils/Testing.h>
#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>

class Zoltan2NoMardynTest : public utils::Test {
	TEST_SUITE(Zoltan2NoMardynTest);
	TEST_METHOD(multiJaggedTest);
	TEST_METHOD(multiJaggedTest2);
	TEST_SUITE_END();

public:
	Zoltan2NoMardynTest() = default;
	~Zoltan2NoMardynTest() override = default;

	static void multiJaggedTest();

	static void multiJaggedTest2();
private:
	static void printArray(const std::string& info, double* array, size_t count, int rank);

	template<typename inputAdapter_t>
	static void printDecomposition(int rank, int nprocs, int dim,
								   Zoltan2::PartitioningProblem<inputAdapter_t>* problem){
		auto view = problem->getSolution().getPartBoxesView();

		if (rank == 0) {
			for (int rankid = 0; rankid < nprocs; ++rankid) {
				auto lmins = view[rankid].getlmins();
				auto lmaxs = view[rankid].getlmaxs();
				std::cout << "rank " << rankid << ": ";
				for (int i = 0; i < dim; ++i) {
					std::cout << "[" << lmins[i] << ", " << lmaxs[i] << (i < dim - 1 ? "] x " : "]");
				}
				std::cout << std::endl;
			}
		}
	}

	template<typename inputAdapter_t>
	static void checkMetric(int rank, double tolerance, Teuchos::ParameterList& params,
	                        inputAdapter_t& ia,
							Zoltan2::PartitioningProblem<inputAdapter_t>& problem){
		Zoltan2::EvaluatePartition<inputAdapter_t> metricObject (&ia, &params, problem.getComm(),
		                                                           &problem.getSolution());
		if (rank == 0){
			metricObject.printMetrics(std::cout);
			double imb = metricObject.getWeightImbalance(0);
			if (imb <= tolerance)
				std::cout << "pass: " << imb << std::endl;
			else
				std::cout << "fail: " << imb << std::endl;
			std::cout << std::endl;
		}
	}
};
