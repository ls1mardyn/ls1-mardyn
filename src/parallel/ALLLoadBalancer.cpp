/**
 * @file ALLLoadBalancer.cpp
 * @author seckler
 * @date 04.06.19
 */

#include "ALLLoadBalancer.h"
ALLLoadBalancer::ALLLoadBalancer(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, double gamma,
								 MPI_Comm comm, const std::array<size_t, 3> &globalSize,
								 const std::array<size_t, 3> &localCoordinates,
								 const std::array<double, 3> &minimalPartitionSize)
	: _all(3 /*dim*/, gamma), _minimalPartitionSize(minimalPartitionSize) {
	// convert input into non-const vector because that is what ALL expects
	std::vector<Point> points {
		{3, boxMin.data()},
		{3, boxMax.data()},
	};
	_all.set_vertices(points);
	// convert input into non-const int arrays because that is what ALL expects
	std::array<int, 3> globalSizeIntArray{static_cast<int>(globalSize[0]), static_cast<int>(globalSize[1]),
								   static_cast<int>(globalSize[2])};
	std::array<int, 3> coords{static_cast<int>(localCoordinates[0]), static_cast<int>(localCoordinates[1]),
							  static_cast<int>(localCoordinates[2])};
	_all.set_proc_grid_params(coords.data(), globalSizeIntArray.data());
	_all.set_communicator(comm);

	_coversWholeDomain = {globalSizeIntArray[0] == 1, globalSizeIntArray[1] == 1, globalSizeIntArray[2] == 1};
}
std::tuple<std::array<double, 3>, std::array<double, 3>> ALLLoadBalancer::rebalance(double work) {
	_all.set_work(work);
	_all.setup(ALL_LB_t::STAGGERED);
	_all.set_min_domain_size(ALL_LB_t::STAGGERED, _minimalPartitionSize.data());
	_all.balance(ALL_LB_t::STAGGERED);
	auto resultVertices = _all.get_result_vertices();
	_all.set_vertices(resultVertices);
	const std::array<double, 3> boxMin{resultVertices[0].x(0), resultVertices[0].x(1), resultVertices[0].x(2)};
	const std::array<double, 3> boxMax{resultVertices[1].x(0), resultVertices[1].x(1), resultVertices[1].x(2)};
	return std::make_tuple(boxMin, boxMax);
}
