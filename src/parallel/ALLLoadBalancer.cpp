/**
 * @file ALLLoadBalancer.cpp
 * @author seckler
 * @date 04.06.19
 */

#include "ALLLoadBalancer.h"
ALLLoadBalancer::ALLLoadBalancer(std::array<double, 3> boxMin, std::array<double, 3> boxMax, double gamma,
								 MPI_Comm comm, std::array<size_t, 3> globalSize,
								 std::array<size_t, 3> localCoordinates, double minimalPartitionSize)
	: _all(3 /*dim*/, gamma) {
	std::vector<Point> points;
	points.emplace_back(3, boxMin.data());
	points.emplace_back(3, boxMax.data());
	_all.set_vertices(points);
	std::array<int, 3> global_size{static_cast<int>(globalSize[0]), static_cast<int>(globalSize[1]),
								   static_cast<int>(globalSize[2])};
	std::array<int, 3> coords{static_cast<int>(localCoordinates[0]), static_cast<int>(localCoordinates[1]),
							  static_cast<int>(localCoordinates[2])};
	_all.set_proc_grid_params(coords.data(), global_size.data());
	_all.set_communicator(comm);

	_coversWholeDomain = {globalSize[0] == 1, global_size[1] == 1, global_size[2] == 1};

	_minimalPartitionSize = minimalPartitionSize;
}
std::tuple<std::array<double, 3>, std::array<double, 3>> ALLLoadBalancer::rebalance(double work) {
	_all.set_work(work);
	_all.setup(ALL_LB_t::STAGGERED);
	std::array<double, 3> minDomainSize{_minimalPartitionSize, _minimalPartitionSize, _minimalPartitionSize};
	_all.set_min_domain_size(ALL_LB_t::STAGGERED, minDomainSize.data());
	_all.balance(ALL_LB_t::STAGGERED);
	auto resultVertices = _all.get_result_vertices();
	std::array<double, 3> boxMin{resultVertices[0].x(0), resultVertices[0].x(1), resultVertices[0].x(2)};
	std::array<double, 3> boxMax{resultVertices[1].x(0), resultVertices[1].x(1), resultVertices[1].x(2)};
	_all.set_vertices(resultVertices);
	return std::make_tuple(boxMin, boxMax);
}
