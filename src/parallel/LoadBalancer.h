/**
 * @file LoadBalancer.h
 * @author seckler
 * @date 04.06.19
 */

#pragma once
#include <array>
#include <tuple>

/**
 * LoadBalancer class for the usage of arbitrary load balancing classes that are handled by GeneralDomainDecomposition.
 */
class LoadBalancer {
public:
	/**
	 * Virtual destructor.
	 */
	virtual ~LoadBalancer() = default;

	/**
	 * The rebalancing call.
	 * @param work Arbitrary unit of work, e.g., time for the current process
	 * @return New decomposition. First entry is the new boxMin, second the new boxMax.
	 */
	virtual std::tuple<std::array<double, 3>, std::array<double, 3>> rebalance(double work) = 0;
};
