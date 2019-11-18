/**
 * @file LoadBalancer.h
 * @author seckler
 * @date 04.06.19
 */

#pragma once
#include <array>
#include <tuple>
#include "utils/xmlfileUnits.h"

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
	 * Based on the current domain and the work for that domain this function determines a new
	 * domain decomposition that provides a better load balancing.
	 * This call will normally include communication and exchange of information with other processes.
	 * @param work Arbitrary unit of work, e.g., time for the current process
	 * @return New domain boundaries for the current process. First entry is the new boxMin,
	 * second the new boxMax.
	 */
	virtual std::tuple<std::array<double, 3>, std::array<double, 3>> rebalance(double work) = 0;

	/**
	 * Read Config file
	 * @param xmlconfig
	 */
	virtual void readXML(XMLfileUnits& xmlconfig) = 0;
};
