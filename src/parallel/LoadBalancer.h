/**
 * @file LoadBalancer.h
 * @author seckler, Georg von Bismarck
 * @date 29.09.2026
 */

#pragma once
#include <array>
#include "utils/xmlfileUnits.h"
#include "parallel/DomainDecompMPIBase.h"

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
	 * @param localSupdomain The local Supdomain of the Prozess. It does not have to correspond exactly to the actually implemented subdomain.  
	 * @param work Arbitrary unit of work, e.g., time for the current process
	 * @return New domain boundaries for the current process.
	 */
	virtual DomainBox rebalance(DomainBox localBox, double work) = 0;

	/**
	 * Read Config file
	 * @param xmlconfig
	 */
	virtual void readXML(XMLfileUnits& xmlconfig) = 0;

	/**
	 * Indicates if the current process / MPI rank spans the full length of a dimension.
	 * @return Array of bools, for each dimension one value: true, iff the process spans the entire domain along this dimension.
	 */
	virtual std::array<bool, 3> getCoversWholeDomain() = 0;
};
