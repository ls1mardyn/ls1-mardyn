/**
 * @file Zoltan2LoadBalancer.h
 * @author seckler
 * @date 14.11.19
 */
#include <mpi.h>
#include <Teuchos_ParameterList.hpp>
#include "LoadBalancer.h"

#pragma once

/**
 * Load balancer to interface the Zoltan2 package.
 * Uses the MultiJagged algorithm provided by Zoltan2, which is similar to the multisection method.
 */
class Zoltan2LoadBalancer : public LoadBalancer {
public:
	/**
	 * Constructor
	 */
	Zoltan2LoadBalancer(std::array<double, 3> boxMin, std::array<double, 3> boxMax, MPI_Comm comm,
						double minimalDomainSize, std::array<double, 3> domainLength);

	/**
	 * default destructor.
	 */
	~Zoltan2LoadBalancer() override = default;

	// doc see base class.
	std::tuple<std::array<double, 3>, std::array<double, 3>> rebalance(double work) override;

	void readXML(XMLfileUnits& xmlconfig) override;

private:
	Teuchos::ParameterList _params;

private:
	MPI_Comm _comm;
	std::array<double, 3> _boxMin;
	std::array<double, 3> _boxMax;
	std::array<double, 3> _domainLength;
	int _numRanks{1};
	int _rank{0};
};
