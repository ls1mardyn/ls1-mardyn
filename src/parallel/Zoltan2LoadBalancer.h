/**
 * @file Zoltan2LoadBalancer.h
 * @author seckler
 * @date 14.11.19
 */
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
	Zoltan2LoadBalancer();

	/**
	 * default destructor.
	 */
	~Zoltan2LoadBalancer() override = default;

	// doc see base class.
	std::tuple<std::array<double, 3>, std::array<double, 3>> rebalance(double work) override;
};
