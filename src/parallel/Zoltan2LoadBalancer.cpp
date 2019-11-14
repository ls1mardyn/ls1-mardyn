/**
 * @file Zoltan2LoadBalancer.cpp
 * @author seckler
 * @date 14.11.19
 */

#include "Zoltan2LoadBalancer.h"

#include <Zoltan2_Adapter.hpp>
#include <Zoltan2_MultiJagged_ReductionOps.hpp>

Zoltan2LoadBalancer::Zoltan2LoadBalancer() = default;

std::tuple<std::array<double, 3>, std::array<double, 3>> Zoltan2LoadBalancer::rebalance(double work) {
	return std::tuple<std::array<double, 3>, std::array<double, 3>>();
}
