/**
 * @file LoadBalancer.h
 * @author seckler
 * @date 04.06.19
 */

#pragma once
#include <array>
#include <tuple>

class LoadBalancer {
public:
	virtual ~LoadBalancer() = default;
	virtual std::tuple<std::array<double, 3>, std::array<double, 3>> rebalance(double work) = 0;
};
