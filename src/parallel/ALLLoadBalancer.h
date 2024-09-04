/**
 * @file ALLLoadBalancer.h
 * @author seckler
 * @date 04.06.19
 */

#pragma once
#ifdef ENABLE_ALLLBL
#include <ALL.hpp>
#include "LoadBalancer.h"

class ALLLoadBalancer : public LoadBalancer {
public:
	ALLLoadBalancer(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, double gamma, MPI_Comm comm,
					const std::array<size_t, 3>& globalSize, const std::array<size_t, 3>& localCoordinates,
					const std::array<double, 3>& minimalPartitionSize);

	~ALLLoadBalancer() override = default;
	std::tuple<std::array<double, 3>, std::array<double, 3>> rebalance(double work) override;
	void readXML(XMLfileUnits& xmlconfig) override {
		// nothing yet.
	}

	const std::array<bool, 3>& getCoversWholeDomain() const override { return _coversWholeDomain; }

private:
	ALL<double, double> _all;
	using Point = ALL_Point<double>;
	std::array<double, 3> _minimalPartitionSize{};
	std::array<bool, 3> _coversWholeDomain{};
};
#endif
