/**
 * @file ALLLoadBalancer.h
 * @author seckler, Georg von Bismarck
 * @date 23.09.2026
 */

#pragma once
#ifdef ENABLE_ALLLBL
#include <ALL.hpp>
#include "LoadBalancer.h"
#include "parallel/DomainDecompMPIBase.h"

#include <tuple>
class ALLLoadBalancer : public LoadBalancer {
public:
	ALLLoadBalancer(double gamma, MPI_Comm comm, DomainGridPoint globalSize, std::vector<double> minimalPartitionSize);

	~ALLLoadBalancer() override = default;
	DomainBox rebalance(DomainBox localBox, double work) override;
	void readXML(XMLfileUnits& xmlconfig) override;

	std::array<bool, 3> getCoversWholeDomain() override { return _coversWholeDomain; }

private:
	std::unique_ptr<ALL::ALL<double, double>> _all;
	MPI_Comm _comm;
	double _gamma;

	std::vector<double> _minimalPartitionSize{};
	std::array<bool, 3> _coversWholeDomain{};
};
#endif
