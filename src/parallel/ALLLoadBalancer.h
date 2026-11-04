/**
 * @file ALLLoadBalancer.h
 * @author seckler, Georg von Bismarck
 * @date 4.11.2026
 */

#pragma once
#ifdef ENABLE_ALLLBL
#include <ALL.hpp>
#include "LoadBalancer.h"

/** 
 * Integration of A Loadbalancing Library (ALL) into LS1 for use through GeneralDomainDecomposition
 */ 
class ALLLoadBalancer : public LoadBalancer {
public:

	/**
     * Constructs from two corner points.
     * @param comm !Cartesian-Communicator! 
	 * @param globalSize global Domain Size
	 * @param minimalDomainSize
     */
	ALLLoadBalancer(MPI_Comm comm, DomainGridPoint globalSize, std::vector<double> minimalDomainSize) noexcept;
	~ALLLoadBalancer() override = default;

	/** 
	 * documentation see base class (LoadBalancer.h)
	 */ 
	DomainBox rebalance(DomainBox localBox, double work) override;

	/** 
	 * documentation see base class (LoadBalancer.h)
	 */ 
	void readXML(XMLfileUnits& xmlconfig) override;

	/** 
	 * documentation see base class (LoadBalancer.h)
	 */ 
	std::array<bool, 3> getCoversWholeDomain() noexcept override { return _coversWholeDomain; }

private:
	std::unique_ptr<ALL::ALL<double, double>> _all;
	MPI_Comm _comm;

	std::vector<double> _minimalDomainSize{};
	std::array<bool, 3> _coversWholeDomain{};
};
#endif
