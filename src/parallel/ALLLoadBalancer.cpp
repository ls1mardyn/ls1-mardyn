/**
 * @file ALLLoadBalancer.cpp
 * @author seckler, Georg von Bismarck
 * @date 4.11.2026
 */

#include "ALLLoadBalancer.h"
#include <string>
#include "ALL.hpp"
#include "utils/mardyn_assert.h"

ALLLoadBalancer::ALLLoadBalancer(MPI_Comm comm, DomainGridPoint globalSize, std::vector<double> minimalDomainSize) noexcept {
	_comm = comm;
	_minimalDomainSize = minimalDomainSize;
	
	_coversWholeDomain = {globalSize[0] == 1, globalSize[1] == 1, globalSize[2] == 1};;
}

void ALLLoadBalancer::readXML(XMLfileUnits& xmlconfig){
	ALL::LB_t mode = ALL::LB_t::UNIMPLEMENTED;
	double gamma = 4;

	std::string loadBalancer("STAGGERED");
	xmlconfig.getNodeValue("mode", loadBalancer);
	
	if (loadBalancer == "STAGGERED") { //default
		mode = ALL::LB_t::STAGGERED;
	} else if (loadBalancer == "TENSOR") {
		mode = ALL::LB_t::TENSOR;
	} else if (loadBalancer == "FORCEBASED") {
		mode = ALL::LB_t::FORCEBASED; // has not been fully tested in LS1-Mardyn and may produce unexpected performance results
	} else if (loadBalancer == "ALL_VORONOI_ACTIVE") {
		#ifdef ALL_VORONOI_ACTIVE
			mode = ALL::LB_t::VORONOI; // has not been fully tested in LS1-Mardyn and may produce unexpected performance results
		#else
			std::ostringstream error_message;
			error_message << "ALLLoadBalancer: ALL libery has VORONOI not active. Aborting! Please select a valid option!";
			MARDYN_EXIT(error_message.str());
		#endif
	} else if (loadBalancer == "HISTOGRAM") {
		mode = ALL::LB_t::HISTOGRAM;  // has not been fully tested in LS1-Mardyn and may produce unexpected performance results
	} else if (loadBalancer == "TENSOR_MAX") {
		mode = ALL::LB_t::TENSOR_MAX;  // has not been fully tested in LS1-Mardyn and may produce unexpected performance results
	} else {
		std::ostringstream error_message;
		error_message << "ALLLoadBalancer: Unsupported load balancer " << loadBalancer << " was selected. Aborting! Please select a valid option!";
		MARDYN_EXIT(error_message.str());
	}

	//gamma is used as a correction factor if required by the selected balancing method
	xmlconfig.getNodeValue("gamma", gamma);

	Log::global_log->info() << "ALLLoadBalancer: using the " << loadBalancer << " load balancer with gamma of " << gamma << std::endl;

	_all = std::make_unique<ALL::ALL<double, double>>(ALL::TENSOR, DIMgeom, gamma);
	_all->setCommunicator(_comm);
	_all->setMinDomainSize(_minimalDomainSize);
    _all->setup();
}

DomainBox ALLLoadBalancer::rebalance(DomainBox localBox, double work) {
	std::vector<ALL::Point<double>> domain(2, ALL::Point<double>(DIMgeom));

	domain[0][0] = localBox.min[0];
	domain[0][1] = localBox.min[1];
	domain[0][2] = localBox.min[2];

	domain[1][0] = localBox.max[0];
	domain[1][1] = localBox.max[1];
	domain[1][2] = localBox.max[2];

	_all->setVertices(domain);
	_all->setWork(work);

	std::vector<ALL::Point<double>> updatedVertices = _all->balance();
	
	return DomainBox(updatedVertices[0][0], updatedVertices[0][1], updatedVertices[0][2], updatedVertices[1][0], updatedVertices[1][1], updatedVertices[1][2]);
}
