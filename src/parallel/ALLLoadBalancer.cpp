/**
 * @file ALLLoadBalancer.cpp
 * @author seckler, Georg von Bismarck
 * @date 23.09.2026
 */

#include "ALLLoadBalancer.h"
#include <string>
#include "ALL.hpp"
#include "parallel/DomainDecompMPIBase.h"

ALLLoadBalancer::ALLLoadBalancer(double gamma, MPI_Comm comm, DomainGridPoint globalSize, std::vector<double> minimalPartitionSize) {
	_comm = comm;
	_gamma = gamma;
	_minimalPartitionSize = minimalPartitionSize;
	
	_coversWholeDomain = {globalSize[0] == 1, globalSize[1] == 1, globalSize[2] == 1};;
}

void ALLLoadBalancer::readXML(XMLfileUnits& xmlconfig){
	ALL::LB_t mode = ALL::LB_t::UNIMPLEMENTED;

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

	Log::global_log->info() << "ALLLoadBalancer: using the " << loadBalancer << " load balancer" << std::endl;

	_all = std::make_unique<ALL::ALL<double, double>>(ALL::TENSOR, DIMgeom, _gamma);
	_all->setCommunicator(_comm);
	_all->setMinDomainSize(_minimalPartitionSize);
    _all->setup();
}

DomainBox ALLLoadBalancer::rebalance(DomainBox localBox, double work) {
	std::vector<ALL::Point<double>> domain(2, ALL::Point<double>(DIMgeom));

	for (int i = 0; i < DIMgeom; ++i) {
		domain[0][i] = localBox[0][i];
		domain[1][i] = localBox[1][i];
	}

	_all->setVertices(domain);
	_all->setWork(work);
	_all->balance();

	std::vector<ALL::Point<double>> updatedVertices = _all->getVertices();
	DomainBox newlocalBox;

	for (int i = 0; i < DIMgeom; ++i) {
		newlocalBox[0][i] = updatedVertices[0][i];
		newlocalBox[1][i] = updatedVertices[1][i];
	}

	return newlocalBox;
}
