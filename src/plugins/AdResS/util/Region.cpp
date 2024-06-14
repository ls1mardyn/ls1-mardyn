//
// Created by alex on 05.04.24.
//
#include "plugins/AdResS/util/Region.h"
#include "Simulation.h"
#include "Domain.h"

void Resolution::FPRegion::readXML(XMLfileUnits &xmlconfig) {
	xmlconfig.getNodeValue("lowX", _low[0]);
	xmlconfig.getNodeValue("lowY", _low[1]);
	xmlconfig.getNodeValue("lowZ", _low[2]);
	xmlconfig.getNodeValue("highX", _high[0]);
	xmlconfig.getNodeValue("highY", _high[1]);
	xmlconfig.getNodeValue("highZ", _high[2]);
	xmlconfig.getNodeValue("hybridDimX", _hybridDims[0]);
	xmlconfig.getNodeValue("hybridDimY", _hybridDims[1]);
	xmlconfig.getNodeValue("hybridDimZ", _hybridDims[2]);

	init();
}

void Resolution::FPRegion::init() {
	for(int d = 0; d < 3; d++) {
		_lowHybrid[d] = _low[d] - _hybridDims[d];
		_highHybrid[d] = _high[d] + _hybridDims[d];
		_center[d] = (_high[d] - _low[d])/2 + _low[d];
		_dim[d] = _high[d] - _low[d];
	}
}

std::array<double, 3>
Resolution::FPRegion::computeIntersection(std::array<double, 3> point, Resolution::FPRegion::Intersection inter) const {
	double scale = std::max(
			std::max(std::abs(point[0] - _center[0])/(_dim[0]+2*_hybridDims[0]*inter)*2.0,
					 std::abs(point[1] - _center[1])/(_dim[1]+2*_hybridDims[1]*inter)*2.0),
			std::abs(point[2] - _center[2])/(_dim[2]+2*_hybridDims[2]*inter)*2.0
	);
	std::array<double,3> result{0,0,0};
	for(int d = 0; d < 3; d++) result[d] = _center[d] + (point[d] - _center[d]) / scale;
	return result;
}

bool
Resolution::FPRegion::isInnerPoint(std::array<double, 3> point, std::array<double, 3> low, std::array<double, 3> high) {
	bool res = true;
	for(int d = 0; d < 3; d++) res &= point[d] >= low[d];
	for(int d = 0; d < 3; d++) res &= point[d] < high[d];
	return res;
}

bool Resolution::FPRegion::isRegionInBox(std::array<double, 3> low, std::array<double, 3> high) {
	return _low[0] <= high[0] && _low[1] <= high[1] && _low[2] <= high[2] &&
		   _high[0] >= low[0] && _high[1] >= low[1] && _high[2] >= low[2];
}
