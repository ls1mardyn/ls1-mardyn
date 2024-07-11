//
// Created by alex on 05.04.24.
//
#include "plugins/PotentialOfMeanForce/Region.h"
#include "Simulation.h"

void FPRegion::readXML(XMLfileUnits &xmlconfig) {
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

void FPRegion::init() {
	for(int d = 0; d < 3; d++) {
		_lowHybrid[d] = _low[d] - _hybridDims[d];
		_highHybrid[d] = _high[d] + _hybridDims[d];
		_center[d] = (_high[d] - _low[d])/2 + _low[d];
		_dim[d] = _high[d] - _low[d];
	}
}

std::array<double, 3>
FPRegion::computeIntersection(std::array<double, 3> point, FPRegion::Intersection inter) const {
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
FPRegion::isInnerPoint(std::array<double, 3> point, std::array<double, 3> low, std::array<double, 3> high) {
	bool res = true;
	for(int d = 0; d < 3; d++) res &= point[d] >= low[d];
	for(int d = 0; d < 3; d++) res &= point[d] < high[d];
	return res;
}

bool FPRegion::isInnerPointDomain(Domain *domain, ResolutionType region,
											  std::array<double, 3> point) const {
	std::array<double,3> globLen{0};
	for(int d = 0; d < 3; d++) {
		globLen[d] = domain->getGlobalLength(d);
	}

	std::array<bool, 3> inDim {false};
	bool inDomain = true;
	for(int d = 0; d < 3; d++) {
		inDim[d] = _lowHybrid[d] >= 0 && _highHybrid[d] <= globLen[d];
		inDomain &= inDim[d];
	}

	//special case for interface situation
	std::array<double,3> borderOffset{0,0,0};
	for(int d = 0; d < 3; d++) {
		borderOffset[d] = (_hybridDims[d] == 0) ? _simulation.getcutoffRadius() : 0;
	}

	//region is within domain bounds
	if(inDomain && region == FullParticle) return isInnerPoint(point, {_low[0] - borderOffset[0], _low[1] - borderOffset[1], _low[2] - borderOffset[2]}, {_high[0] + borderOffset[0], _high[1] + borderOffset[1], _high[2] + borderOffset[2]});
	if(inDomain && region == Hybrid) return isInnerPoint(point, {_lowHybrid[0] - borderOffset[0], _lowHybrid[1] - borderOffset[1], _lowHybrid[2] - borderOffset[2]}, {_highHybrid[0] + borderOffset[0], _highHybrid[1] + borderOffset[1], _highHybrid[2] + borderOffset[2]});

	//region crosses bound
	bool checkOuter = region == Hybrid;
	std::array<std::array<double,2>,3> deltaLowHighDim{};
	for(int d = 0; d < 3; d++) {
		deltaLowHighDim[d][0] = std::max(-(_low[d] - checkOuter * _hybridDims[d]), 0.);
		deltaLowHighDim[d][1] = std::max((_high[d] + checkOuter * _hybridDims[d]) - globLen[d], 0.);
	}

	// check all regions that wrap around due to periodic bounds
	// we only need to check for each dim twice if it actually wraps around in that dimension
	// in the arrays we create the indices of the boxes depending on the viewed case
	// for x y or z: if they are 0 then we do not check a wrap around in that dimension
	bool result = false;
	for(int x = 0; x <= !inDim[0]; x++) {
		for(int y = 0; y <= !inDim[1]; y++) {
			for(int z = 0; z <= !inDim[2]; z++) {
				std::array<double, 3> checkLow {
						(1 - x) * std::max(_low[0] - checkOuter * _hybridDims[0], 0.) + x * ((deltaLowHighDim[0][0] != 0) * (globLen[0] - deltaLowHighDim[0][0]) + (deltaLowHighDim[0][1] != 0) * 0),
						(1 - y) * std::max(_low[1] - checkOuter * _hybridDims[1], 0.) + y * ((deltaLowHighDim[1][0] != 0) * (globLen[1] - deltaLowHighDim[1][0]) + (deltaLowHighDim[1][1] != 0) * 0),
						(1 - z) * std::max(_low[2] - checkOuter * _hybridDims[2], 0.) + z * ((deltaLowHighDim[2][0] != 0) * (globLen[2] - deltaLowHighDim[2][0]) + (deltaLowHighDim[2][1] != 0) * 0) };

				std::array<double, 3> checkHigh {
						(1 - x) * std::min(_high[0] + checkOuter * _hybridDims[0], globLen[0]) + x * ((deltaLowHighDim[0][0] != 0) * (globLen[0]) + (deltaLowHighDim[0][1] != 0) * deltaLowHighDim[0][1]),
						(1 - y) * std::min(_high[1] + checkOuter * _hybridDims[1], globLen[1]) + y * ((deltaLowHighDim[1][0] != 0) * (globLen[1]) + (deltaLowHighDim[1][1] != 0) * deltaLowHighDim[1][1]),
						(1 - z) * std::min(_high[2] + checkOuter * _hybridDims[2], globLen[2]) + z * ((deltaLowHighDim[2][0] != 0) * (globLen[2]) + (deltaLowHighDim[2][1] != 0) * deltaLowHighDim[2][1]) };
				result |= isInnerPoint(point, checkLow, checkHigh);
			}
		}
	}
	return result;
}

bool FPRegion::isRegionInBox(std::array<double, 3> low, std::array<double, 3> high) {
	return _low[0] <= high[0] && _low[1] <= high[1] && _low[2] <= high[2] &&
		   _high[0] >= low[0] && _high[1] >= low[1] && _high[2] >= low[2];
}

bool FPRegion::IsInsideResolutionRegion(std::array<double,3> point, ResolutionType resolution){

	if(resolution == ResolutionType::FullParticle){
		std::array<double,3> shifted_low = this->_low;
		std::array<double,3> shifted_high = this->_high;

		shifted_low[1]-= _simulation.getcutoffRadius();
		shifted_low[2]-= _simulation.getcutoffRadius();


		shifted_high[1]+= _simulation.getcutoffRadius();
		shifted_high[2]+= _simulation.getcutoffRadius();

		if(_low[0]==0.0){
			shifted_low[0] -= _simulation.getcutoffRadius();
		}


		if(_high[0]==_simulation.getDomain()->getGlobalLength(0)){
			shifted_high[0] += _simulation.getcutoffRadius();
		}

		return isInnerPoint(point, shifted_low, shifted_low);
	}

	if(resolution == ResolutionType::Hybrid){

		std::array<double,3> shifted_low = this->_lowHybrid;
		std::array<double,3> shifted_high = this->_highHybrid;

		shifted_low[1]-= _simulation.getcutoffRadius();
		shifted_low[2]-= _simulation.getcutoffRadius();


		shifted_high[1]+= _simulation.getcutoffRadius();
		shifted_high[2]+= _simulation.getcutoffRadius();

			//TODO:properly shift these coordinates on the x axis
			return isInnerPoint(point, shifted_low, shifted_high);
		
	}

	if(resolution == ResolutionType::Hybrid){
		std::array<double,3> shifted_low {0.0,0.0,0.0};
		std::array<double,3> shifted_high {_simulation.getDomain()->getGlobalLength(0),_simulation.getDomain()->getGlobalLength(1),_simulation.getDomain()->getGlobalLength(2)};

		shifted_low[1]-= _simulation.getcutoffRadius();
		shifted_low[2]-= _simulation.getcutoffRadius();


		shifted_high[1]+= _simulation.getcutoffRadius();
		shifted_high[2]+= _simulation.getcutoffRadius();

		if(!_low[0]==0.0){
			shifted_low[0] -= _simulation.getcutoffRadius();
		}


		if(!_high[0]==_simulation.getDomain()->getGlobalLength(0)){
			shifted_high[0] += _simulation.getcutoffRadius();
		}

		return isInnerPoint(point, shifted_low, shifted_low);
	}
	Log::global_log->info()<<"We should not reach this point"<<std::endl;
	Simulation::exit(670);
	
}
