/*
 * DropletGenerator.cpp
 *
 * @Date: 29.07.2010
 * @Author: eckhardw
 */

#include "DropletPlacement.h"
#include <cmath>
#include <memory>

#include "utils/Logger.h"


#define SIZE 100

DropletPlacement::DropletPlacement(double fluidVolume, double maxSphereVolume, int numSphereSizes, std::shared_ptr<Log::Logger> logger)
 : _fluidVolume(fluidVolume / 100.), _maxSphereRadius(pow((3.0*maxSphereVolume / 100.)/(4.0 * M_PI), 1.0/3.0)),
   _numSphereSizes(numSphereSizes), _numOccupied(0), _logger(logger)
{
	initFields(SIZE);
	_logger -> debug() << "Instantiated Droplet Generator for parameter set " <<
			" fluidVolume " << _fluidVolume << " maxSphereVolume " << (maxSphereVolume/100) << " numSphereSizes " << _numSphereSizes << std::endl;
	_logger -> debug() << " _maxSphereRadius is " << _maxSphereRadius << std::endl;
}

DropletPlacement::~DropletPlacement() {
}

void DropletPlacement::initFields(int size){
	_occupiedFields.resize(size);
	for(int ix=0; ix<size; ix++){
		_occupiedFields[ix].resize(size);
		for(int iy=0; iy<size; iy++){
			_occupiedFields[ix][iy].resize(size);
			for(int iz=0; iz<size; iz++){
				_occupiedFields[ix][iy][iz] = false;
			}
		}
	}
}

double DropletPlacement::calcDistance(std::vector<double> pos1, double pos2[3]){
	double temp = 0;
	for(int i=0; i<3; i++){
		temp += pow(pos1[i]-pos2[i],2.0);
	}
	return sqrt(temp);
}


void DropletPlacement::placeSphereRandomly(double radius, std::vector<Droplet>& droplets){
	double center[3];
	center[0] = rand()/(1.0*RAND_MAX);
	center[1] = rand()/(1.0*RAND_MAX);
	center[2] = rand()/(1.0*RAND_MAX);

	//  int debug = 0;
	//  std::cout << "Debug " << debug << std::endl; debug++;
	std::vector<int> centerCell;
	centerCell.resize(3);
	centerCell[0] = (int) floor(center[0]*SIZE);
	centerCell[1] = (int) floor(center[1]*SIZE);
	centerCell[2] = (int) floor(center[2]*SIZE);

	int rInCells = (int) ceil(radius*SIZE);
	std::vector<double> cellPos;
	cellPos.resize(3);

	for(int ix=centerCell[0]-rInCells; ix<centerCell[0]+rInCells; ix++){
		for(int iy=centerCell[1]-rInCells; iy<centerCell[1]+rInCells; iy++){
			for(int iz=centerCell[2]-rInCells; iz<centerCell[2]+rInCells; iz++){
				// calculate position of the cell center
				cellPos[0] = (ix+0.5)/SIZE;
				cellPos[1] = (iy+0.5)/SIZE;
				cellPos[2] = (iz+0.5)/SIZE;
				// Distance of cell center to sphere center
				double distance = calcDistance(cellPos, center);
				if(distance<radius){
					if(_occupiedFields[(ix+SIZE)%SIZE][(iy+SIZE)%SIZE][(iz+SIZE)%SIZE]==false){
						_occupiedFields[(ix+SIZE)%SIZE][(iy+SIZE)%SIZE][(iz+SIZE)%SIZE]=true;
						_numOccupied++;
					}
				}
			}
		}
	}

	Droplet droplet(center, radius);
	droplets.push_back(droplet);
}

std::vector<DropletPlacement::Droplet> DropletPlacement::generateDroplets() {

	// in case there is only one droplet to be generated, place it in the middle!
	double maxDropVolume = (4.0 / 3.0 * M_PI * pow(_maxSphereRadius,3));
	if (_numSphereSizes == 1 && maxDropVolume >= _fluidVolume) {
		std::vector<Droplet> droplets;
		double center[3] = {0.5, 0.5, 0.5};
		Droplet droplet(center, _maxSphereRadius);
		droplets.push_back(droplet);
		return droplets;
	}

	// sphereclass 1:
	double currentRadius = _maxSphereRadius;
	double percentageOccupied = 0.0;
	double percentageNeeded = 0.0;
	std::vector<Droplet> droplets;
	for (int sphereclass = 0; sphereclass < _numSphereSizes; sphereclass++) {
		percentageNeeded = (sphereclass + 1.0) / _numSphereSizes * _fluidVolume;
		percentageOccupied = _numOccupied / pow(SIZE, 3.0);
		int count = 0;
		while (percentageOccupied < percentageNeeded) {
			placeSphereRandomly(currentRadius, droplets);
			count++;
			percentageOccupied = _numOccupied / pow(SIZE, 3.0);
		}
		_logger->debug() << "Created " << count << " spheres with radius " << currentRadius << std::endl;
		currentRadius -= _maxSphereRadius / 10.0;
	}
	_logger->debug() << "PercentOccupied: " << (100 * _numOccupied / pow(SIZE, 3.0)) << " %" << std::endl;
	return droplets;
}

Log::Logger& operator<<(Log::Logger& str, DropletPlacement::Droplet& droplet) {
	str << " Droplet: center [" << droplet._center[0] << "," << droplet._center[1] << "," << droplet._center[2] << "] r:" << droplet._radius << std::endl;
	return str;
}
