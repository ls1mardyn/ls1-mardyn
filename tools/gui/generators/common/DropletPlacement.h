/*
 * DropletGenerator.h
 *
 * @Date: 29.07.2010
 * @Author: Wolfgang Eckhardt, Martin Buchholz
 */

#ifndef DROPLETPLACEMENT_H_
#define DROPLETPLACEMENT_H_

#include "utils/Logger.h"
#include <vector>
#include <memory>

/**
 *  Places drops of different sizes in the domain, so that a given percentage
 *  of the volume (@see maxSphereVolume) is covered by fluid.
 */
class DropletPlacement {

public:

	/**
	 * A droplet is basically a sphere with center and radius.
	 */
	struct Droplet {
		Droplet(double c[3], double r) :
		 _radius(r){
			_center[0] = c[0];
			_center[1] = c[1];
			_center[2] = c[2];
		}

		double _center[3];
		double _radius;
	};

	/**
	 *  @param fluidVolume value in ]0;1[, percentage of volume covered by fluid, i.e. drops
	 *  @param maxSphereVolume value in ]0;1[, percentage of volume covered by the largest drop
	 *  @param numSphereSizes determines how many different sizes of spheres exist, with
	 *         - each class covering the same volume in total
	 *         - the size of a sphere is determined by pow(0.9, i) * maxSphereRadius; i in [1,maxSphereSize]
	 */
	DropletPlacement(double fluidVolume, double maxSphereVolume, int numSphereSizes, std::shared_ptr<Log::Logger> logger);

	/**
	 * Generates droplets with sizes as specified by numSphereSizes, fluidVolume
	 * and maxSphereVolume.
	 */
	std::vector<Droplet> generateDroplets();

	virtual ~DropletPlacement();

private:
	double _fluidVolume;
	double _maxSphereRadius;
	int _numSphereSizes;

	int _numOccupied;
	std::vector<std::vector<std::vector<bool> > > _occupiedFields;

	std::shared_ptr<Log::Logger> _logger;

	/**
	 * Initialize vector _occupiedFields
	 */
	void initFields(int size);

	/**
	 * calculate Euclidean distance
	 */
	double calcDistance(std::vector<double> pos1, double pos2[3]);

	/**
	 * Places a droplet randomly in the domain (spheres may overlap), adds it to droplets
	 * and marks occupiedFields to calculate the percentage of volume covered.
	 */
	void placeSphereRandomly(double radius, std::vector<Droplet>& droplets);
};

Log::Logger& operator<<(Log::Logger& str, DropletPlacement::Droplet& droplet);

#endif /* DropletPlacement_H_ */
