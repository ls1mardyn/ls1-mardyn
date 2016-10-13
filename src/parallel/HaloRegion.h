/*
 * HaloRegion.h
 *
 *  Created on: Oct 13, 2016
 *      Author: seckler
 */

#ifndef SRC_PARALLEL_HALOREGION_H_
#define SRC_PARALLEL_HALOREGION_H_

struct HaloRegion {
	double rmin[3]; // lower corner
	double rmax[3]; // higher corner
	int offset[3]; // offset (direction) of the halo region
};

#endif /* SRC_PARALLEL_HALOREGION_H_ */
