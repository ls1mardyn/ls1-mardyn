/*
 * HaloRegion.h
 *
 *  Created on: Oct 13, 2016
 *      Author: seckler
 */

#pragma once

struct HaloRegion {
	double rmin[3]; // lower corner
	double rmax[3]; // higher corner
	int offset[3]; // offset (direction) of the halo region
	double width; // Halo width (e.g. one cutoff)
};
