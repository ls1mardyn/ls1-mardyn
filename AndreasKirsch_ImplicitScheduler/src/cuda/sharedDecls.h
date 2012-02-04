/*
 * sharedDecls.h
 *
 *  Created on: Jun 21, 2011
 *      Author: andreas
 */

#ifndef SHAREDDECLS_H_
#define SHAREDDECLS_H_

#include <host_defines.h>

#include "config.h"

struct LockStorage {
	volatile uint lockValue;
};

typedef char Molecule_ComponentType;

struct CellStatsStorage {
	volatile floatType potential;
	volatile floatType virial;
};

struct Matrix3x3Storage {
	floatType3 rows[3];
};

struct QuaternionStorage {
	floatType w, x, y, z;
};

struct LJParameters {
	floatType epsilon;
	floatType sigma;
};

struct ComponentDescriptor {
	int numLJCenters;
	int numCharges;
	int numDipoles;
	int numQuadrupoles;

	// double3 has a 4byte alignment in gcc x86 and an 8 byte alignment in nvcc
	// it is defined in vector_types.h and its alignment can't be changed easily.
	// FIXME: manual padding required!

	struct Site {
		floatType3 relativePosition;
	};

	struct LJCenter : Site {
		LJParameters ljParameters;
	};

	struct OrientedSite : Site {
		floatType3 relativeE;
	};

	struct Dipole : OrientedSite {
		floatType absMu;
	};

	struct Charge : Site {
		floatType q;
	};

	struct Quadrupole : OrientedSite {
		floatType absQ;
	};

#if MAX_NUM_LJCENTERS > 0
	LJCenter ljCenters[MAX_NUM_LJCENTERS];
#endif

#if MAX_NUM_CHARGES > 0
	Charge charges[MAX_NUM_CHARGES];
#endif

#if MAX_NUM_DIPOLES > 0
	Dipole dipoles[MAX_NUM_DIPOLES];
#endif

#if MAX_NUM_QUADRUPOLES > 0
	Quadrupole quadrupoles[MAX_NUM_QUADRUPOLES];
#endif
};

typedef ComponentDescriptor ComponentDescriptors[MAX_NUM_COMPONENTS];
typedef floatType ComponentMixCoefficients[MAX_NUM_COMPONENTS][MAX_NUM_COMPONENTS];


#endif /* SHAREDDECLS_H_ */
