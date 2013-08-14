#include "parallel/ParticleData.h"

#include <mpi.h>

#include "ensemble/EnsembleBase.h"
#include "molecules/Molecule.h"
#include "Simulation.h"
#include "utils/Logger.h"


void ParticleData::setMPIType(MPI_Datatype &sendPartType) {
	int blocklengths[] = { 1, 1, 13 }; // 1 unsLong value (id), 1 int value (cid), 13 double values (3r, 3v, 4q, 3D)
	MPI_Datatype types[] = { MPI_UNSIGNED_LONG, MPI_INT, MPI_DOUBLE };

	MPI_Aint displacements[3];
	ParticleData pdata_dummy;
#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 0
	MPI_CHECK( MPI_Get_address(&pdata_dummy, displacements) );
	MPI_CHECK( MPI_Get_address(&pdata_dummy.cid, displacements + 1) );
	MPI_CHECK( MPI_Get_address(&pdata_dummy.r[0], displacements + 2) );
#else
	MPI_CHECK( MPI_Address(&pdata_dummy, displacements) );
	MPI_CHECK( MPI_Address(&pdata_dummy.cid, displacements + 1) );
	MPI_CHECK( MPI_Address(&pdata_dummy.r[0], displacements + 2) );
#endif
	MPI_Aint base = displacements[0];
	for (int i = 0; i < 3; i++)
		displacements[i] -= base;

#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 0
	MPI_CHECK( MPI_Type_create_struct(3, blocklengths, displacements, types, &sendPartType) );
#else
MPI_CHECK( MPI_Type_struct(3, blocklengths, displacements, types, &sendPartType) );
#endif
	MPI_CHECK( MPI_Type_commit(&sendPartType) );
}

void ParticleData::MoleculeToParticleData(ParticleData &particleStruct, Molecule &molecule) {
	particleStruct.id = molecule.id();
	particleStruct.cid = molecule.componentid();
	particleStruct.r[0] = molecule.r(0);
	particleStruct.r[1] = molecule.r(1);
	particleStruct.r[2] = molecule.r(2);
	particleStruct.v[0] = molecule.v(0);
	particleStruct.v[1] = molecule.v(1);
	particleStruct.v[2] = molecule.v(2);
	particleStruct.q[0] = molecule.q().qw();
	particleStruct.q[1] = molecule.q().qx();
	particleStruct.q[2] = molecule.q().qy();
	particleStruct.q[3] = molecule.q().qz();
	particleStruct.D[0] = molecule.D(0);
	particleStruct.D[1] = molecule.D(1);
	particleStruct.D[2] = molecule.D(2);
}

void ParticleData::ParticleDataToMolecule(ParticleData &particleStruct, Molecule **molecule) {
	Component* component = _simulation.getEnsemble()->component(particleStruct.cid);
	*molecule = new Molecule(particleStruct.id, component,
								particleStruct.r[0], particleStruct.r[1], particleStruct.r[2],
								particleStruct.v[0], particleStruct.v[1], particleStruct.v[2],
								particleStruct.q[0], particleStruct.q[1], particleStruct.q[2], particleStruct.q[3],
								particleStruct.D[0], particleStruct.D[1], particleStruct.D[2]
	);
}

#ifndef NDEBUG
ParticleData::ParticleData() : id(0), cid(-1) {
	for (int i = 0; i < 3; i++ ) {
		r[i] = 0.0;
		v[i] = 0.0;
		D[i] = 0.0;
		q[i] = 0.0;
	}
	q[3] = 0.0;
}
#endif
