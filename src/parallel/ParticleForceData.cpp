#include "ParticleForceData.h"

#include <mpi.h>

#include "molecules/Molecule.h"
#include "Simulation.h"
#include "utils/Logger.h"


void ParticleForceData::getMPIType(MPI_Datatype &sendPartType) {
	int blocklengths[] = { 1, 12 }; // 1 unsLong value (id), 12 double values (3r, 3F, 3M, 3Vi)
	MPI_Datatype types[] = { MPI_UNSIGNED_LONG, MPI_DOUBLE };

	MPI_Aint displacements[2];
	ParticleForceData pdata_dummy;

	//if the following statements are not true, then the 12 double values do not follow one after the other!
	mardyn_assert(&(pdata_dummy.r[0]) + 3 == &(pdata_dummy.F[0]));
	mardyn_assert(&(pdata_dummy.r[0]) + 6 == &(pdata_dummy.M[0]));
	mardyn_assert(&(pdata_dummy.r[0]) + 9 == &(pdata_dummy.Vi[0]));

	MPI_CHECK( MPI_Get_address(&pdata_dummy.id, displacements) );
	MPI_CHECK( MPI_Get_address(&pdata_dummy.r[0], displacements + 1) );

	MPI_Aint base;
	MPI_CHECK( MPI_Get_address(&pdata_dummy, &base) );
	for (int i = 0; i < 2; i++){
		displacements[i] -= base;
	}

	MPI_CHECK( MPI_Type_create_struct(2, blocklengths, displacements, types, &sendPartType) );

	MPI_CHECK( MPI_Type_commit(&sendPartType) );
}

void ParticleForceData::MoleculeToParticleData(ParticleForceData &particleStruct, Molecule &molecule) {
	particleStruct.id = molecule.getID();
	particleStruct.r[0] = molecule.r(0);
	particleStruct.r[1] = molecule.r(1);
	particleStruct.r[2] = molecule.r(2);
	particleStruct.F[0] = molecule.F(0);
	particleStruct.F[1] = molecule.F(1);
	particleStruct.F[2] = molecule.F(2);
	particleStruct.M[0] = molecule.M(0);
	particleStruct.M[1] = molecule.M(1);
	particleStruct.M[2] = molecule.M(2);
	particleStruct.Vi[0] = molecule.Vi(0);
	particleStruct.Vi[1] = molecule.Vi(1);
	particleStruct.Vi[2] = molecule.Vi(2);
	particleStruct.Vi[3] = molecule.Vi(3);
	particleStruct.Vi[4] = molecule.Vi(4);
	particleStruct.Vi[5] = molecule.Vi(5);
	particleStruct.Vi[6] = molecule.Vi(6);
	particleStruct.Vi[7] = molecule.Vi(7);
	particleStruct.Vi[8] = molecule.Vi(8);
}

void ParticleForceData::AddParticleForceDataToMolecule(ParticleForceData &particleStruct, Molecule &molecule) {
		molecule.Fadd(particleStruct.F);
		molecule.Madd(particleStruct.M);
		molecule.ViaddAll(particleStruct.Vi);
}
