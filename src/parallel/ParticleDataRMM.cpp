#include "ParticleDataRMM.h"

#include <mpi.h>
#include <typeinfo>

#include "ensemble/EnsembleBase.h"
#include "molecules/Molecule.h"
#include "Simulation.h"
#include "utils/Logger.h"


void ParticleDataRMM::getMPIType(MPI_Datatype &sendPartType) {
	int blocklengths[] = { 1, 6 }; // 1 unsLong value (id), 6 double values (3r, 3v)


	MPI_Datatype types[2];
	types[0] = MPI_UNSIGNED_LONG;


	ParticleDataRMM pdata_dummy;

	// ensure, that the types of v[0] and r[0] match!:
	mardyn_assert(typeid(pdata_dummy.v[0])==typeid(pdata_dummy.r[0]));

	// get the sizes of r and v - is it double or single precision?
	if (sizeof(pdata_dummy.r[0]) == 8) {  // 8 bytes for double
		types[1] = MPI_DOUBLE;
	} else if (sizeof(pdata_dummy.r[0]) == 4) {  // 4 bytes for single
		types[1] = MPI_FLOAT;
	} else {
		global_log->error() << "invalid size of vcp_real_calc";
		Simulation::exit(4852);
	}

	//if the following statement is not true, then the 6 double values do not follow one after the other.
	mardyn_assert(&(pdata_dummy.r[0]) + 3 == &(pdata_dummy.v[0]));

	MPI_Aint displacements[3];
#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 0
	MPI_CHECK( MPI_Get_address(&pdata_dummy.id, displacements) );
	MPI_CHECK( MPI_Get_address(&pdata_dummy.r[0], displacements + 1) );
#else
	MPI_CHECK( MPI_Address(&pdata_dummy.id, displacements) );
	MPI_CHECK( MPI_Address(&pdata_dummy.r[0], displacements + 1) );
#endif
	MPI_Aint base;
	MPI_CHECK( MPI_Get_address(&pdata_dummy, &base) );
	for (int i = 0; i < 2; i++)
		displacements[i] -= base;

#if MPI_VERSION >= 2 && MPI_SUBVERSION >= 0
	MPI_CHECK( MPI_Type_create_struct(2, blocklengths, displacements, types, &sendPartType) );
#else
	MPI_CHECK( MPI_Type_struct(2, blocklengths, displacements, types, &sendPartType) );
#endif
	MPI_CHECK( MPI_Type_commit(&sendPartType) );
}

void ParticleDataRMM::MoleculeToParticleData(ParticleDataRMM &particleStruct, Molecule &molecule) {
	particleStruct.id = molecule.getID();
	particleStruct.r[0] = molecule.r(0);
	particleStruct.r[1] = molecule.r(1);
	particleStruct.r[2] = molecule.r(2);
	particleStruct.v[0] = molecule.v(0);
	particleStruct.v[1] = molecule.v(1);
	particleStruct.v[2] = molecule.v(2);
}

void ParticleDataRMM::ParticleDataToMolecule(const ParticleDataRMM &particleStruct, Molecule &molecule) {
	Component* component = _simulation.getEnsemble()->getComponent(0);
	molecule = Molecule(particleStruct.id, component,
						particleStruct.r[0], particleStruct.r[1], particleStruct.r[2],
						particleStruct.v[0], particleStruct.v[1], particleStruct.v[2],
						1., 0., 0., 0.,
						0., 0., 0.
	);
}
