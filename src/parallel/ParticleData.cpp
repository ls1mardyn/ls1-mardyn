#include "parallel/ParticleData.h"

#include <mpi.h>

#include "ensemble/EnsembleBase.h"
#include "molecules/Molecule.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "parallel/CollectiveCommunication.h"

void ParticleData::setMPIType(MPI_Datatype &sendPartType) {
	int blocklengths[] = { 2, 1, 58 }; // 2 unsLong values (id, cluster), 1 int value (cid), 13 double values (3r, 3v, 4q, 3D)
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

void ParticleData::MoleculeToParticleData(ParticleData &particleStruct, Molecule &molecule)
{
   MoleculeToParticleData(particleStruct, molecule, 0);
}
void ParticleData::MoleculeToParticleData(ParticleData &particleStruct, Molecule &molecule, unsigned long cluster_id)
{
	particleStruct.id = molecule.id();
        particleStruct.cluster = cluster_id;
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
	particleStruct.rOld[0] = molecule.rOld(0);
	particleStruct.rOld[1] = molecule.rOld(1);
	particleStruct.rOld[2] = molecule.rOld(2);
	particleStruct.rOldStress[0] = molecule.rOldStress(0);
	particleStruct.rOldStress[1] = molecule.rOldStress(1);
	particleStruct.rOldStress[2] = molecule.rOldStress(2);
	particleStruct.vdir[0] = (double)molecule.getDirectedVelocity(0);
	particleStruct.vdir[1] = (double)molecule.getDirectedVelocity(1);
	particleStruct.vdir[2] = (double)molecule.getDirectedVelocity(2);
	particleStruct.vdirSlab[0] = (double)molecule.getDirectedVelocitySlab(0);
	particleStruct.vdirSlab[1] = (double)molecule.getDirectedVelocitySlab(1);
	particleStruct.vdirSlab[2] = (double)molecule.getDirectedVelocitySlab(2);
	particleStruct.vdirStress[0] = (double)molecule.getDirectedVelocityStress(0);
	particleStruct.vdirStress[1] = (double)molecule.getDirectedVelocityStress(1);
	particleStruct.vdirStress[2] = (double)molecule.getDirectedVelocityStress(2);
	particleStruct.vdirConfinement[0] = (double)molecule.getDirectedVelocityConfinement(0);
	particleStruct.vdirConfinement[1] = (double)molecule.getDirectedVelocityConfinement(1);
	particleStruct.vdirConfinement[2] = (double)molecule.getDirectedVelocityConfinement(2);
	particleStruct.vdirAvConf[0] = (double)molecule.getDirectedVelocityAverageConfinement(0);
	particleStruct.vdirAvConf[1] = (double)molecule.getDirectedVelocityAverageConfinement(1);
	particleStruct.vdirAvConf[2] = (double)molecule.getDirectedVelocityAverageConfinement(2);
	particleStruct.vdirAvStress[0] = (double)molecule.getDirectedVelocityAverageStress(0);
	particleStruct.vdirAvStress[1] = (double)molecule.getDirectedVelocityAverageStress(1);
	particleStruct.vdirAvStress[2] = (double)molecule.getDirectedVelocityAverageStress(2);
	particleStruct.pressVir[0] = (double)molecule.getPressureVirial(0);
	particleStruct.pressVir[1] = (double)molecule.getPressureVirial(1);
	particleStruct.pressVir[2] = (double)molecule.getPressureVirial(2);
	particleStruct.pressKin[0] = (double)molecule.getPressureKin(0);
	particleStruct.pressKin[1] = (double)molecule.getPressureKin(1);
	particleStruct.pressKin[2] = (double)molecule.getPressureKin(2);
	particleStruct.pressVirConf[0] = (double)molecule.getPressureVirialConfinement(0);
	particleStruct.pressVirConf[1] = (double)molecule.getPressureVirialConfinement(1);
	particleStruct.pressVirConf[2] = (double)molecule.getPressureVirialConfinement(2);
	particleStruct.pressKinConf[0] = (double)molecule.getPressureKinConfinement(0);
	particleStruct.pressKinConf[1] = (double)molecule.getPressureKinConfinement(1);
	particleStruct.pressKinConf[2] = (double)molecule.getPressureKinConfinement(2);
	particleStruct.F[0] = (double)molecule.F(0);
	particleStruct.F[1] = (double)molecule.F(1);
	particleStruct.F[2] = (double)molecule.F(2);
	particleStruct.conPotHeat[0] = (double)molecule.getConvectivePotHeatflux(0);
	particleStruct.conPotHeat[1] = (double)molecule.getConvectivePotHeatflux(1);
	particleStruct.conPotHeat[2] = (double)molecule.getConvectivePotHeatflux(2);
	particleStruct.conPotHeatStress[0] = (double)molecule.getConvectivePotHeatfluxStress(0);
	particleStruct.conPotHeatStress[1] = (double)molecule.getConvectivePotHeatfluxStress(1);
	particleStruct.conPotHeatStress[2] = (double)molecule.getConvectivePotHeatfluxStress(2);
	
}

unsigned long ParticleData::ParticleDataToMolecule(ParticleData &particleStruct, Molecule **molecule) {
	Component* component = _simulation.getEnsemble()->component(particleStruct.cid);
	*molecule = new Molecule(particleStruct.id, component,
								particleStruct.r[0], particleStruct.r[1], particleStruct.r[2],
								particleStruct.v[0], particleStruct.v[1], particleStruct.v[2],
								particleStruct.q[0], particleStruct.q[1], particleStruct.q[2], particleStruct.q[3],
								particleStruct.D[0], particleStruct.D[1], particleStruct.D[2],
								particleStruct.rOld[0], particleStruct.rOld[1], particleStruct.rOld[2],
								particleStruct.rOldStress[0], particleStruct.rOldStress[1], particleStruct.rOldStress[2],
								particleStruct.vdir[0], particleStruct.vdir[1], particleStruct.vdir[2],
								particleStruct.vdirSlab[0], particleStruct.vdirSlab[1], particleStruct.vdirSlab[2],
								particleStruct.vdirStress[0], particleStruct.vdirStress[1], particleStruct.vdirStress[2],
								particleStruct.vdirConfinement[0], particleStruct.vdirConfinement[1], particleStruct.vdirConfinement[2],
								particleStruct.vdirAvConf[0], particleStruct.vdirAvConf[1], particleStruct.vdirAvConf[2],
								particleStruct.vdirAvStress[0], particleStruct.vdirAvStress[1], particleStruct.vdirAvStress[2],
								particleStruct.pressVir[0], particleStruct.pressVir[1], particleStruct.pressVir[2],
								particleStruct.pressKin[0], particleStruct.pressKin[1], particleStruct.pressKin[2],
								particleStruct.pressVirConf[0], particleStruct.pressVirConf[1], particleStruct.pressVirConf[2],
								particleStruct.pressKinConf[0], particleStruct.pressKinConf[1], particleStruct.pressKinConf[2],
								particleStruct.F[0], particleStruct.F[1], particleStruct.F[2],
								particleStruct.conPotHeat[0], particleStruct.conPotHeat[1], particleStruct.conPotHeat[2],
								particleStruct.conPotHeatStress[0], particleStruct.conPotHeatStress[1], particleStruct.conPotHeatStress[2]
	);
        return particleStruct.cluster;
}

#ifndef NDEBUG
ParticleData::ParticleData()
{
        id = 0;
        cluster = 0;
        cid = -1;
	for (int i = 0; i < 3; i++ ) {
		r[i] = 0.0;
		v[i] = 0.0;
		D[i] = 0.0;
		q[i] = 0.0;
		rOld[i] = 0.0;
		rOldStress[i] = 0.0;
		vdir[i] = 0.0;
		vdirSlab[i] = 0.0;
		vdirStress[i] = 0.0;
		vdirConfinement[i] = 0.0;
		vdirAvConf[i] = 0.0;
		vdirAvStress[i] = 0.0;
		pressVir[i] = 0.0;
		pressKin[i] = 0.0;
		pressVirConf[i] = 0.0;
		pressKinConf[i] = 0.0;
		F[i] = 0.0;
		conPotHeat[i] = 0.0;
		conPotHeatStress[i] = 0.0;
	}
	q[3] = 0.0;
}
#endif
