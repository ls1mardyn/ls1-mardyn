#ifndef ParticleForceData_H_
#define ParticleForceData_H_

#include "mpi.h"

#include "molecules/MoleculeForwardDeclaration.h"

//! @brief Class to represent the particle force data that is necessary for the exchange between processes
class ParticleForceData {
public:
	//! @brief defines a MPI datatype which can be used to transfer a MacroscopicData object
	static void getMPIType(MPI_Datatype &sendPartType);

	//! @brief copy data from object of class Molecule to object of class ParticleForceData
	static void MoleculeToParticleData(ParticleForceData &particleStruct, Molecule &molecule);

	//! @brief add data from object of class ParticleForceData to object of class Molecule
	static void AddParticleForceDataToMolecule(ParticleForceData &particleStruct, Molecule &molecule);

	unsigned long id; //< ID of the particle
	int cid;  //< ID of the component of the particle
	double r[3];  //< position
	double F[3]; //< force
	double M[3]; //< torsional moment
	double Vi[3]; //< virial tensor
};

#endif /* ParticleForceData_H_ */
