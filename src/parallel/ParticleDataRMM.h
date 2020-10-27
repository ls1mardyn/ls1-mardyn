#ifndef ParticleDataRMM_H_
#define ParticleDataRMM_H_

#include "mpi.h"

#include "molecules/MoleculeForwardDeclaration.h"
#include "particleContainer/adapter/vectorization/SIMD_TYPES.h"

//! @brief class to represent that particle data that is necessary for the exchange between processes
//! @author Martin Buchholz
//!
//! The data that has to be transfered for each particle consists of several values of different
//! types. As MPI-commands are expensive, it is desirable to transfer all data with a single
//! command. Therefore, two things have to be done:
//! - The data has to be "packed together", so all data at succeeding memory locations
//! - A MPI-datatype has to be defined.
//! This class achieves both. It is a class without constructor and destructor and only static
//! methods, which means that only the member variables use memory, so it can be used like a
//! C stuct. The static method setMPIType sets the given type to represent all the member variables.
class ParticleDataRMM {
public:
	//! @brief defines a MPI datatype which can be used to transfer a MacroscopicData object
	static void getMPIType(MPI_Datatype &sendPartType);

	//! @brief copy data from object of class Molecule to object of class ParticleDataRMM
	static void MoleculeToParticleData(ParticleDataRMM &particleStruct, Molecule &molecule);

	//! @brief copy data from object of class class ParticleDataRMM to object of class Molecule
	static void ParticleDataToMolecule(const ParticleDataRMM &particleStruct, Molecule &molecule);

	unsigned long id;
	vcp_real_calc r[3];
	vcp_real_calc v[3];
};

#endif /* ParticleDataRMM_H_ */
