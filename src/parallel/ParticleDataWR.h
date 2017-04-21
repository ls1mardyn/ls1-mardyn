#ifndef ParticleDataWR_H_
#define ParticleDataWR_H_

#include <mpi.h>

#include "molecules/MoleculeForwardDeclaration.h"

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
class ParticleDataWR {
public:
	//! @brief defines a MPI datatype which can be used to transfer a MacroscopicData object
	static void getMPIType(MPI_Datatype &sendPartType);

	//! @brief copy data from object of class Molecule to object of class ParticleDataWR
	static void MoleculeToParticleData(ParticleDataWR &particleStruct, Molecule &molecule);

	//! @brief copy data from object of class class ParticleDataWR to object of class Molecule
	static void ParticleDataToMolecule(ParticleDataWR &particleStruct, Molecule &molecule);

	unsigned long id;
	double r[3];
	double v[3];
};

#endif /* ParticleDataWR_H_ */
