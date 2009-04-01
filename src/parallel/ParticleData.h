#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "molecules/Molecule.h"

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
class ParticleData{
 public:
  //! @brief defines a MPI datatype which can be used to transfer a MacroscopicData object
  static void setMPIType(MPI_Datatype &sendPartType){
    int blocklengths[] = {1,1,13}; // 1 unsLong value (id), 1 int value (cid), 13 double values (3r, 3v, 4q, 3D)
    MPI_Aint displacements[] = {0, sizeof(unsigned long),sizeof(int)+sizeof(unsigned long)};
    MPI_Datatype types[] = {MPI_UNSIGNED_LONG, MPI_INT,MPI_DOUBLE};
    MPI_Type_create_struct(3, blocklengths, displacements, types, &sendPartType);
  	MPI_Type_commit(&sendPartType);
  }

  //! @brief copy data from object of class Molecule to object of class ParticleData
  static void setParticleData(ParticleData &particleStruct, Molecule &molecule){
    particleStruct.id = molecule.id();
    particleStruct.cid= molecule.componentid();
    particleStruct.rx = molecule.r(0);
    particleStruct.ry = molecule.r(1);
    particleStruct.rz = molecule.r(2);
    particleStruct.vx = molecule.v(0);
    particleStruct.vy = molecule.v(1);
    particleStruct.vz = molecule.v(2);
    particleStruct.qw = molecule.q().qw();
    particleStruct.qx = molecule.q().qx();
    particleStruct.qy = molecule.q().qy();
    particleStruct.qz = molecule.q().qz();
    particleStruct.Dx = molecule.D(0);
    particleStruct.Dy = molecule.D(1);
    particleStruct.Dz = molecule.D(2);
  }

  unsigned long id;
  int cid;
  double rx, ry, rz;
  double vx, vy, vz;
  double qw, qx, qy, qz;
  double Dx, Dy, Dz;
};

#endif /*PARTICLE_H_*/
