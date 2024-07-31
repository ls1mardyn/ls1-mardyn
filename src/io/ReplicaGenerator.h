#ifndef REPLICA_GENERATOR_H
#define REPLICA_GENERATOR_H

#include "io/InputBase.h"
#include "molecules/Molecule.h"
#include "molecules/Component.h"
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <cstdint>
#include <memory>


enum SystemTypes : uint8_t {
	ST_UNKNOWN = 0,
	ST_HOMOGENEOUS = 1,
	ST_HETEROGENEOUS_VAPOR_LIQUID_VAPOR = 2,
	ST_HETEROGENEOUS_LIQUID_VAPOR = 3,
};

struct SubDomain {
	std::string strFilePathHeader;
	std::string strFilePathData;
	std::vector<Molecule> vecParticles;
	uint64_t numParticles;
	std::array<double, 3> arrBoxLength;
	double dVolume;
	double dDensity;
};

/** @brief Generator of VLE scenario by replicating equilibrated liquid and vapor system.
 *
 * Description
 */
class DomainDecompBase;

class MoleculeDataReader;

class ReplicaGenerator : public InputBase {
public:
	ReplicaGenerator();

	~ReplicaGenerator();

	void readPhaseSpaceHeader(Domain* /*domain*/, double /*timestep*/) {}

	unsigned long readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp);

	/** @brief Read in XML configuration for MkesferaGenerator and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <generator name="ReplicaGenerator">
	   </generator>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);

private:
	void init();

	void readReplicaPhaseSpaceHeader(SubDomain& subDomain);

	void readReplicaPhaseSpaceData(SubDomain& subDomain, DomainDecompBase* domainDecomp);

private:
	std::vector<SubDomain> _vecSubDomains;
	uint64_t _numParticlesTotal;
	uint32_t _numBlocksXZ;
	uint32_t _numBlocksLiqY;
	uint32_t _numBlocksVapY;
	uint32_t _nIndexLiqBeginY;
	uint32_t _nIndexLiqEndY;
	uint32_t _nMoleculeFormat;
	std::unique_ptr<MoleculeDataReader> _moleculeDataReader;
	double _dMoleculeDiameter;
	double _fspY[6];  // free space positions
	uint8_t _nSystemType;
	std::vector<uint32_t> _vecChangeCompIDsVap;
	std::vector<uint32_t> _vecChangeCompIDsLiq;
};

class MoleculeDataReader {
protected:
	MoleculeDataReader() {}

public:
	virtual ~MoleculeDataReader() {};

	virtual void read(std::ifstream& ifs, Molecule& mol, std::vector<Component>& components) = 0;
};

class MoleculeDataReaderICRVQD : public MoleculeDataReader {
public:
	MoleculeDataReaderICRVQD() {}

public:
	virtual ~MoleculeDataReaderICRVQD() {};

	virtual void read(std::ifstream& ifs, Molecule& mol, std::vector<Component>& components) {
		double rx, ry, rz, vx, vy, vz, q0, q1, q2, q3, Dx, Dy, Dz;
		rx = ry = rz = vx = vy = vz = q0 = q1 = q2 = q3 = Dx = Dy = Dz = 0.;
		uint64_t id = 0;
		uint32_t cid = 0;

		ifs.read(reinterpret_cast<char*> (&id), 8);
		ifs.read(reinterpret_cast<char*> (&cid), 4);
		ifs.read(reinterpret_cast<char*> (&rx), 8);
		ifs.read(reinterpret_cast<char*> (&ry), 8);
		ifs.read(reinterpret_cast<char*> (&rz), 8);
		ifs.read(reinterpret_cast<char*> (&vx), 8);
		ifs.read(reinterpret_cast<char*> (&vy), 8);
		ifs.read(reinterpret_cast<char*> (&vz), 8);
		ifs.read(reinterpret_cast<char*> (&q0), 8);
		ifs.read(reinterpret_cast<char*> (&q1), 8);
		ifs.read(reinterpret_cast<char*> (&q2), 8);
		ifs.read(reinterpret_cast<char*> (&q3), 8);
		ifs.read(reinterpret_cast<char*> (&Dx), 8);
		ifs.read(reinterpret_cast<char*> (&Dy), 8);
		ifs.read(reinterpret_cast<char*> (&Dz), 8);

		Component* component = &components.at(cid - 1);
		Molecule tmp(id, component,
					 rx, ry, rz,
					 vx, vy, vz,
					 q0, q1, q2, q3,
					 Dx, Dy, Dz);
		mol = tmp;
	}
};

class MoleculeDataReaderICRV : public MoleculeDataReader {
public:
	MoleculeDataReaderICRV() {}

public:
	virtual ~MoleculeDataReaderICRV() {};

	virtual void read(std::ifstream& ifs, Molecule& mol, std::vector<Component>& components) {
		double rx, ry, rz, vx, vy, vz, q0, q1, q2, q3, Dx, Dy, Dz;
		rx = ry = rz = vx = vy = vz = q0 = q1 = q2 = q3 = Dx = Dy = Dz = 0.;
		uint64_t id = 0;
		uint32_t cid = 0;

		ifs.read(reinterpret_cast<char*> (&id), 8);
		ifs.read(reinterpret_cast<char*> (&cid), 4);
		ifs.read(reinterpret_cast<char*> (&rx), 8);
		ifs.read(reinterpret_cast<char*> (&ry), 8);
		ifs.read(reinterpret_cast<char*> (&rz), 8);
		ifs.read(reinterpret_cast<char*> (&vx), 8);
		ifs.read(reinterpret_cast<char*> (&vy), 8);
		ifs.read(reinterpret_cast<char*> (&vz), 8);

		Component* component = &components.at(cid - 1);
		Molecule tmp(id, component,
					 rx, ry, rz,
					 vx, vy, vz,
					 q0, q1, q2, q3,
					 Dx, Dy, Dz);
		mol = tmp;
	}
};

class MoleculeDataReaderIRV : public MoleculeDataReader {
public:
	MoleculeDataReaderIRV() {}

public:
	virtual ~MoleculeDataReaderIRV() {};

	virtual void read(std::ifstream& ifs, Molecule& mol, std::vector<Component>& components) {
		double rx, ry, rz, vx, vy, vz, q0, q1, q2, q3, Dx, Dy, Dz;
		rx = ry = rz = vx = vy = vz = q0 = q1 = q2 = q3 = Dx = Dy = Dz = 0.;
		uint64_t id = 0;
		uint32_t cid = 0;

		ifs.read(reinterpret_cast<char*> (&id), 8);
		ifs.read(reinterpret_cast<char*> (&rx), 8);
		ifs.read(reinterpret_cast<char*> (&ry), 8);
		ifs.read(reinterpret_cast<char*> (&rz), 8);
		ifs.read(reinterpret_cast<char*> (&vx), 8);
		ifs.read(reinterpret_cast<char*> (&vy), 8);
		ifs.read(reinterpret_cast<char*> (&vz), 8);

		Component* component = &components.at(0);
		Molecule tmp(id, component,
					 rx, ry, rz,
					 vx, vy, vz,
					 q0, q1, q2, q3,
					 Dx, Dy, Dz);
		mol = tmp;
	}
};

#endif /* REPLICA_GENERATOR_H */

