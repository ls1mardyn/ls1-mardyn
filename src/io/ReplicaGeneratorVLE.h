#ifndef REPLICAGENERATORVLE_H
#define REPLICAGENERATORVLE_H

#include "io/InputBase.h"
#include "molecules/Molecule.h"
#include "molecules/Component.h"
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <cstdint>

using namespace std;

/** @brief Generator of VLE scenario by replicating equilibrated liquid and vapor system.
 *
 * Description
 */
class DomainDecompBase;
class MoleculeDataReader;
class ReplicaGeneratorVLE : public InputBase
{
public:
	ReplicaGeneratorVLE();
	~ReplicaGeneratorVLE();

	void setPhaseSpaceFile(std::string /*filename*/){}
	void setPhaseSpaceHeaderFile(std::string /*filename*/){}

	void readPhaseSpaceHeader(Domain* /*domain*/, double /*timestep*/){}
	unsigned long readPhaseSpace(ParticleContainer* particleContainer, std::list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp);

	/** @brief Read in XML configuration for MkesferaGenerator and all its included objects.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <generator name="ReplicaGeneratorVLE">
	   </generator>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig);

private:
	void init();
	void readReplicaPhaseSpaceHeader(const std::string& strFilePathHeader, uint64_t& numParticles, double& dBoxLengthXYZ);
	void readReplicaPhaseSpaceData(const std::string& strFilePathData, const uint64_t& numParticles, std::vector<Molecule>& vecParticles);

private:
	std::vector<Molecule> _vecParticlesLiq;
	std::vector<Molecule> _vecParticlesVap;
	uint64_t _numParticlesLiq;
	uint64_t _numParticlesVap;
	uint64_t _numParticlesTotal;
	uint32_t _numBlocksXZ;
	uint32_t _numBlocksLiqY;
	uint32_t _numBlocksVapY;
	uint32_t _nIndexLiqBeginY;
	uint32_t _nIndexLiqEndY;
	std::string _strFilePathHeaderLiq;
	std::string _strFilePathDataLiq;
	std::string _strFilePathHeaderVap;
	std::string _strFilePathDataVap;
	double _dBoxLengthLiqXYZ;
	double _dBoxLengthVapXYZ;
	double _dBoxLengthXYZ;
	uint32_t _nMoleculeFormat;
	MoleculeDataReader* _moleculeDataReader;
	uint64_t _nMaxID;
	double _dMoleculeDiameter;
	double _fspY[6];  // free space positions
	double _dDensityLiq;
	double _dBoxVolumeLiq;
	bool _bCreateHomogenous;
};

class MoleculeDataReader
{
protected:
	MoleculeDataReader() {}
public:
	virtual ~MoleculeDataReader() {};
	virtual void read(std::ifstream& ifs, Molecule& mol, std::vector<Component>& components) = 0;
};

class MoleculeDataReaderICRVQD : public MoleculeDataReader
{
public:
	MoleculeDataReaderICRVQD() {}
public:
	virtual ~MoleculeDataReaderICRVQD() {};
	virtual void read(std::ifstream& ifs, Molecule& mol, std::vector<Component>& components)
	{
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

		Component* component = &components.at(cid-1);
		Molecule tmp(id, component,
				rx, ry ,rz,
				vx, vy, vz,
				q0, q1, q2, q3,
				Dx, Dy, Dz);
		mol = tmp;
	}
};

class MoleculeDataReaderICRV : public MoleculeDataReader
{
public:
	MoleculeDataReaderICRV() {}
public:
	virtual ~MoleculeDataReaderICRV() {};
	virtual void read(std::ifstream& ifs, Molecule& mol, std::vector<Component>& components)
	{
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

		Component* component = &components.at(cid-1);
		Molecule tmp(id, component,
				rx, ry ,rz,
				vx, vy, vz,
				q0, q1, q2, q3,
				Dx, Dy, Dz);
		mol = tmp;
	}
};

class MoleculeDataReaderIRV : public MoleculeDataReader
{
public:
	MoleculeDataReaderIRV() {}
public:
	virtual ~MoleculeDataReaderIRV() {};
	virtual void read(std::ifstream& ifs, Molecule& mol, std::vector<Component>& components)
	{
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
				rx, ry ,rz,
				vx, vy, vz,
				q0, q1, q2, q3,
				Dx, Dy, Dz);
		mol = tmp;
	}
};

#endif /* REPLICAGENERATORVLE_H */

