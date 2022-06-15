#ifndef VELOCITYEXCHANGE_H_
#define VELOCITYEXCHANGE_H_

#include <string>
#include <map>
#include <unordered_map>
#include <list>
#include <cstdint>
#include <vector>
#include <array>
#include <string>
#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "plugins/PluginBase.h"
#include "molecules/MoleculeForwardDeclaration.h"

// struct pacIDtype // particle and component ID
// {
// 	uint64_t pid;  // particle ID
// 	uint32_t cid;  // component ID unity based
// };

class XMLfileUnits;
class Domain;
class ParticleContainer;
class DomainDecompBase;

class VelocityExchange: public PluginBase
{
private:

	struct TimestepControl
	{
		uint64_t start;
		uint64_t freq;
		uint64_t stop;
	};

public:
	// constructor and destructor
	VelocityExchange();
	~VelocityExchange();

  	/** @brief Read in XML configuration for VelocityExchange 
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
		<plugin name="VelocityExchange">
			<control>
				<start>0</start>
				<frequency>2</frequency>
				<stop>1000000000</stop>
			</control>
			<coldrange>  <!-- region with lower temperature -->
				<xmin>0</xmin> <xmax>100</xmax>  <!-- range x-axis -->
				<ymin>140</ymin> <ymax>160</ymax>  <!-- range y-axis -->
				<zmin>0</zmin> <zmax>100</zmax>  <!-- range z-axis -->
			</coldrange>
			<warmrange>  <!-- region with higher temperature -->
				<symmetric>1</symmetric>  <!-- 0: no symmetry; 1: symmetry in y direction -->
				<xmin>0</xmin> <xmax>100</xmax>  <!-- range x-axis -->
				<ymin>20</ymin> <ymax>30</ymax>  <!-- range y-axis -->
				<zmin>0</zmin> <zmax>100</zmax>  <!-- range z-axis -->
			</warmrange>
		</plugin>
	   \endcode
	 */

	void readXML(XMLfileUnits& xmlconfig) override;

	void init(ParticleContainer *particleContainer,
			  DomainDecompBase *domainDecomp, Domain *domain) override;

    /** @brief Method beforeForces will be called before forcefields have been applied
     * no alterations w.r.t. Forces shall be made here
     *
     */
    void beforeForces(
            ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */,
            unsigned long /* simstep */
    ) override;

    virtual void endStep(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            Domain* domain, unsigned long simstep) override {}

	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain) override {}

	std::string getPluginName() override {return std::string("VelocityExchange");}
	static PluginBase* createInstance() {return new VelocityExchange();}

private:
	void exchangeVelocities(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep);

private:
	TimestepControl _control;
  struct C_Range {double xmin, xmax, ymin, ymax, zmin, zmax;} _cold_range;
  struct W_Range {double xmin, xmax, ymin, ymax, zmin, zmax;} _warm_range;

	template<typename T>
	void resizeExactly(std::vector<T>& v, unsigned int numElements) const {
		v.reserve(numElements);
		v.resize(numElements);
	}

  bool symmetry;
  double _boxLength[3];

  std::vector<double>  v_c_abs_local;
  std::vector<double>  v_w_abs_local;
  std::vector<double>  v_c_abs_global;
  std::vector<double>  v_w_abs_global;

  std::vector<uint32_t>  cold_mol_local;
  std::vector<uint32_t>  warm_mol_local;
  std::vector<uint32_t>  cold_mol_global;
  std::vector<uint32_t>  warm_mol_global;

  std::vector<double>  v_cold_local;
  std::vector<double>  v_warm_local;
  std::vector<double>  D_cold_local;
  std::vector<double>  D_warm_local;

  std::vector<double>  v_cold_global;
  std::vector<double>  v_warm_global;
  std::vector<double>  D_cold_global;
  std::vector<double>  D_warm_global;

};

#endif /*VELOCITYEXCHANGE_H_*/