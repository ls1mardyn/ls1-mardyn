/*
 * DensityControl.h
 *
 *  Created on: 16.11.2021
 *      Author: mheinen
 */

#ifndef DENSITYCONTROL_H_
#define DENSITYCONTROL_H_

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
#include "utils/Random.h"

struct pacIDtype // particle and component ID
{
	uint64_t pid;  // particle ID
	uint32_t cid;  // component ID unity based
};

class XMLfileUnits;
class Domain;
class ParticleContainer;
class DomainDecompBase;

class DensityControl : public PluginBase
{
private:

	struct TimestepControl
	{
		uint64_t start;
		uint64_t freq;
		uint64_t stop;
	};
	
	struct DensityTargetType
	{
		std::vector<double> density;
		std::vector<uint64_t> count;
	};

public:
	// constructor and destructor
	DensityControl();
	~DensityControl();

	/** @brief Read in XML configuration for DensityControl and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
		<plugin name="DensityControl">
			<control>
				<start>INT</start>           <!-- start time step -->
				<frequency>INT</frequency>   <!-- frequency of checking -->
				<stop>INT</stop>             <!-- stop time step -->
			</control>
			<range>
				<inclusive>BOOL</inclusive>  <!-- enable checking inside range (default): 1 | outside range: 0 -->
				<xmin>FLOAT</xmin> <xmax>FLOAT</xmax>  <!-- range x-axis -->
				<ymin>FLOAT</ymin> <ymax>FLOAT</ymax>  <!-- range y-axis -->
				<zmin>FLOAT</zmin> <zmax>FLOAT</zmax>  <!-- range z-axis -->
			</range>
			<targets>
				<target cid="INT">   <!-- cid: component id of target particles -->
					<density>DOUBLE</density>            <!-- target density -->
				</target>
			</targets>
			<priority>INT,INT,INT</priority>  <!-- Usally: Molecule size sorted in descending order -->
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

	std::string getPluginName() override {return std::string("DensityControl");}
	static PluginBase* createInstance() {return new DensityControl();}

private:
	bool moleculeInsideRange(std::array<double,3>& r);
	void controlDensity(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep);
	void updateBalanceVector(std::vector<int64_t>& vecBalance, std::map<int, std::vector<uint64_t> >& pidMap, uint32_t& numComponents);
	uint32_t tokenize_int_list(std::vector<uint32_t>& vec, std::string str, std::string del = ",");

private:
	Random _rnd;
	TimestepControl _control;
	DensityTargetType _densityTarget;
	struct Range {double xmin, xmax, ymin, ymax, zmin, zmax, volume;} _range;
	std::vector<uint32_t> _vecPriority;
#ifdef ENABLE_MPI
	MPI_Datatype pacID_mpi_type;
#endif
};

#endif /*DENSITYCONTROL_H_*/
