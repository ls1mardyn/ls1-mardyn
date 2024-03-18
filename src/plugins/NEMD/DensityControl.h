/*
 * DensityControl.h
 *
 *  Created on: 16.11.2021
 *	  Author: mheinen
 */

#pragma once

class DensityControlTest;

#include <array>
#include <cstdint>
#include <list>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>
#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "molecules/MoleculeForwardDeclaration.h"
#include "particleContainer/ParticleContainer.h"
#include "plugins/PluginBase.h"
#include "utils/CommVar.h"
#include "utils/Random.h"

class XMLfileUnits;
class Domain;
class ParticleContainer;
class DomainDecompBase;

class DensityControl : public PluginBase {
private:
	friend DensityControlTest;

	// particle and component ID
	struct pacIDtype {
		uint64_t pid;  // particle ID
		uint32_t cid;  // component ID unity based
	};

	struct TimestepControl {
		uint64_t start;
		uint64_t freq;
		uint64_t stop;
	};

	struct DensityTargetType {
		std::vector<double> density;
		std::vector<uint64_t> count;
		std::vector<bool> specified;
	};

public:
	// constructor and destructor
	DensityControl();
	~DensityControl() override;

	/** @brief Read in XML configuration for DensityControl and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
		<plugin name="DensityControl">
			<control>
				<start>INT</start>            <!-- start time step -->
				<frequency>INT</frequency>    <!-- frequency of checking -->
				<stop>INT</stop>              <!-- stop time step -->
			</control>
			<range>
				<inclusive>BOOL</inclusive>  <!-- enable checking inside range (default): true | outside range: false -->
				<xmin>FLOAT</xmin> <xmax>FLOAT</xmax>  <!-- range x-axis -->
				<ymin>FLOAT</ymin> <ymax>FLOAT</ymax>  <!-- range y-axis -->
				<zmin>FLOAT</zmin> <zmax>FLOAT</zmax>  <!-- range z-axis -->
			</range>
			<targets>
				<target cid="INT">   <!-- cid: component id of target particles -->
					<density>DOUBLE</density>          <!-- target density -->
				</target>
			</targets>
			<priority>INT,INT,INT</priority>  <!-- Comma separated list of component IDs (unity based) to specify the
	 priority for component exchange. \ Usually, IDs should be sorted in descending order of corresponding molecule size
	 -->
		</plugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;

	/**
	 * Method beforeForces will be called before forcefields have been applied
	 * no alterations w.r.t. Forces shall be made here
	 */
	void beforeForces(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */,
					  unsigned long /* simstep */
					  ) override;

	void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
				 unsigned long simstep) override {}

	void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override {}

	std::string getPluginName() override { return std::string("DensityControl"); }
	static PluginBase* createInstance() { return new DensityControl(); }

private:
	/**
	 * Uses the position \p r of a molecule to check whether it is in the specified range \p _range in which the
	 * density has to be controlled by this plugin.
	 * @param r Position of a molecule
	 * @ return True, if molecule is located inside the specified range \p _range,		otherwise False.
	 */
	bool moleculeInsideRange(const std::array<double, 3>& r);
	/**
	 * Establishes the specified target partial densities of all components, stored in \p _densityTarget . The
	 * target total density is calculated by the sum of all target partial densities. If the target value for a
	 * component is not specified, it is assumed that the target value equals the initial value for this component at
	 *the specified start timestep \p _control.start , with which the plugin starts to work. First, an attempt is made
	 *to produce the desired composition by component exchange. Therefore, the order (priority) of the component IDs,
	 *which is stored in \p _vecPriority , is taken into account. Usually, the IDs should be sorted in descending order
	 *of molecule size. This results in as many larger molecules as possible being exchanged for smaller ones.
	 * @param particleContainer Data structure containing the molecule objects.
	 * @param domainDecomp Active domain decomposition.
	 * @param simstep Actual simulation timestep
	 */
	void controlDensity(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep);
	/**
	 * Calculates the difference between target and actual value of the molecule count for each component and
	 * stores the information in the vector \p vecBalance .
	 * @param[out] vecBalance Vector containing the difference between target and actual value of the molecule count
	 * for each component, where indices refer to the unity based component id. The difference for the total molecule
	 * count is stored at index 0. Too many molecules are indicated by a positive	value, too few by a negative value.
	 * @param pidMap Map containing particle (molecule) IDs of all molecules within specified range \p _range ,
	 * sorted with respect to the respective unity based component IDs.
	 * @param numComponents Number of components present in simulation.
	 */
	void updateBalanceVector(std::vector<int64_t>& vecBalance, std::map<int, std::vector<uint64_t> >& pidMap,
							 uint32_t& numComponents);
	/**
	 * Extracts (tokenizes) numbers from the string \p str , where those numbers are separated by delimiter \p del .
	 * @param[out] vec List (vector) containing the tokenized numbers
	 * @param[in] str String containing numbers separated by delimiter \p del
	 * @param del Delimiter, default: ','
	 * @return Number of tokens stored in the list (vector) \p vec
	 */
	uint32_t tokenize_int_list(std::vector<uint32_t>& vec, std::string str, std::string del = ",");

	void initTargetValues(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp);

private:
	Random _rnd;
	TimestepControl _control;
	DensityTargetType _densityTarget;
	struct Range {
		double xmin, xmax, ymin, ymax, zmin, zmax, volume;
	} _range;
	std::vector<uint32_t> _vecPriority;
#ifdef ENABLE_MPI
	MPI_Datatype pacID_mpi_type;
#endif
};
