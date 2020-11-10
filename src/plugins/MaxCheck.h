/*
 * MaxCheck.h
 *
 *  Created on: 28.05.2018
 *      Author: mheinen
 */

#ifndef MAXCHECK_H_
#define MAXCHECK_H_

#include "PluginBase.h"

#include <string>
#include <map>
#include <unordered_map>
#include <list>
#include <cstdint>
#include <vector>

#include "molecules/MoleculeForwardDeclaration.h"

class MaxCheck : public PluginBase
{
private:
	enum MaxCheckMethods
	{
		MCM_UNKNOWN = 0,
		MCM_LIMIT_TO_MAX_VALUE = 1,
		MCM_LIMIT_TO_MAX_VALUE_OVERLAPS = 2,
		MCM_DELETE_PARTICLES = 3
	};

	struct TimestepControl
	{
		uint64_t start;
		uint64_t freq;
		uint64_t stop;
	};

	struct MaxVals
	{
		double F;
		double F2;
		double v;
		double v2;
		uint32_t method;
	};
	typedef std::unordered_map<uint32_t, MaxVals> maxvals_map;

public:
	// constructor and destructor
	MaxCheck();
	~MaxCheck();

	/** @brief Read in XML configuration for MaxCheck and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
		<plugin name="MaxCheck">
			<control>
				<start>INT</start>           <!-- start time step -->
				<frequency>INT</frequency>   <!-- frequency of checking -->
				<stop>INT</stop>             <!-- stop time step -->
			</control>
			<yrange> <min>FLOAT</min> <max>FLOAT</max> </yrange>   <!-- range (y axis) within checking is active -->
			<targets>
				<target cid="INT" method="INT">   <!-- cid: component id of target particles
												  <!-- method: handling particles failing the check, 1:limit to max value | 2:limit to max value if force is too high (overlaps) | 3:delete particle -->
					<Fmax>FLOAT</Fmax>            <! force limit -->
					<vmax>FLOAT</vmax>            <! velocity limit -->
				</target>
			</targets>
		</plugin
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	void init(ParticleContainer *particleContainer,
			  DomainDecompBase *domainDecomp, Domain *domain) override;

    /** @brief Method will be called first thing in a new timestep. */
	void beforeEventNewTimestep(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) override;

    /** @brief Method siteWiseForces will be called before forcefields have been applied
     *  alterations to sitewise forces and fullMolecule forces can be made here
     */

	void siteWiseForces(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) override {}

	/** @brief Method afterForces will be called after forcefields have been applied
	 *
	 * make pure Virtual ?
	 */
	void afterForces(
			ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
			unsigned long simstep
	) override;

	void endStep(
			ParticleContainer *particleContainer,
			DomainDecompBase *domainDecomp, Domain *domain,
			unsigned long simstep
	) override {}

	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain) override {}

	std::string getPluginName() override {return std::string("MaxCheck");}
	static PluginBase* createInstance() {return new MaxCheck();}

private:
	double calcSquaredVectorLength(double* vec) {return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);}
	void checkMaxVals(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep);

private:
	TimestepControl _control;
	maxvals_map _maxVals;
	std::vector<Molecule*> _deletions;
	struct YRange {double min, max;} _yrange;
};

#endif /*MAXCHECK_H_*/
