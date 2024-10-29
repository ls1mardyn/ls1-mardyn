#ifdef MAMICO_COUPLING

#pragma once

#include <coupling/interface/impl/ls1/LS1RegionWrapper.h>
#include <coupling/services/CouplingCellService.h>
#include "PluginBase.h"

/**
 * @brief Allows execution of MaMiCo code to enable coupling with MaMiCo.
 * @author Amartya Das Sharma
 *
 * When MaMiCo is coupled with any simulation software, it needs to control the simulation from both the outside and
 * inside. From the outside, MaMiCo needs to be able to start and stop the simulation, and keep unique simulations
 * distinct. From the inside, MaMiCo code needs to be executed at specific points within the same simulation step. This
 * plugin enables the "inside" behaviour.
 *
 * MaMiCo coupling requires the MAMICO_COUPLING flag to be set, and the MAMICO_SRD_DIR flag to point to MaMiCo source
 * files. With these enabled, ls1 compiles as a library. To make sure the program compiles and works with tests, the
 * relevant MaMiCo portions for the code are put in #ifdef regions.
 *
 * \code{.xml}
 * <plugin name="MamicoCoupling">
 * </plugin>
 * \endcode
 */
class MamicoCoupling : public PluginBase {
public:
	MamicoCoupling() = default;
	~MamicoCoupling() override = default;

	void init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) override;

	/**
	 * No XML tags defined for this plugin, so does nothing
	 */
	void readXML(XMLfileUnits &) {}

	void beforeEventNewTimestep(ParticleContainer *, DomainDecompBase *, unsigned long) {}

	/**
	 * @brief Takes coupling steps such as particle insertion, to make sure they are
	 * accounted for before forces are calculated.
	 *
	 * Following steps are taken, if coupling is switched on:
	 * - Iterate over cells to average values like momentum and mass, to pass to macroscopic solvers
	 * - Distribute incoming mass from macroscopic solver by inserting perticles (if enabled)
	 * - Run the MaMiCo thermostat cell by cell
	 *
	 * The distributeMass method calls the updateParticleContainerAndDecomposition() function at the end, so we end up
	 * with a container with halo particles present. Hence a manual halo clearance is done to restore the state of the
	 * container.
	 */
	void beforeForces(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
					  unsigned long simstep) override;

	/**
	 * @brief Performs adjustments after force calculation
	 *
	 * Following steps are taken, if coupling is switched on:
	 * - Distribute incoming momentum among affected cells
	 * - Apply boundary force to molecules near microscopic domain boundary
	 */
	void afterForces(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
					 unsigned long simstep) override;

	void endStep(ParticleContainer *, DomainDecompBase *, Domain *, unsigned long) {}

	void finish(ParticleContainer *, DomainDecompBase *, Domain *) {}

	std::string getPluginName() override { return "MamicoCoupling"; }

	static PluginBase *createInstance() { return new MamicoCoupling(); }

	/**
	 * @brief Sets the macroscopicCellService object that controls the inner coupling logic.
	 *
	 * MaMiCo extracts the MamicoCoupling plugin object from the simulation object after initialization and uses this
	 * function to set the macroscopicCellCervice.
	 *
	 * The code for this object can be found in
	 * https://github.com/HSU-HPC/MaMiCo/blob/master/coupling/services/CouplingCellService.cpph
	 */
	void setMamicoCouplingCellService(coupling::services::CouplingCellService<3> *couplingCellService) {
		_couplingCellService =
			static_cast<coupling::services::CouplingCellServiceImpl<ls1::LS1RegionWrapper, 3> *>(couplingCellService);
	}

	/**
	 * @brief Enables coupling logic, allowing coupling steps to run while simulation is running.
	 *
	 * Set from within MaMiCo, check
	 * https://github.com/HSU-HPC/MaMiCo/blob/master/coupling/interface/MDSimulationFactory.h, class LS1MDSimulation
	 */
	void switchOnCoupling() { _couplingEnabled = true; }

	/**
	 * @brief Disables coupling logic, not allowing coupling steps to run. Typically done
	 * when equilibrating the simulation initially.
	 *
	 * Set from within MaMiCo, check
	 * https://github.com/HSU-HPC/MaMiCo/blob/master/coupling/interface/MDSimulationFactory.h,
	 * class LS1MDSimulation
	 */
	void switchOffCoupling() { _couplingEnabled = false; }

	/**
	 * @brief Getter method for the coupling state.
	 *
	 * Not used currently, but may be used in future
	 */
	bool getCouplingState() const { return _couplingEnabled; }

private:
	coupling::services::CouplingCellServiceImpl<ls1::LS1RegionWrapper, 3> *_couplingCellService = nullptr;
	bool _couplingEnabled = false;
};

#endif  // MAMICO_COUPLING
