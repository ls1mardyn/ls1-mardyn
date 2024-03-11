//
// Created by Alex Hocks on 05.03.24.
//

#ifndef MARDYN_LANGEVIN_H
#define MARDYN_LANGEVIN_H


#include "Integrator.h"

#include <array>
#include <vector>

/**
 * New time integrator that deviates from the standard Verlet integration scheme.
 * Equations of motion are changed to follow Langevin equations.
 * This is a stochastic integration scheme, which is incorporated as random "kicks" from the connected heat bath.
 *
 * This is to be used with the TemperatureObserver, which disables the normal velocity scaling and defines
 * the regions in which the Langevin Thermostat should be active.
 * */
class Langevin : public Integrator {
public:
	/**
	 * XML Format:
	 * <integrator type="Langevin">
	 *     <timestep>DOUBLE</timestep>
	 *     <friction>DOUBLE</friction>
	 *     <ActiveRegions>
	 *         <region>
	 *             <lower><x>DOUBLE</x> <y>DOUBLE</y> <z>DOUBLE</z></lower>
	 *             <upper><x>DOUBLE</x> <y>DOUBLE</y> <z>DOUBLE</z></upper>
	 *         </region>
	 *     </ActiveRegions>
	 * </integrator>
	 * */
	void readXML(XMLfileUnits &xmlconfig) override;

	void init() override;

	void eventNewTimestep(ParticleContainer *moleculeContainer, Domain *domain) override;

	void eventForcesCalculated(ParticleContainer *moleculeContainer, Domain *domain) override;

private:
	using d3 = std::array<double, 3>;

	struct box_t {
		box_t() = default;

		box_t(const d3 &l, const d3 &h) : low(l), high(h) {}

		d3 low, high;
	};

	//! @brief friction strength
	double _gamma;

	//! @brief time step length halved
	double _dt_half;

	//! @brief regions in which friction is not set to 0
	std::vector<box_t> _stochastic_regions;

	//! @brief We need to check at least twice if a TemperatureObserver exists
	bool _checkFailed;

	/**
	 * Sample from Gaussian with technically 0 mean and sigma**2 = 2 * m * friction * k_b * T_target / delta_t.
	 * Is already adapted to match the integration scheme.
	 * @param m mass
	 * @param T temp target
	 * */
	d3 sampleRandomForce(double m, double T);

	/**
	 * Adds the parts of the equation, that only exist in the Langevin Equation's of motion
	 * */
	 void addLangevinContribution(ParticleContainer* particleContainer);
};


#endif //MARDYN_LANGEVIN_H
