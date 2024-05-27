//
// Created by alex on 05.04.24.
//

#ifndef MARDYN_FTH_H
#define MARDYN_FTH_H

#include "Resolution.h"
#include "../Interpolation.h"
#include "FTH_Grid.h"
#include "../density/Sampler.h"
#include "plugins/AdResS/density/DensityProfile3D.h"

#include <vector>

namespace FTH {
	/**
 	* Struct containing all fields and properties for FTH Handler
 	* */
	struct Config {
		//! @brief True iff F_th functionality should be used
		bool _enableThermodynamicForce;

		//! @brief True iff F_th must be sampled first
		bool _createThermodynamicForce;

		//! @brief enables logging of FTH to file
		bool _logFTH;

		//! @brief enables logging of the current simulation densities to file
		bool _logDensities;

		//! @brief Simulation iterations between each F_th recalculation
		int _thermodynamicForceSampleGap;

		//! @brief Iterations since last sampling
		int _thermodynamicForceSampleCounter;

		//! @brief Describes a accuracy measure in terms of: max|rho(x)-rho_target(x)|/rho_target(x) <= threshold
		//! used to find convergence of F_th
		double _convergenceThreshold;

		//! @brief Controls speed of convergence for F_th
		double _convergenceFactor;

		//! @brief initial and target density
		double _rho0;

		//! @brief Thermodynamic force used to correct the density difference created by plain AdResS
		Interpolation::Function _thermodynamicForce;

		//! @brief Gradient of density distribution, used for convergence checking
		Interpolation::Function _lastGradient;

		//! @brief Pointer to grid, if one is used (can be null)
		grid_t *_grid;

		//! @brief Pointer to the active density sampler (never null)
		SamplerBase *_density_sampler;
	};

	/**
	 * Computes a FTH and applies it within all hybrid regions.
	 * */
	class Handler {
	public:
		Handler() = default;

		/**
		 * Initializes this FTH Handler based on the provided config object.
		 * */
		void init(const Config& config);

		/**
		 * Recomputes F_th in _thermodynamicForce by using the current density profile and interpolating the gradient.
     	 * According to: F_k+1(x) = F_k(x) - c * d'(x)
		 * */
		void computeIteration(ParticleContainer& container, const Resolution::FPRegions_t &regions);

		/**
     	* Checks if the current simulation density is at most _convergenceThreshold apart from _targetDensity.
     	* @returns true iff converged
     	* */
		bool checkConvergence();

		/**
		 * Should be called once every simulation step (before forces are computed). Will update internal counters. Based on configured intervals,
		 * will check for convergence or compute fth iteration.
		 * Afterwards there is no need to manually call checkConvergence or computeIteration!
		 * */
		void step(ParticleContainer& container, const Resolution::FPRegions_t &regions);

		/**
     	* Applies the thermodynamic force to all molecules in hybrid regions.
     	* @param container container of considered molecules
     	* @param regions iterable of all FPRegions
     	* */
		void apply(ParticleContainer& container, const Resolution::FPRegions_t& regions, const Resolution::CompResMap_t& compResMap);

		/**
		 * Initiate writing of density and fth logs depending on the configuration.
		 * If both flags are off then nothing is written
		 * */
		void writeLogs(ParticleContainer &particleContainer, DomainDecompBase &domainDecomp, Domain &domain, unsigned long simstep);

		/**
		 * Write the current version of the FTH regardless of the logging state.
		 * It is assumed that this is called on the finish event of the AdResS plugin.
		 * */
		void writeFinalFTH();

	private:
		//! @brief all FTH properties
		Config _config;

		//! @brief Class to create 3D density profiles
		DensityProfile3D _densityProfiler;
	};
};
#endif //MARDYN_FTH_H
