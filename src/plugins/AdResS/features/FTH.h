//
// Created by alex on 05.04.24.
//

#ifndef MARDYN_FTH_H
#define MARDYN_FTH_H

#include "Resolution.h"
#include "../Interpolation.h"
#include "FTH_Grid.h"
#include "../density/Sampler.h"

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

		virtual ~Handler() = default;

		/**
		 * Initializes this FTH Handler based on the provided config object.
		 * */
		virtual void init(const Config& config) = 0;

		/**
		 * Recomputes F_th in _thermodynamicForce by using the current density profile and interpolating the gradient.
     	 * According to: F_k+1(x) = F_k(x) - c * d'(x)
		 * */
		virtual void updateForce(ParticleContainer& container, const Resolution::FPRegions_t &regions) = 0;

		/**
     	* Checks if the current simulation density is at most _convergenceThreshold apart from _targetDensity.
     	* @returns true iff converged
     	* */
		virtual bool checkConvergence() = 0;

		/**
		 * Should be called once every simulation step (before forces are computed). Will update internal counters. Based on configured intervals,
		 * will check for convergence or compute fth iteration.
		 * Afterwards there is no need to manually call checkConvergence or updateForce!
		 * */
		void computeSingleIteration(ParticleContainer& container, const Resolution::FPRegions_t &regions);

		/**
     	* Applies the thermodynamic force to all molecules in hybrid regions.
     	* @param container container of considered molecules
     	* @param regions iterable of all FPRegions
     	* */
		virtual void applyForce(ParticleContainer& container, const Resolution::FPRegions_t& regions, const Resolution::CompResMap_t& compResMap) = 0;

		/**
		 * Initiate writing of density and fth logs depending on the configuration.
		 * If both flags are off then nothing is written
		 * */
		virtual void writeLogs(ParticleContainer &particleContainer, DomainDecompBase &domainDecomp, Domain &domain, unsigned long simstep) = 0;

		/**
		 * Write the current version of the FTH regardless of the logging state.
		 * It is assumed that this is called on the finish event of the AdResS plugin.
		 * */
		virtual void writeFinalFTH() = 0;

	protected:
		//! @brief all FTH properties
		Config _config;
	};

	/**
	 * Computes FTH for a density sampler, that uses a 3D grid
	 * */
	class Grid3DHandler : public Handler {
	public:
		void init(const Config &config) override;

		void updateForce(ParticleContainer &container, const Resolution::FPRegions_t &regions) override;

		bool checkConvergence() override;

		void applyForce(ParticleContainer &container, const Resolution::FPRegions_t &regions, const Resolution::CompResMap_t &compResMap) override;

		void writeLogs(ParticleContainer &particleContainer, DomainDecompBase &domainDecomp, Domain &domain, unsigned long simstep) override;

		void writeFinalFTH() override;
	};

	/**
	 * Computes FTH for a density sampler, that uses a 1D grid
	 * */
	class Grid1DHandler : public Handler {
	public:
		void init(const Config &config) override;

		void updateForce(ParticleContainer &container, const Resolution::FPRegions_t &regions) override;

		bool checkConvergence() override;

		void applyForce(ParticleContainer &container, const Resolution::FPRegions_t &regions, const Resolution::CompResMap_t &compResMap) override;

		void writeLogs(ParticleContainer &particleContainer, DomainDecompBase &domainDecomp, Domain &domain, unsigned long simstep) override;

		void writeFinalFTH() override;

	private:
		//! @brief Thermodynamic force used to correct the density difference created by plain AdResS
		Interpolation::Function _thermodynamicForce;
	};

	/**
	 * Computes FTH for a density sampler, that uses a 3D function
	 * */
	class Function3DHandler : public Handler {
	public:
		void init(const Config &config) override;

		void updateForce(ParticleContainer &container, const Resolution::FPRegions_t &regions) override;

		bool checkConvergence() override;

		void applyForce(ParticleContainer &container, const Resolution::FPRegions_t &regions, const Resolution::CompResMap_t &compResMap) override;

		void writeLogs(ParticleContainer &particleContainer, DomainDecompBase &domainDecomp, Domain &domain, unsigned long simstep) override;

		void writeFinalFTH() override;

	private:
		//! @brief Thermodynamic force used to correct the density difference created by plain AdResS
		Interpolation::Function3D _thermodynamicForce;
	};

	/**
	 * Computes FTH for a density sampler, that uses a 1D function
	 * */
	class Function1DHandler : public Handler {
	public:
		void init(const Config &config) override;

		void updateForce(ParticleContainer &container, const Resolution::FPRegions_t &regions) override;

		bool checkConvergence() override;

		void applyForce(ParticleContainer &container, const Resolution::FPRegions_t &regions, const Resolution::CompResMap_t &compResMap) override;

		void writeLogs(ParticleContainer &particleContainer, DomainDecompBase &domainDecomp, Domain &domain, unsigned long simstep) override;

		void writeFinalFTH() override;

	private:
		//! @brief Thermodynamic force used to correct the density difference created by plain AdResS
		Interpolation::Function _thermodynamicForce;
	};
};
#endif //MARDYN_FTH_H
