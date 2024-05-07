//
// Created by alex on 05.04.24.
//

#ifndef MARDYN_RESOLUTION_H
#define MARDYN_RESOLUTION_H

#include "../util/Region.h"
#include "molecules/Component.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"

#include <vector>

namespace Resolution {
	using CompResMap_t = std::vector<ResolutionType>;
	using FPRegions_t = std::vector<FPRegion>;

	/**
 	* Struct containing all fields and properties for Resolution Handler
 	* */
	struct Config {
		/**
		 * Reference to all saved components in the current simulation.
		 * This resource is owned by the active domain instance.
		 * */
		std::vector<Component> *components;

		/**
		 * @brief Active Domain instance of the simulation
		 * */
		Domain *domain;

		/**
		 * @brief Maps each component id to its resolution for faster component switching.
		 * */
		CompResMap_t comp_to_res;

		/**
		 * @brief Container for all areas of interest for AdResS. The ares are currently boxes.
		 * */
		FPRegions_t fpRegions;
	};

	/**
 	* Checks all Resolution of all molecules in the specified regions and sets the correct component accordingly.
 	* Each Molecule must be represented by 3 components FP_M, H_M and CG_M.
 	* These must occur in this sequence in the configuration xml file.
 	* Molecules in FP region get the FP_M component assigned and so forth.
 	* */
	class Handler {
	public:
		Handler() = default;

		/**
		 * Initializes this Resolution Handler based on the provided config object.
		 * */
		void init(Config& config);

		/**
		 * Loops over all molecules in the provided particle container and checks their respective resolution.
		 * Will set the component if change is needed.
		 * @param particleContainer particle container
		 * */
		void checkResolution(ParticleContainer &particleContainer);

		/**
		 * @returns an immutable mapping from component IDs to their resolution
		 * */
		[[nodiscard]] const CompResMap_t& getCompResMap() const;

		/**
		 * @returns an immutable iterable of all FPRegions
		 * */
		[[nodiscard]] const FPRegions_t& getRegions() const;

	private:
		/**
     	* @brief checks the component of @param molecule and sets it to the correct LOD depending the @param targetRes.
     	* */
		void checkMoleculeLOD(Molecule& molecule, ResolutionType targetRes);

		//! @brief all resolution properties
		Config _config;
	};
};


#endif //MARDYN_RESOLUTION_H
