//
// Created by alex on 27.05.24.
//

#ifndef MARDYN_FTH_GRID_H
#define MARDYN_FTH_GRID_H

#include "../Interpolation.h"

namespace FTH {
	/**
	 * Data of each node of the FE scheme for density and fth computation
	 * */
	struct NodeData {
		NodeData() : particles(0), density(0), fth({0, 0, 0}), grad({0, 0, 0}) {}
		/// number of particles "around" this node
		int particles;
		/// density normalized by volume
		double density;
		/// approximated fth at this node
		std::array<double, 3> fth;
		/// gradient of density: drho/dx, drho/dy, drho/dz
		std::array<double, 3> grad;
	};
	/// grid data type for density and fth
	using grid_t = Interpolation::FE::Grid<NodeData>;
	using node_t = Interpolation::FE::Node<NodeData>;
	using d3 = Interpolation::FE::d3;
	using i3 = Interpolation::FE::i3;
	using idx_t = Interpolation::FE::idx_t;

	/**
	 * Interpolates the fth based on the provided grid.
	 * */
	d3 interpolateGridFTH(grid_t& grid, const d3& point);

	/**
	 * Writes all fth values of the grid into a file.
	 * */
	[[maybe_unused]] void writeGridFTH(grid_t& grid, const std::string& filename, int simstep=-1);

	/**
	 * Loads all fth values into grid.
	 * */
	[[maybe_unused]] void loadGridFTH(grid_t& grid, const std::string& filename);
}

#endif //MARDYN_FTH_GRID_H