//
// Created by alex on 29.05.24.
//
#include "FTH_Grid.h"
#include "Simulation.h"
#include "parallel/DomainDecompBase.h"

FTH::d3 FTH::interpolateGridFTH(FTH::grid_t &grid, const FTH::d3 &point) {
	const d3& grid_lower = grid.getLower();
	const d3& grid_upper = grid.getUpper();
	d3 p = point;
	for (int dim = 0; dim < 3; dim++) p[dim] = std::clamp(p[dim], grid_lower[dim], grid_upper[dim]);

	auto& element =  grid.getElementOf(p);
	d3 relative_position = element.getRelativeLocalPosition(p);

	d3 fth { 0 };
	static const std::array<d3, 8> coefficients {
			d3 {-1, -1, -1},
			d3 { 1, -1, -1},
			d3 { 1,  1, -1},
			d3 {-1,  1, -1},
			d3 {-1, -1,  1},
			d3 { 1, -1,  1},
			d3 { 1,  1,  1},
			d3 {-1,  1,  1}
	};

	const auto& element_grid_idx = element.getNodes();
	for (int n = 0; n < 8; n++) {
		const auto idx = element_grid_idx[n];
		const d3& node_fth = grid.getNodes().at(idx).data().fth;
		const double factor = 0.125 * (1.0 + coefficients[n][0] * relative_position[0])
							  * (1.0 + coefficients[n][1] * relative_position[1])
							  * (1.0 + coefficients[n][2] * relative_position[2]);
		for (int dim = 0; dim < 3; dim++) {
			fth[dim] += factor * node_fth[dim];
		}
	}

	return fth;
}

void FTH::writeGridFTH(FTH::grid_t &grid, const std::string &filename, int simstep) {
	int rank = 0;
#ifdef ENABLE_MPI
	rank = _simulation.domainDecomposition().getRank();
#endif

	try {
		std::stringstream ss;
		if (simstep < 0) ss << filename << "." << rank << ".txt";
		else ss << filename << "_" << simstep << "." << rank << ".txt";

		std::ofstream file {ss.str()};
		const d3& elem_dims = grid.getElements().getElementSize();
		file << "# Distance between each node: " << elem_dims[0] << "\t" << elem_dims[1] << "\t" << elem_dims[2] << "\n";
		file << "x\ty\tz\tfx\tfy\tfz\n";

		const i3& shape = grid.getNodes().shape();
		for (idx_t z = 0; z < shape[2]; z++) {
			for (idx_t y = 0; y < shape[1]; y++) {
				for (idx_t x = 0; x < shape[0]; x++) {
					file << x << "\t" << y << "\t" << z << "\t";
					const d3& fth = grid.getNodes().at(x, y, z).data().fth;
					file << fth[0] << "\t" << fth[1] << "\t" << fth[2] << "\n";
				}
			}
		}

		file.flush();
		file.close();
	} catch (std::ifstream::failure& e) {
		Log::global_log->error() << "[FTH] Failed to write Interpolation function.\n" << e.what() << std::endl;
		_simulation.exit(-1);
	}
}

void FTH::loadGridFTH(FTH::grid_t &grid, const std::string &filename) {
	int rank = 0;
#ifdef ENABLE_MPI
	rank = _simulation.domainDecomposition().getRank();

	try {
		std::stringstream ss;
		auto ext_pos = filename.find_last_of('.');
		if (ext_pos == std::string::npos) throw std::runtime_error("bad file name - should be filename.txt without rank id");
		ss << filename.substr(0, ext_pos) << "." << rank << ".txt";

		std::ifstream file {ss.str()};
		// skip first 2 lines
		file.ignore(2048, '\n');
		file.ignore(2048, '\n');

		// read data
		std::array<idx_t, 3> coord {};
		while (!file.eof() && file.good()) {
			file >> coord[0] >> coord[1] >> coord[2];
			auto& fth = grid.getNodes().at(coord[0], coord[1], coord[2]).data().fth;
			file >> fth[0] >> fth[1] >> fth[2];
		}

		file.close();
	} catch (std::ifstream::failure& e) {
		Log::global_log->error() << "[FTH] Failed to load Interpolation function.\n" << e.what() << std::endl;
		_simulation.exit(-1);
	}
#endif
}
