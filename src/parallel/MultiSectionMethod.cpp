/**
 * @file MultiSectionMethod.cpp
 * @author seckler
 * @date 11.04.19
 */

#include "MultiSectionMethod.h"
#include "Domain.h"
#include "NeighbourCommunicationScheme.h"

MultiSectionMethod::MultiSectionMethod(double cutoffRadius, Domain* domain)
	: _boxMin{0.},
	  _boxMax{0.},
	  _gridSize{0},
	  _gridCoords{0},
	  _cutoffRadius{cutoffRadius},
	  _step{0},
	  _updateFrequency{10} {
	std::array<double, 3> domainLength = {domain->getGlobalLength(0), domain->getGlobalLength(1),
										  domain->getGlobalLength(2)};

	_gridSize = getOptimalGrid(domainLength, this->getNumProcs());
	_gridCoords = getCoordsFromRank(_gridSize, _rank);
	std::tie(_boxMin, _boxMax) = initializeRegularGrid(domainLength, _gridSize, _gridCoords);
	global_log->info() << "MultiSectionMethod initial box: [" << _boxMin[0] << ", " << _boxMax[0] << "] x ["
					   << _boxMin[1] << ", " << _boxMax[1] << "] x [" << _boxMin[2] << ", " << _boxMax[2] << "]"
					   << std::endl;
}

double MultiSectionMethod::getBoundingBoxMin(int dimension, Domain* /*domain*/) { return _boxMin[dimension]; }

MultiSectionMethod::~MultiSectionMethod() = default;

double MultiSectionMethod::getBoundingBoxMax(int dimension, Domain* /*domain*/) { return _boxMax[dimension]; }

bool MultiSectionMethod::queryRebalancing(size_t step, size_t updateFrequency, double /*lastTraversalTime*/) {
	return step % updateFrequency == 0;
}

void MultiSectionMethod::balanceAndExchange(double lastTraversalTime, bool forceRebalancing,
											ParticleContainer* moleculeContainer, Domain* domain) {
	bool rebalance = queryRebalancing(_step, _updateFrequency, lastTraversalTime) or forceRebalancing;
	if (_step == 0) {
		// ensure that there are no outer particles
		moleculeContainer->deleteOuterParticles();
		// init communication partners
		initCommPartners(moleculeContainer, domain);
	} else {
		if (rebalance) {
			// first transfer leaving particles
			DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_ONLY);

			// ensure that there are no outer particles
			moleculeContainer->deleteOuterParticles();

			// rebalance
			std::array<double, 3> newBoxMin{0.}, newBoxMax{0.};
			std::tie(newBoxMin, newBoxMax) = doRebalancing(lastTraversalTime, moleculeContainer, domain);

			// migrate the particles, this will rebuild the moleculeContainer!
			migrateParticles(domain, moleculeContainer, newBoxMin, newBoxMax);

			// set new boxMin and boxMax
			_boxMin = newBoxMin;
			_boxMax = newBoxMax;

			// init communication partners
			initCommPartners(moleculeContainer, domain);
		} else {
			if (sendLeavingWithCopies()) {
				DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_AND_HALO_COPIES);
			} else {
				DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_ONLY);
				moleculeContainer->deleteOuterParticles();
				DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);
			}
		}
	}
	++_step;
}

std::tuple<std::array<double, 3>, std::array<double, 3>> MultiSectionMethod::doRebalancing(
	double lastTraversalTime, ParticleContainer* particleContainer, Domain* domain) {
	return std::make_tuple(_boxMin, _boxMax);
}

void MultiSectionMethod::migrateParticles(Domain* domain, ParticleContainer* particleContainer, array<double, 3> newMin,
										  array<double, 3> newMax) {
	std::array<double, 3> oldBoxMin{particleContainer->getBoundingBoxMin(0), particleContainer->getBoundingBoxMin(1),
								 particleContainer->getBoundingBoxMin(2)};
	std::array<double, 3> oldBoxMax{particleContainer->getBoundingBoxMax(0), particleContainer->getBoundingBoxMax(1),
								 particleContainer->getBoundingBoxMax(2)};


}

void MultiSectionMethod::initCommPartners(ParticleContainer* moleculeContainer,
										  Domain* domain) {  // init communication partners
	for (int d = 0; d < DIMgeom; ++d) {
		// this needs to be updated for proper initialization of the neighbours
		_neighbourCommunicationScheme->setCoverWholeDomain(d, _gridSize[d] == 1);
	}
	_neighbourCommunicationScheme->initCommunicationPartners(_cutoffRadius, domain, this, moleculeContainer);
	exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);
}

void MultiSectionMethod::readXML(XMLfileUnits& xmlconfig) {
	DomainDecompMPIBase::readXML(xmlconfig);
	global_log->info() << "The MultiSectionMethod is enforcing the direct-pp neighbor scheme using fs, so setting it."
					   << std::endl;
	setCommunicationScheme("direct-pp", "fs");

	xmlconfig.getNodeValue("updateFrequency", _updateFrequency);
	global_log->info() << "MultiSectionMethod update frequency: " << _updateFrequency << endl;
}

/**
 * Get the ordering of the input data.
 * The ordering will contain indices of elements of data, starting with the smallest going to the biggest.
 * e.g.:
 * returns {2, 0, 3, 1} for data={5, 16, 4, 7}
 * @tparam ArrayType should be some array type, e.g., vector of int
 * @param data the data to define the order
 * @return The ordering.
 */
template <typename ArrayType>
std::vector<size_t> getOrdering(const ArrayType& data) {
	std::vector<size_t> index(data.size(), 0);
	for (int i = 0; i != index.size(); i++) {
		index[i] = i;
	}
	std::sort(index.begin(), index.end(), [&](const size_t& a, const size_t& b) { return (data[a] < data[b]); });
	return index;
}

std::array<size_t, 3> MultiSectionMethod::getOptimalGrid(const std::array<double, 3>& domainLength, int numProcs) {
	// generate default grid
	std::array<int, 3> gridSize{0};
	MPI_CHECK(MPI_Dims_create(numProcs, 3, gridSize.data()));
	std::sort(gridSize.begin(), gridSize.end());

	// we want the grid to actually resemble the domain a bit (long sides -> many processes!)
	std::array<size_t, 3> grid{0};
	auto ordering = getOrdering(domainLength);
	for (size_t i = 0; i < 3; ++i) {
		grid[ordering[i]] = gridSize[i];
	}
	return grid;
}

std::array<size_t, 3> MultiSectionMethod::getCoordsFromRank(const std::array<size_t, 3>& gridSize, int rank) {
	auto yzSize = gridSize[1] * gridSize[2];
	auto zSize = gridSize[2];
	auto x = rank / yzSize;
	auto y = (rank - x * yzSize) / zSize;
	auto z = (rank - x * yzSize - y * zSize);
	return {x, y, z};
}

std::tuple<std::array<double, 3>, std::array<double, 3>> MultiSectionMethod::initializeRegularGrid(
	const std::array<double, 3>& domainLength, const std::array<size_t, 3>& gridSize,
	const std::array<size_t, 3>& gridCoords) {
	std::array<double, 3> boxMin{0.};
	std::array<double, 3> boxMax{0.};
	// initialize it as regular grid!
	for (size_t dim = 0; dim < 3; ++dim) {
		boxMin[dim] = gridCoords[dim] * domainLength[dim] / gridSize[dim];
		boxMax[dim] = (gridCoords[dim] + 1) * domainLength[dim] / gridSize[dim];
		if (gridCoords[dim] == gridSize[dim] - 1) {
			// ensure that the upper domain boundaries match.
			// lower domain boundaries always match, because they are 0.
			boxMax[dim] = domainLength[dim];
		}
	}
	return std::make_tuple(boxMin, boxMax);
}
