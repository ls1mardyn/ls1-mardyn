/**
 * @file MultiSectionMethod.cpp
 * @author seckler
 * @date 11.04.19
 */

#include "MultiSectionMethod.h"
#include "Domain.h"

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

MultiSectionMethod::MultiSectionMethod(double cutoffRadius, Domain* domain)
	: _boxMin{0.}, _boxMax{0.}, _gridSize{0}, _coords{0}, _cutoffRadius{cutoffRadius} {
	std::array<double, 3> domainLength = {domain->getGlobalLength(0), domain->getGlobalLength(1),
										  domain->getGlobalLength(2)};
	_gridSize = getOptimalGrid(domainLength, this->getNumProcs());
}

double MultiSectionMethod::getBoundingBoxMin(int dimension, Domain* /*domain*/) { return _boxMin[dimension]; }

MultiSectionMethod::~MultiSectionMethod() = default;

double MultiSectionMethod::getBoundingBoxMax(int dimension, Domain* /*domain*/) { return _boxMax[dimension]; }

void MultiSectionMethod::balanceAndExchange(double lastTraversalTime, bool forceRebalancing,
											ParticleContainer* moleculeContainer, Domain* domain) {}

void MultiSectionMethod::readXML(XMLfileUnits& xmlconfig) {
	DomainDecompMPIBase::readXML(xmlconfig);
	setCommunicationScheme("direct-pp", "fs");
}
