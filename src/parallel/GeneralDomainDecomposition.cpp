/**
 * @file GeneralDomainDecomposition.cpp
 * @author seckler
 * @date 11.04.19
 */

#include "GeneralDomainDecomposition.h"
#include "ALLLoadBalancer.h"
#include "Domain.h"
#include "NeighbourCommunicationScheme.h"

GeneralDomainDecomposition::GeneralDomainDecomposition(double cutoffRadius, Domain* domain)
	: _boxMin{0.}, _boxMax{0.}, _cutoffRadius{cutoffRadius} {
	std::array<double, 3> domainLength = {domain->getGlobalLength(0), domain->getGlobalLength(1),
										  domain->getGlobalLength(2)};

	auto gridSize = getOptimalGrid(domainLength, this->getNumProcs());
	auto gridCoords = getCoordsFromRank(gridSize, _rank);
	_coversWholeDomain = {gridSize[0] == 1, gridSize[1] == 1, gridSize[2] == 1};
	std::tie(_boxMin, _boxMax) = initializeRegularGrid(domainLength, gridSize, gridCoords);
	_loadBalancer =
		std::make_unique<ALLLoadBalancer>(_boxMin, _boxMax, 4 /*gamma*/, this->getCommunicator(), gridSize, gridCoords);

	global_log->info() << "GeneralDomainDecomposition initial box: [" << _boxMin[0] << ", " << _boxMax[0] << "] x ["
					   << _boxMin[1] << ", " << _boxMax[1] << "] x [" << _boxMin[2] << ", " << _boxMax[2] << "]"
					   << std::endl;
}

double GeneralDomainDecomposition::getBoundingBoxMin(int dimension, Domain* /*domain*/) { return _boxMin[dimension]; }

GeneralDomainDecomposition::~GeneralDomainDecomposition() = default;

double GeneralDomainDecomposition::getBoundingBoxMax(int dimension, Domain* /*domain*/) { return _boxMax[dimension]; }

void GeneralDomainDecomposition::balanceAndExchange(double lastTraversalTime, bool forceRebalancing,
													ParticleContainer* moleculeContainer, Domain* domain) {
	bool doRebalancing = forceRebalancing or ((_steps % _rebuildFrequency == 0 or _steps <= 1));
	if (doRebalancing) {
		if (_steps > 0) {
			DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_ONLY);
		}
		moleculeContainer->deleteOuterParticles();

		std::array<double, 3> newBoxMin{};
		std::array<double, 3> newBoxMax{};
		std::tie(newBoxMin, newBoxMax) = _loadBalancer->rebalance(lastTraversalTime);

        migrateParticles(newBoxMin, newBoxMax);

		std::tie(_boxMin, _boxMax) = std::tie(newBoxMin, newBoxMax);

		_neighbourCommunicationScheme->initCommunicationPartners(_cutoffRadius, domain, this, moleculeContainer);
		DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);
	} else {
		for (int d = 0; d < DIMgeom; ++d) {
			_neighbourCommunicationScheme->setCoverWholeDomain(d, _coversWholeDomain[d]);
		}
		_neighbourCommunicationScheme->initCommunicationPartners(_cutoffRadius, domain, this, moleculeContainer);
		if (sendLeavingWithCopies()) {
			DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_AND_HALO_COPIES);
		} else {
			DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_ONLY);
			moleculeContainer->deleteOuterParticles();
			DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);
		}
	}

	_steps++;
}

void GeneralDomainDecomposition::readXML(XMLfileUnits& xmlconfig) {
	DomainDecompMPIBase::readXML(xmlconfig);
	global_log->info()
		<< "The GeneralDomainDecomposition is enforcing the direct-pp neighbor scheme using fs, so setting it."
		<< std::endl;
	setCommunicationScheme("direct-pp", "fs");
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

std::array<size_t, 3> GeneralDomainDecomposition::getOptimalGrid(const std::array<double, 3>& domainLength,
																 int numProcs) {
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

std::array<size_t, 3> GeneralDomainDecomposition::getCoordsFromRank(const std::array<size_t, 3>& gridSize, int rank) {
	auto yzSize = gridSize[1] * gridSize[2];
	auto zSize = gridSize[2];
	auto x = rank / yzSize;
	auto y = (rank - x * yzSize) / zSize;
	auto z = (rank - x * yzSize - y * zSize);
	return {x, y, z};
}

std::tuple<std::array<double, 3>, std::array<double, 3>> GeneralDomainDecomposition::initializeRegularGrid(
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