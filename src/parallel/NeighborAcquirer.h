/**
 * @file NeighborAcquirer.h
 * @author seckler
 * @date 06.05.19
 */

#pragma once
#include <vector>
#include <array>
#include "CommunicationPartner.h"

class Domain;
class HaloRegion;

class NeighborAcquirer {
public:
	/**
	 * Acquire the needed neighbors defined through the specific desired HaloRegions.
	 *
	 * @param domain The domain object.
	 * @param ownRegion The region of the own process.
	 * @param desiredRegions This is a vector of the desired regions. Partners are generated if at least parts of the
	 * desiredRegions lie outside of ownRegion.
	 * @param partners01 Vector of communication partners that contain domains outside of ownRegion.
	 * @param partners02 Vector of communication partners that contain domains inside of ownRegion.
	 * @param comm the mpi communicator
	 * @param excludeOwnRank Specifies to not include CommunicationPartners communicating with the own rank.
	 * @return A tuple of 2 vectors: The first vector represents the partners NOT owning the haloDomain, while the
	 * second vector will own the particles.
	 */
	static std::tuple<std::vector<CommunicationPartner>, std::vector<CommunicationPartner>> acquireNeighbors(
		const std::array<double, 3>& globalDomainLength, HaloRegion* ownRegion, const std::vector<HaloRegion>& desiredRegions,
		const MPI_Comm& comm, bool excludeOwnRank=true);

	static std::vector<CommunicationPartner> squeezePartners(const std::vector<CommunicationPartner>& partners);

private:
	static bool isIncluded(HaloRegion* myRegion, HaloRegion* inQuestion);

	static HaloRegion overlap(const HaloRegion& myRegion, const HaloRegion& inQuestion);

	/**
	 * Calculates all possible shifted regions from the given region.
	 * The shifted regions correspond to the initial region, but wrapped around the periodic boundaries.
	 * All mirror images, as well as, the orginal region are returned.
	 * @param domainLength The total length of the domain.
	 * @param region The initial region.
	 * @return A pair of two vectors of equal length. The first vector contains all shifted regions. The second vector
	 * contains the according shifts.
	 */
	static std::pair<std::vector<HaloRegion>, std::vector<std::array<double, 3>>> getPotentiallyShiftedRegions(
							 const std::array<double, 3>& domainLength, const HaloRegion& region);

    friend class NeighborAcquirerTest;
};
