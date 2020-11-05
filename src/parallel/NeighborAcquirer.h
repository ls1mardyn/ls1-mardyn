/**
 * @file NeighborAcquirer.h
 * @author seckler
 * @date 06.05.19
 */

#pragma once
#include <vector>
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
	 * @return A tuple of 2 vectors: The first vector represents the partners NOT owning the haloDomain, while the
	 * second vector will own the particles.
	 */
	static std::tuple<std::vector<CommunicationPartner>, std::vector<CommunicationPartner>> acquireNeighbors(
		const std::array<double,3>& globalDomainLength, HaloRegion* ownRegion, std::vector<HaloRegion>& desiredRegions, double skin, const MPI_Comm& comm);

	static std::vector<CommunicationPartner> squeezePartners(const std::vector<CommunicationPartner>& partners);

private:
	static bool isIncluded(HaloRegion* myRegion, HaloRegion* inQuestion);

	static HaloRegion overlap(const HaloRegion& myRegion, const HaloRegion& inQuestion);

	static HaloRegion getPotentiallyShiftedRegion(const std::array<double,3>& domainLength, const HaloRegion& region,
												  double* shiftArray, double skin);

	/**
	 * Get all possible combinations of halo regions and shifts, where the given halo region nonShiftedRegion can get
	 * particles from.
	 * For that all possible combinations of the shifted and non-shifted regions will be taken.
	 * If, e.g., the shift happens in the dimensions 1(+) and 2(-) and no shift happens in the direction 0, the returned halo regions will have the shifts:
	 * (0,0,0), (0,0,-), (0,+,0), (0,+,-).
	 * If the shift has 0 non-zero entries, the length of the returned vectors is 1=2^0.
	 * If the shift has 1 non-zero entries, the length of the returned vectors is 2=2^1.
	 * If the shift has 2 non-zero entries, the length of the returned vectors is 4=2^2.
	 * If the shift has 3 non-zero entries, the length of the returned vectors is 8=2^3.
	 * @param nonShiftedRegion
	 * @param shiftedRegion
	 * @param shift
	 * @return tuple of vectors of halo regions with the according shifts.
	 */
	static std::tuple<std::vector<HaloRegion>, std::vector<std::array<double, 3>>>
	getAllShiftedAndNonShiftedRegionsAndShifts(HaloRegion nonShiftedRegion, HaloRegion shiftedRegion,
											  std::array<double, 3> shift);

    friend class NeighborAcquirerTest;
};
