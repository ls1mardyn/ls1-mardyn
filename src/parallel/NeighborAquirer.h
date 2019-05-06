/**
 * @file NeighborAquirer.h
 * @author seckler
 * @date 06.05.19
 */

#pragma once
#include <vector>
#include "CommunicationPartner.h"

class Domain;
class HaloRegion;

class NeighborAquirer {
public:
	static void aquireNeighbours(Domain* domain, HaloRegion* haloRegion, std::vector<HaloRegion>& desiredRegions,
								 std::vector<CommunicationPartner>& partners01,
								 std::vector<CommunicationPartner>& partners02);
	static std::vector<CommunicationPartner> squeezePartners(const std::vector<CommunicationPartner>& partners);

private:
	static bool isIncluded(HaloRegion* myRegion, HaloRegion* inQuestion);
	static void overlap(HaloRegion* myRegion, HaloRegion* inQuestion);
	static void shiftIfNeccessary(const double* domainLength, HaloRegion* region, double* shiftArray);

	friend class NeighbourCommunicationSchemeTest;
};
