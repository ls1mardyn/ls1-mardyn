/*
 * File:   NeighborAcquirerTest.h
 * Author: bierth, seckler
 *
 * Created on February 27, 2018, 5:01 PM
 */
#include <mpi.h>

#include "NeighborAcquirerTest.h"
#include "parallel/NeighborAcquirer.h"

using namespace std;

TEST_SUITE_REGISTRATION(NeighborAcquirerTest);

NeighborAcquirerTest::NeighborAcquirerTest() {
	_fullShell = new FullShell();
}

NeighborAcquirerTest::~NeighborAcquirerTest() { delete _fullShell; }

void NeighborAcquirerTest::testShiftIfNecessary() {
	HaloRegion region; // rmin, rmax, offset, width
	std::array<double,3> domainLength = {10.0, 10.0, 10.0};

	// region does not need to be shifted
	for(int i = 0; i < 3; i++) region.rmax[i] = 5.0;
	for(int i = 0; i < 3; i++) region.rmin[i] = 3.0;

	{
		auto regionsShiftsPair =
			NeighborAcquirer::getPotentiallyShiftedRegions(domainLength, region, 0.);  // region is within domain box

		ASSERT_EQUAL(regionsShiftsPair.first.size(), 1ul);

		auto haloRegion = regionsShiftsPair.first[0];
		auto shift = regionsShiftsPair.second[0];

		// Sanity to check that shift is 0.
		for (int i = 0; i < 3; i++) ASSERT_EQUAL(shift[i], 0.0);
	}

	for(int i = 0; i < 3; i++) region.rmax[i] = 12.0;
	for(int i = 0; i < 3; i++) region.rmin[i] = 11.0;
	{
		auto regionsShiftsPair =
			NeighborAcquirer::getPotentiallyShiftedRegions(domainLength, region, 0.);  // region is within domain box

		ASSERT_EQUAL(regionsShiftsPair.first.size(), 1ul);

		auto haloRegion = regionsShiftsPair.first[0];
		auto shift = regionsShiftsPair.second[0];

		for (int i = 0; i < 3; i++) {
			ASSERT_EQUAL(shift[i], -10.0);
			shift[i] = 0.0;
		}
	}
	
	for(int i = 0; i < 3; i++) region.rmax[i] = 0.0;
	for(int i = 0; i < 3; i++) region.rmin[i] = -1.0;

	{
		auto regionsShiftsPair =
			NeighborAcquirer::getPotentiallyShiftedRegions(domainLength, region, 0.);  // region is within domain box

		ASSERT_EQUAL(regionsShiftsPair.first.size(), 1ul);

		auto haloRegion = regionsShiftsPair.first[0];
		auto shift = regionsShiftsPair.second[0];

		for(int i = 0; i < 3; i++) { ASSERT_EQUAL(shift[i], 10.0); shift[i] = 0.0; }
	}


	
	
}

void NeighborAcquirerTest::testOverlap() { // assume this one works for now, because you thought about it long and hard.
	HaloRegion region01;
	HaloRegion region02;
	
	for(int i = 0; i < 3; i++) { 
		region01.rmax[i] = 4.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 6.0;
		region02.rmin[i] = 3.0;
	}

	auto overlap = NeighborAcquirer::overlap(region01, region02);
	
	for(int i = 0; i < 3; i++) {
		ASSERT_EQUAL(overlap.rmax[i], 4.0);
		ASSERT_EQUAL(overlap.rmin[i], 3.0);
	}
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 6.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 5.0;
		region02.rmin[i] = 3.0;
	}

	overlap = NeighborAcquirer::overlap(region01, region02);
	
	for(int i = 0; i < 3; i++) {
		ASSERT_EQUAL(overlap.rmax[i], 5.0);
		ASSERT_EQUAL(overlap.rmin[i], 3.0);
	}
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 4.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 3.0;
		region02.rmin[i] = 1.0;
	}
	
	overlap = NeighborAcquirer::overlap(region01, region02);
	
	for(int i = 0; i < 3; i++) {
		ASSERT_EQUAL(overlap.rmax[i], 3.0);
		ASSERT_EQUAL(overlap.rmin[i], 2.0);
	}
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 4.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 6.0;
		region02.rmin[i] = 1.0;
	}
	
	overlap = NeighborAcquirer::overlap(region01, region02);
	
	for(int i = 0; i < 3; i++) {
		ASSERT_EQUAL(overlap.rmax[i], 4.0);
		ASSERT_EQUAL(overlap.rmin[i], 2.0);
	}
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 6.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 6.0;
		region02.rmin[i] = 2.0;
	}
	
	overlap = NeighborAcquirer::overlap(region01, region02);
	
	for(int i = 0; i < 3; i++) {
		ASSERT_EQUAL(overlap.rmax[i], 6.0);
		ASSERT_EQUAL(overlap.rmin[i], 2.0);
	}
}

void NeighborAcquirerTest::testIOwnThis() { // i own a part of this
	HaloRegion region01;
	HaloRegion region02;
	
	for(int i = 0; i < 3; i++) { 
		region01.rmax[i] = 4.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 6.0;
		region02.rmin[i] = 3.0;
	}
	
	ASSERT_EQUAL(NeighborAcquirer::isIncluded(&region01, &region02), true);
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 6.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 5.0;
		region02.rmin[i] = 3.0;
	}
	
	ASSERT_EQUAL(NeighborAcquirer::isIncluded(&region01, &region02), true);
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 4.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 3.0;
		region02.rmin[i] = 1.0;
	}
	
	ASSERT_EQUAL(NeighborAcquirer::isIncluded(&region01, &region02), true);
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 4.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 6.0;
		region02.rmin[i] = 1.0;
	}
	
	ASSERT_EQUAL(NeighborAcquirer::isIncluded(&region01, &region02), true);
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 6.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 6.0;
		region02.rmin[i] = 2.0;
	}
	
	ASSERT_EQUAL(NeighborAcquirer::isIncluded(&region01, &region02), true);
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 5.0;
		region01.rmin[i] = 1.0;
	}
	
	region02.rmax[0] = 6.0;
	region02.rmax[1] = 7.0;
	region02.rmax[2] = 5.0;
	
	region02.rmin[0] = 5.0;
	region02.rmin[1] = 6.0;
	region02.rmin[2] = 5.0;
	
	ASSERT_EQUAL(NeighborAcquirer::isIncluded(&region01, &region02), false);
	
	
}

void NeighborAcquirerTest::testCorrectNeighborAcquisition() {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int numRanks;
	MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

	if (numRanks != 2) {
		std::cout << "SKIPPED: requires two processes, but run with "<< numRanks <<"." << std::endl;
		// test is only meant for two processes!
		return;
	}
	double cutoff = 2.5;

	std::array<HaloRegion, 2> ownRegionArray{
		HaloRegion{32.8465, 140.73, 66.3138, 66.2862, 296.984, 99.5571, 0, 0, 0, cutoff},
		HaloRegion{65.9145, 0, 99.5571, 99.4885, 142.386, 132.6, 0, 0, 0, cutoff}};
	auto ownRegion = ownRegionArray[rank];
	auto otherRegion = ownRegionArray[1 - rank];

	std::array<double, 3> globalDomainLength{132.6, 591.891, 132.6};

	std::array<bool, 3> coversWholeDomain{false, false, false};

	std::vector<HaloRegion> leavingRegions =
		_fullShell->getLeavingExportRegions(ownRegion, cutoff, coversWholeDomain.data());

	std::vector<CommunicationPartner> leavingExportNeighbours;
	std::vector<CommunicationPartner> leavingImportNeighbours;
	std::tie(leavingExportNeighbours, leavingImportNeighbours) =
		NeighborAcquirer::acquireNeighbors(globalDomainLength, &ownRegion, leavingRegions, 0., MPI_COMM_WORLD);
	// p1 notes reply, p2 notes owned as leaving import

	for (auto& neighbor : leavingExportNeighbours) {
		for (auto& haloRegion : neighbor._haloInfo) {
			for (auto i = 0; i < 3; ++i) {
				std::stringstream ss;
				ss << "in leavingExport: " << std::endl;
				neighbor.print(ss);
				ASSERT_TRUE_MSG(ss.str(),haloRegion._leavingLow[i] >= otherRegion.rmin[i]);
				ASSERT_TRUE_MSG(ss.str(),haloRegion._leavingHigh[i] <= otherRegion.rmax[i]);
			}
		}
	}
	for (auto& neighbor : leavingImportNeighbours) {
		for (auto& haloRegion : neighbor._haloInfo) {
			for (auto i = 0; i < 3; ++i) {
				std::stringstream ss;
				ss << "in leavingImport: " << std::endl;
				neighbor.print(ss);
				ASSERT_TRUE_MSG(ss.str(),haloRegion._leavingLow[i] >= ownRegion.rmin[i]);
				ASSERT_TRUE_MSG(ss.str(),haloRegion._leavingHigh[i] <= ownRegion.rmax[i]);
			}
		}
	}
}