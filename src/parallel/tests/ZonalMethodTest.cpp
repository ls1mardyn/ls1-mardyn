/**
 * @file ZonalMethodTest.cpp
 * @author seckler
 * @date 02.10.18
 */

#include "ZonalMethodTest.h"
#include "parallel/ZonalMethods/EighthShell.h"

TEST_SUITE_REGISTRATION(ZonalMethodTest);

void ZonalMethodTest::testES() {
	EighthShell es;
	HaloRegion initialRegion={};
	for (unsigned short d = 0; d < 3; ++d) {
		initialRegion.rmin[d] = 0.;
		initialRegion.rmax[d] = d + 1.;
		initialRegion.offset[d] = 0;
	}
	double cutoffRadius = 0.1;
	initialRegion.width = cutoffRadius;
	bool coversWholeDomain[3] = {false, false, false};
	double cellLength[3] = {0.2, 0.2, 0.2};

	auto haloExportForceImportRegions = es.getHaloExportForceImportRegions(initialRegion, cutoffRadius, coversWholeDomain, cellLength);
	// 7 neighbors
	ASSERT_EQUAL(haloExportForceImportRegions.size(), 7ul);

	for(auto region : haloExportForceImportRegions){
		ASSERT_TRUE(region.offset[0] + region.offset[1] + region.offset[2] < 0 );

		for(unsigned short dim = 0; dim < 3; ++dim){
			ASSERT_TRUE(region.offset[dim] <= 0);
			if (region.offset[dim] == -1) {
				ASSERT_EQUAL(initialRegion.rmin[dim], region.rmin[dim]);
				ASSERT_EQUAL(initialRegion.rmin[dim] + cutoffRadius, region.rmax[dim]);
			} else if (region.offset[dim] == 0) {
				ASSERT_EQUAL(initialRegion.rmin[dim], region.rmin[dim]);
				ASSERT_EQUAL(initialRegion.rmax[dim], region.rmax[dim]);
			} else if (region.offset[dim] == 1) {
				ASSERT_EQUAL(initialRegion.rmax[dim] - cutoffRadius, region.rmin[dim]);
				ASSERT_EQUAL(initialRegion.rmax[dim], region.rmax[dim]);
			}
		}
	}

	auto haloImportForceExportRegions =
		es.getHaloImportForceExportRegions(initialRegion, cutoffRadius, coversWholeDomain, cellLength);
	std::cout << std::endl << "haloImport:" << std::endl;
	for(auto region : haloImportForceExportRegions){
		ASSERT_TRUE(region.offset[0] + region.offset[1] + region.offset[2] > 0 );
		for(unsigned short dim = 0; dim < 3; ++dim){
			ASSERT_TRUE(region.offset[dim] >= 0);
			if (region.offset[dim] == -1) {
				ASSERT_EQUAL(initialRegion.rmin[dim] - cutoffRadius, region.rmin[dim]);
				ASSERT_EQUAL(initialRegion.rmin[dim], region.rmax[dim]);
			} else if (region.offset[dim] == 0) {
				ASSERT_EQUAL(initialRegion.rmin[dim], region.rmin[dim]);
				ASSERT_EQUAL(initialRegion.rmax[dim], region.rmax[dim]);
			} else if (region.offset[dim] == 1) {
				ASSERT_EQUAL(initialRegion.rmax[dim], region.rmin[dim]);
				ASSERT_EQUAL(initialRegion.rmax[dim] + cutoffRadius, region.rmax[dim]);
			}
		}
	}
}
