#include "VelocityExchangeTest.h"

#include "utils/Testing.h"

TEST_SUITE_REGISTRATION(VelocityExchangeTest);

VelocityExchangeTest::VelocityExchangeTest() {}

VelocityExchangeTest::~VelocityExchangeTest() {}

void VelocityExchangeTest::testExchangeVelocities() {

    double delta = 1e-6;  // Tolerate deviation between expected and actual value

    // Region size in x,y,z
    // Domain (cubic) size is 134.266123
    std::array<double, 3> region_warm_min = {  1.2,       2.5,   3.0};
    std::array<double, 3> region_warm_max = {134.266123, 24.2, 111.1};
    std::array<double, 3> region_cold_min = {  7.1,      94.0,   4.0};
    std::array<double, 3> region_cold_max = {131.3,     120.4, 134.266123};

    // Before exchange
    // Values of warmest particle in coldest region for 5 components in x,y,z
    std::array<std::array<double, 3>, 5> velo_warmP;
    std::array<unsigned long, 5> id_warmP;
    // Values of coldest particle in warmest region for 5 components in x,y,z
    std::array<std::array<double, 3>, 5> velo_coldP;
    std::array<unsigned long, 5> id_coldP;

    id_warmP = {408, 457, 531, 136, 90};
    // Taken from input phase space file
    velo_warmP[0] = {-0.0373619231342141450 , 0.186007567735908291, -0.3135313684611428231}; // component 1 (MolID = 408)
    velo_warmP[1] = {0.10432906770468329538, -0.379971100055187416, -0.0262100426827961939}; // component 2 (MolID = 457)
    velo_warmP[2] = {0.25269403194633649479,  0.118841672990568372, -0.3145168493884979987}; // component 3 (MolID = 531)
    velo_warmP[3] = {-0.2215348644272507050, -0.408456081692780692, -0.0992252165419602871}; // component 4 (MolID = 136)
    velo_warmP[4] = {0.02914532040931098327, -0.300353650733340471, -0.0325277836095406794}; // component 5 (MolID = 90)

    id_coldP = {424, 214, 188, 86, 410};
    // Taken from input phase space file
    velo_coldP[0] = {-0.0229580846486864449, -0.079405073547656968,  0.0229694431430437773}; // component 1 (MolID = 424)
    velo_coldP[1] = {0.07933333365111071289,  0.126724820133078786, -0.0899963300407165545}; // component 2 (MolID = 214)
    velo_coldP[2] = {0.01521933222961417177, -0.080162985567112016, -0.0854176829106185892}; // component 3 (MolID = 188)
    velo_coldP[3] = {0.03526942999678589485, -0.034933389801121858, -0.0118687816941973386}; // component 4 (MolID = 86)
    velo_coldP[4] = {-0.0104287815918101215, -0.107184198951575590,  0.0236373245476162070}; // component 5 (MolID = 410)

      // Cutoff might be too small, but has no effect anyway
    std::unique_ptr<ParticleContainer> container{initializeFromFile(ParticleContainerFactory::LinkedCell, "VectorizationMultiComponentMultiPotentials.inp", 6.0)};

    XMLfileUnits inp(getTestDataDirectory()+"/VelocityExchange.xml");

    std::unique_ptr<VelocityExchange> plugin {new VelocityExchange()};

    plugin->init(container.get(), _domainDecomposition, _domain);
    plugin->readXML(inp);
    plugin->beforeForces(container.get(), _domainDecomposition, 0);

    // Check if size of regions is correct (e.g. "box" in xml)
    double min_val;
    double max_val;
    for (unsigned short d = 0; d < 3; d++) {
        // Warm region
        plugin->getRegionCoords(min_val, max_val, d, true);
        ASSERT_DOUBLES_EQUAL_MSG("Box size (min) in warm region in direction "+std::to_string(d)+" not as expected", region_warm_min[d], min_val, delta);
        ASSERT_DOUBLES_EQUAL_MSG("Box size (max) in warm region in direction "+std::to_string(d)+" not as expected", region_warm_max[d], max_val, delta);
        // Cold region
        plugin->getRegionCoords(min_val, max_val, d, false);
        ASSERT_DOUBLES_EQUAL_MSG("Box size (min) in cold region in direction "+std::to_string(d)+" not as expected", region_cold_min[d], min_val, delta);
        ASSERT_DOUBLES_EQUAL_MSG("Box size (max) in cold region in direction "+std::to_string(d)+" not as expected", region_cold_max[d], max_val, delta);
    }

    bool insideBox {false};
    // Search for exchanged particles
    for (auto it = container->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it) {
        const unsigned long pid = it->getID();

#ifndef MARDYN_AUTOPAS
        for (uint32_t cid = 0; cid < 3; cid++) {
            // Formerly warmest particle in cold region
            if (pid == id_warmP[cid]) {
                // Should now be the velocity of the formerly coldest particle
                ASSERT_DOUBLES_EQUAL_MSG("Velocity v_x of formerly warmest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not as expected", velo_coldP[cid][0], it->v(0), delta);
                ASSERT_DOUBLES_EQUAL_MSG("Velocity v_y of formerly warmest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not as expected", velo_coldP[cid][1], it->v(1), delta);
                ASSERT_DOUBLES_EQUAL_MSG("Velocity v_z of formerly warmest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not as expected", velo_coldP[cid][2], it->v(2), delta);
                insideBox = (it->r(0) >= region_warm_min[0]) and (it->r(0) <= region_warm_max[0]);
                ASSERT_EQUAL_MSG("Position r_x of formerly warmest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not in cold region", true, insideBox);
                insideBox = (it->r(1) >= region_warm_min[1]) and (it->r(1) <= region_warm_max[1]);
                ASSERT_EQUAL_MSG("Position r_y of formerly warmest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not in cold region", true, insideBox);
                insideBox = (it->r(2) >= region_warm_min[2]) and (it->r(2) <= region_warm_max[2]);
                ASSERT_EQUAL_MSG("Position r_z of formerly warmest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not in cold region", true, insideBox);
            }

            // Formerly coldest particle in warm region
            if (pid == id_coldP[cid]) {
                // Should now be the velocity of the formerly warmest particle
                ASSERT_DOUBLES_EQUAL_MSG("Velocity v_x of formerly coldest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not as expected", velo_warmP[cid][0], it->v(0), delta);
                ASSERT_DOUBLES_EQUAL_MSG("Velocity v_y of formerly coldest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not as expected", velo_warmP[cid][1], it->v(1), delta);
                ASSERT_DOUBLES_EQUAL_MSG("Velocity v_z of formerly coldest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not as expected", velo_warmP[cid][2], it->v(2), delta);
                insideBox = (it->r(0) >= region_cold_min[0]) and (it->r(0) <= region_cold_max[0]);
                ASSERT_EQUAL_MSG("Position r_x of formerly coldest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not in warm region", true, insideBox);
                insideBox = (it->r(1) >= region_cold_min[1]) and (it->r(1) <= region_cold_max[1]);
                ASSERT_EQUAL_MSG("Position r_y of formerly coldest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not in warm region", true, insideBox);
                insideBox = (it->r(2) >= region_cold_min[2]) and (it->r(2) <= region_cold_max[2]);
                ASSERT_EQUAL_MSG("Position r_z of formerly coldest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not in warm region", true, insideBox);
            }
        }
#else
        // So far, Autopas can only handle one component
        // Warmest and coldest particle is 136 (cid 3) and 86 (cid 3), respectively.
        // Formerly warmest particle in cold region
        if (pid == 136) {
            uint32_t cid = 3;
            // Should now be the velocity of the formerly coldest particle
            ASSERT_DOUBLES_EQUAL_MSG("Velocity v_x of formerly warmest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not as expected", velo_coldP[cid][0], it->v(0), delta);
            ASSERT_DOUBLES_EQUAL_MSG("Velocity v_y of formerly warmest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not as expected", velo_coldP[cid][1], it->v(1), delta);
            ASSERT_DOUBLES_EQUAL_MSG("Velocity v_z of formerly warmest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not as expected", velo_coldP[cid][2], it->v(2), delta);
            insideBox = (it->r(0) >= region_warm_min[0]) and (it->r(0) <= region_warm_max[0]);
            ASSERT_EQUAL_MSG("Position r_x of formerly warmest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not in cold region", true, insideBox);
            insideBox = (it->r(1) >= region_warm_min[1]) and (it->r(1) <= region_warm_max[1]);
            ASSERT_EQUAL_MSG("Position r_y of formerly warmest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not in cold region", true, insideBox);
            insideBox = (it->r(2) >= region_warm_min[2]) and (it->r(2) <= region_warm_max[2]);
            ASSERT_EQUAL_MSG("Position r_z of formerly warmest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not in cold region", true, insideBox);
        }

        // Formerly coldest particle in warm region
        if (pid == 86) {
            uint32_t cid = 3;
            // Should now be the velocity of the formerly warmest particle
            ASSERT_DOUBLES_EQUAL_MSG("Velocity v_x of formerly coldest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not as expected", velo_warmP[cid][0], it->v(0), delta);
            ASSERT_DOUBLES_EQUAL_MSG("Velocity v_y of formerly coldest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not as expected", velo_warmP[cid][1], it->v(1), delta);
            ASSERT_DOUBLES_EQUAL_MSG("Velocity v_z of formerly coldest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not as expected", velo_warmP[cid][2], it->v(2), delta);
            insideBox = (it->r(0) >= region_cold_min[0]) and (it->r(0) <= region_cold_max[0]);
            ASSERT_EQUAL_MSG("Position r_x of formerly coldest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not in warm region", true, insideBox);
            insideBox = (it->r(1) >= region_cold_min[1]) and (it->r(1) <= region_cold_max[1]);
            ASSERT_EQUAL_MSG("Position r_y of formerly coldest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not in warm region", true, insideBox);
            insideBox = (it->r(2) >= region_cold_min[2]) and (it->r(2) <= region_cold_max[2]);
            ASSERT_EQUAL_MSG("Position r_z of formerly coldest particle (ID "+std::to_string(pid)+", cid "+std::to_string(cid+1)+") not in warm region", true, insideBox);
        }
#endif

    }
}
