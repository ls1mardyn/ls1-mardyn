//
// Created by alex on 31.05.23.
//

#include "AdResSTest.h"
#include "plugins/AdResS/AdResS.h"

TEST_SUITE_REGISTRATION(AdResSTest);

AdResSTest::AdResSTest() = default;

AdResSTest::~AdResSTest() = default;

void AdResSTest::computeForcesTest() {
    std::unique_ptr<AdResS> plugin;
    ParticleContainer* container = nullptr;
    FPRegion region({4,4,4}, {6,6,6}, {2,2,2});
    //init
    {
        const char* filename = "AdResS-empty-10x10x10.inp";
        container = initializeFromFile(ParticleContainerFactory::LinkedCell, filename, 2);

        plugin = std::make_unique<AdResS>();
        plugin->_particleContainer = container;

        region.init();
        plugin->_fpRegions.emplace_back(region);
        plugin->_components = _simulation.getEnsemble()->getComponents();
        for(int i = 0; i < 6; i++) plugin->_comp_to_res[i] = static_cast<Resolution>(i % 3);
        plugin->_domain = _simulation.getDomain();
        plugin->weight = plugin->weightNearest;
    }
    ParticlePairs2PotForceAdapter pairsHandler(*_simulation.getDomain());

    /*
     * Since Cutoff is set to 2 and the domain is 10x10x10 with a FPRegion from 4,4,4 to 6,6,6 with hybrid dims 2,2,2
     * the entire domain is handled by the plugin.
     * */

    //check every single cell -> place 2 molecules in one cell
    for(int cz = 1; cz <= 9; cz+=2) {
        for (int cy = 1; cy <= 9; cy+=2) {
            for (int cx = 1; cx <= 9; cx+=2) {
                Molecule m1(0, &plugin->_components->at(2), cx-0.5,cy,cz, 0,0,0, 1,0,0,0);
                Molecule m2(1, &plugin->_components->at(2), cx+0.5,cy,cz, 0,0,0, 1,0,0,0);
                container->addParticle(m1);
                container->addParticle(m2);
                plugin->beforeForces(container, nullptr, 0);
                Molecule* mol1;
                Molecule* mol2;
                for(auto it = container->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it){
                    if(it->getID() == 0) mol1 = it.operator->();
                    if(it->getID() == 1) mol2 = it.operator->();
                }

                mol1->buildOwnSoA();
                mol2->buildOwnSoA();
                std::array<double,3> distanceVector = {};
                double dd = m2.dist2(m1, distanceVector.data());
                if(mol1->componentid() == 0 && 0 == mol2->componentid()) goto clear0;
                if(mol1->componentid() == 2 && 2 == mol2->componentid()) goto clear0;

                pairsHandler.init();
                pairsHandler.processPair(*mol1, *mol2, distanceVector.data(), MOLECULE_MOLECULE, dd, true);
                pairsHandler.finish();
                plugin->siteWiseForces(container, nullptr, 0);
                if(mol1->componentid() == 0 && mol2->componentid() == 1) {
                    double weight = plugin->weight(mol2->r_arr(), region);
                    ASSERT_DOUBLES_EQUAL(-24. * weight, mol1->ljcenter_F(0)[0], 0);
                    ASSERT_DOUBLES_EQUAL(24. * weight, mol2->ljcenter_F(0)[1], 0);
                } else if( mol1->componentid() == 1 && mol2->componentid() == 0) {
                    double weight = plugin->weight(mol1->r_arr(), region);
                    ASSERT_DOUBLES_EQUAL(-24. * weight, mol1->ljcenter_F(0)[1], 0);
                    ASSERT_DOUBLES_EQUAL(24. * weight, mol2->ljcenter_F(0)[0], 0);
                } else {
                    ASSERT_DOUBLES_EQUAL(-24., mol1->ljcenter_F(0)[0] + ((mol1->componentid() == 1) ? mol1->ljcenter_F(1)[0] : 0), 0);
                    ASSERT_DOUBLES_EQUAL(24., mol2->ljcenter_F(0)[0] + ((mol2->componentid() == 1) ? mol2->ljcenter_F(1)[0] : 0), 0);
                }

                clear0:
                mol1->releaseOwnSoA();
                mol2->releaseOwnSoA();
                container->clear();
            }
        }
    }

    //check neighbouring cells
    for(int col = 1; col < 8; col++) {
        for(int cz = 1; cz <= 9; cz+=2) {
            for (int cy = 1; cy <= 9; cy+=2) {
                for (int cx = 1; cx <= 9; cx+=2) {
                    std::array<int,3> offset = threeDimensionalMapping::oneToThreeD(col, {2,2,2});
                    if(cx + offset[0] >= 10 || cy + offset[1] >= 10 || cz + offset[2] >= 10) continue;

                    const double offset_l2 = sqrt(std::pow(offset[0]*2,2) + std::pow(offset[1]*2,2) + std::pow(offset[2]*2,2));
                    const double offset_scale = 2 / offset_l2;

                    Molecule m1(0, &plugin->_components->at(2), cx,cy,cz, 0,0,0, 1,0,0,0);
                    Molecule m2(1, &plugin->_components->at(2), cx+offset[0]*2*offset_scale*0.9,cy+offset[1]*2*offset_scale*0.9,cz+offset[2]*2*offset_scale*0.9, 0,0,0, 1,0,0,0);
                    container->addParticle(m1);
                    container->addParticle(m2);
                    plugin->beforeForces(container, nullptr, 0);
                    Molecule* mol1;
                    Molecule* mol2;
                    for(auto it = container->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it){
                        if(it->getID() == 0) mol1 = it.operator->();
                        if(it->getID() == 1) mol2 = it.operator->();
                    }

                    mol1->buildOwnSoA();
                    mol2->buildOwnSoA();
                    std::array<double,3> distanceVector = {};
                    std::array<double,3> baseForce = {};
                    double dd = m2.dist2(m1, distanceVector.data());
                    if(mol1->componentid() == 0 && 0 == mol2->componentid()) goto clear1;
                    if(mol1->componentid() == 2 && 2 == mol2->componentid()) goto clear1;

                    pairsHandler.init();
                    pairsHandler.processPair(*mol1, *mol2, distanceVector.data(), MOLECULE_MOLECULE, dd, true);
                    pairsHandler.finish();
                    for(int d = 0; d < 3; d++) baseForce[d] += mol1->ljcenter_F(0)[d] + ((mol1->componentid() == 1) ? mol1->ljcenter_F(1)[d] : 0);

                    plugin->siteWiseForces(container, nullptr, 0);
                    if(mol1->componentid() == 1 && mol2->componentid() == 1) {
                        auto f1 = mol1->ljcenter_F(0);
                        auto f2 = mol2->ljcenter_F(0);
                        for(int d = 0; d<3; d++) {
                            f1[d] += mol1->ljcenter_F(1)[d];
                            f2[d] += mol2->ljcenter_F(1)[d];
                            ASSERT_TRUE(f2[d] == -f1[d]);
                        }
                        ASSERT_TRUE(f1[0] != 0 || f1[1] != 0 || f1[2] != 0);
                        ASSERT_TRUE(f1[0] != baseForce[0] || f1[1] != baseForce[1] || f1[2] != baseForce[2]);
                    } else if(mol1->componentid() == 1 && mol2->componentid() != 1) {
                        auto f1 = mol1->ljcenter_F(0);
                        auto f2 = mol2->ljcenter_F(0);
                        for(int d = 0; d<3; d++) {
                            f1[d] += mol1->ljcenter_F(1)[d];
                            ASSERT_TRUE(f2[d] == -f1[d]);
                        }
                        ASSERT_TRUE(f1[0] != 0 || f1[1] != 0 || f1[2] != 0);
                        ASSERT_TRUE(f1[0] != baseForce[0] || f1[1] != baseForce[1] || f1[2] != baseForce[2]);
                    } else if(mol1->componentid() != 1 && mol2->componentid() == 1) {
                        auto f1 = mol1->ljcenter_F(0);
                        auto f2 = mol2->ljcenter_F(0);
                        for(int d = 0; d<3; d++) {
                            f2[d] += mol2->ljcenter_F(1)[d];
                            ASSERT_TRUE(f2[d] == -f1[d]);
                        }
                        ASSERT_TRUE(f1[0] != 0 || f1[1] != 0 || f1[2] != 0);
                        ASSERT_TRUE(f1[0] != baseForce[0] || f1[1] != baseForce[1] || f1[2] != baseForce[2]);
                    } else {
                        auto f1 = mol1->ljcenter_F(0);
                        auto f2 = mol2->ljcenter_F(0);
                        for(int d = 0; d<3; d++) {
                            ASSERT_TRUE(f2[d] == -f1[d]);
                        }
                        ASSERT_TRUE(f1[0] != 0 || f1[1] != 0 || f1[2] != 0);
                        ASSERT_TRUE(f1[0] != baseForce[0] || f1[1] != baseForce[1] || f1[2] != baseForce[2]);
                    }

                    clear1:
                    mol1->releaseOwnSoA();
                    mol2->releaseOwnSoA();
                    container->clear();
                }
            }
        }
    }
}