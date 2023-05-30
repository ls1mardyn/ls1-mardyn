//
// Created by alex on 30.05.23.
//

#include "AdResSTest.h"
#include "plugins/AdResS/AdResS.h"

TEST_SUITE_REGISTRATION(AdResSTest);

AdResSTest::AdResSTest() = default;

AdResSTest::~AdResSTest() = default;

void AdResSTest::computeForcesTest() {
    std::unique_ptr<AdResS> plugin;
    std::unique_ptr<AdResSForceAdapter> forceAdapter;
    FPRegion region({4,4,4}, {6,6,6}, {2,2,2});
    //init
    {
        const char* filename = "AdResS-empty-10x10x10.inp";
        ParticleContainer* container0(initializeFromFile(ParticleContainerFactory::LinkedCell, filename, 2));
        std::vector<Component> compsBackup;
        Comp2Param paramBackup(_simulation.getDomain()->getComp2Params());
        for(Component& c : *_simulation.getEnsemble()->getComponents()) compsBackup.emplace_back(c);
        ParticleContainer* container1(initializeFromFile(ParticleContainerFactory::LinkedCell, filename, 2));
        ParticleContainer* container2(initializeFromFile(ParticleContainerFactory::LinkedCell, filename, 2));
        _simulation.getEnsemble()->getComponents()->clear();
        for(Component& c : compsBackup) _simulation.getEnsemble()->getComponents()->emplace_back(c);
        _simulation.getDomain()->getComp2Params() = paramBackup;

        _simulation.getMoleculeContainers().push_back(container0);
        _simulation.getMoleculeContainers().push_back(container1);
        _simulation.getMoleculeContainers().push_back(container2);

        plugin = std::make_unique<AdResS>();

        region.init();
        plugin->_fpRegions.emplace_back(region);
        plugin->_components = _simulation.getEnsemble()->getComponents();
        for(int i = 0; i < 6; i++) plugin->_comp_to_res[i] = static_cast<Resolution>(i % 3);
        plugin->_domain = _simulation.getDomain();
        plugin->weight = plugin->weightNearest;

        forceAdapter = std::make_unique<AdResSForceAdapter>(*plugin);
        plugin->_forceAdapter = forceAdapter.get();
    }

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
                plugin->_particleContainers[2]->addParticle(m1);
                plugin->_particleContainers[2]->addParticle(m2);
                for(ParticleContainer* ptr : plugin->_particleContainers) {
                    plugin->beforeForces(ptr, nullptr, 0);
                }
                Molecule* mol1;
                Molecule* mol2;
                for(ParticleContainer* ptr : plugin->_particleContainers) {
                    for(auto it = ptr->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it){
                        if(it->getID() == 0) mol1 = it.operator->();
                        if(it->getID() == 1) mol2 = it.operator->();
                    }
                }

                mol1->buildOwnSoA();
                mol2->buildOwnSoA();
                if(mol1->componentid() == mol2->componentid()) goto clear0;

                plugin->siteWiseForces(plugin->_particleContainers[FullParticle], nullptr, 0);
                if(mol1->componentid() == 0 && mol2->componentid() == 1) {
                    double weight = plugin->weight(mol2->r_arr(), region);
                    ASSERT_DOUBLES_EQUAL(-24. * weight, mol1->ljcenter_F(0)[0], 0);
                    ASSERT_DOUBLES_EQUAL(24. * weight, mol2->ljcenter_F(0)[1], 0);
                } else if( mol1->componentid() == 1 && mol2->componentid() == 0) {
                    double weight = plugin->weight(mol1->r_arr(), region);
                    ASSERT_DOUBLES_EQUAL(-24. * weight, mol1->ljcenter_F(0)[1], 0);
                    ASSERT_DOUBLES_EQUAL(24. * weight, mol2->ljcenter_F(0)[0], 0);
                } else {
                    ASSERT_DOUBLES_EQUAL(-24., mol1->ljcenter_F(0)[0] + mol1->componentid() == 1 ? mol1->ljcenter_F(1)[0] : 0, 0);
                    ASSERT_DOUBLES_EQUAL(24., mol2->ljcenter_F(0)[0] + mol2->componentid() == 1 ? mol2->ljcenter_F(1)[0] : 0, 0);
                }

                clear0:
                mol1->releaseOwnSoA();
                mol2->releaseOwnSoA();
                for(ParticleContainer* ptr : plugin->_particleContainers) {
                    ptr->clear();
                }
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
                    plugin->_particleContainers[2]->addParticle(m1);
                    plugin->_particleContainers[2]->addParticle(m2);
                    for(ParticleContainer* ptr : plugin->_particleContainers) {
                        plugin->beforeForces(ptr, nullptr, 0);
                    }
                    Molecule* mol1;
                    Molecule* mol2;
                    for(ParticleContainer* ptr : plugin->_particleContainers) {
                        for(auto it = ptr->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it){
                            if(it->getID() == 0) mol1 = it.operator->();
                            if(it->getID() == 1) mol2 = it.operator->();
                        }
                    }

                    mol1->buildOwnSoA();
                    mol2->buildOwnSoA();
                    if(mol1->componentid() == mol2->componentid()) goto clear1;

                    plugin->siteWiseForces(plugin->_particleContainers[FullParticle], nullptr, 0);
                    if(mol1->componentid() == 1 && mol2->componentid() == 1) {
                        auto f1 = mol1->ljcenter_F(0);
                        auto f2 = mol2->ljcenter_F(0);
                        for(int d = 0; d<3; d++) {
                            f1[d] += mol1->ljcenter_F(1)[d];
                            f2[d] += mol2->ljcenter_F(1)[d];
                            ASSERT_TRUE(f2[d] == -f1[d]);
                        }
                        ASSERT_TRUE(f1[0] != 0 || f1[1] != 0 || f1[2] != 0);
                    } else if(mol1->componentid() == 1 && mol2->componentid() != 1) {
                        auto f1 = mol1->ljcenter_F(0);
                        auto f2 = mol2->ljcenter_F(0);
                        for(int d = 0; d<3; d++) {
                            f1[d] += mol1->ljcenter_F(1)[d];
                            ASSERT_TRUE(f2[d] == -f1[d]);
                        }
                        ASSERT_TRUE(f1[0] != 0 || f1[1] != 0 || f1[2] != 0);
                    } else if(mol1->componentid() != 1 && mol2->componentid() == 1) {
                        auto f1 = mol1->ljcenter_F(0);
                        auto f2 = mol2->ljcenter_F(0);
                        for(int d = 0; d<3; d++) {
                            f2[d] += mol2->ljcenter_F(1)[d];
                            ASSERT_TRUE(f2[d] == -f1[d]);
                        }
                        if(f1[0] == 0 && f1[1] == 0 && f1[2] == 0)
                            std::cout << "yeet";
                        ASSERT_TRUE(f1[0] != 0 || f1[1] != 0 || f1[2] != 0);
                    } else {
                        auto f1 = mol1->ljcenter_F(0);
                        auto f2 = mol2->ljcenter_F(0);
                        for(int d = 0; d<3; d++) {
                            ASSERT_TRUE(f2[d] == -f1[d]);
                        }
                        ASSERT_TRUE(f1[0] != 0 || f1[1] != 0 || f1[2] != 0);
                    }

                    clear1:
                    mol1->releaseOwnSoA();
                    mol2->releaseOwnSoA();
                    for(ParticleContainer* ptr : plugin->_particleContainers) {
                        ptr->clear();
                    }
                }
            }
        }
    }
}
