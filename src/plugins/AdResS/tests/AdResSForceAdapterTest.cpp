//
// Created by alex on 28.05.23.
//

#include "AdResSForceAdapterTest.h"
#include "particleContainer/ParticleContainer.h"
#include "plugins/AdResS/AdResS.h"

TEST_SUITE_REGISTRATION(AdResSForceAdapterTest);

AdResSForceAdapterTest::AdResSForceAdapterTest() = default;

AdResSForceAdapterTest::~AdResSForceAdapterTest() = default;

void AdResSForceAdapterTest::processPairTest() {
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

    std::unique_ptr<AdResS> plugin = std::make_unique<AdResS>();
    FPRegion region({4,4,4}, {6,6,6}, {2,2,2});
    region.init();
    plugin->_fpRegions.emplace_back(region);
    plugin->_components = _simulation.getEnsemble()->getComponents();
    for(int i = 0; i < 6; i++) plugin->_comp_to_res[i] = static_cast<Resolution>(i % 3);
    plugin->_domain = _simulation.getDomain();
    plugin->weight = plugin->weightNearest;

    std::unique_ptr<AdResSForceAdapter> forceAdapter = std::make_unique<AdResSForceAdapter>(*plugin);
    plugin->_forceAdapter = forceAdapter.get();

    //check CG-CG
    {
        Molecule m1(0, &plugin->_components->at(2), 0.5,1,1, 0,0,0, 1,0,0,0);
        Molecule m2(0, &plugin->_components->at(2), 1.5,1,1, 0,0,0, 1,0,0,0);
        m1.buildOwnSoA();
        m2.buildOwnSoA();
        container2->addParticle(m1);
        container2->addParticle(m2);

        double distanceVector[3];
        double dd = m2.dist2(m1, distanceVector);
        ASSERT_EQUAL(dd, 1.);
        forceAdapter->init();
        forceAdapter->processPair(m1, m2, distanceVector, MOLECULE_MOLECULE, dd, true);
        forceAdapter->finish();
        ASSERT_DOUBLES_EQUAL(-24., m1.ljcenter_F(0)[0], 0);
        ASSERT_DOUBLES_EQUAL(24., m2.ljcenter_F(0)[0], 0);
        container2->clear();
        m1.releaseOwnSoA();
        m2.releaseOwnSoA();
    }

    //check CG-H
    {
        Molecule m1(0, &plugin->_components->at(2), 1.5,3,3, 0,0,0, 1,0,0,0);
        Molecule m2(0, &plugin->_components->at(1), 2.5,3,3, 0,0,0, 1,0,0,0);
        m1.buildOwnSoA();
        m2.buildOwnSoA();
        container2->addParticle(m1);
        container2->addParticle(m2);

        double distanceVector[3];
        double dd = m2.dist2(m1, distanceVector);
        ASSERT_EQUAL(dd, 1.);
        forceAdapter->init();
        forceAdapter->processPair(m1, m2, distanceVector, MOLECULE_MOLECULE, dd, true);
        forceAdapter->finish();
        ASSERT_DOUBLES_EQUAL(-24., m1.ljcenter_F(0)[0], 0);
        ASSERT_DOUBLES_EQUAL(24., m2.ljcenter_F(0)[0] + m2.ljcenter_F(1)[0], 0);
        container2->clear();
        m1.releaseOwnSoA();
        m2.releaseOwnSoA();
    }

    //check H-H
    {
        Molecule m1(0, &plugin->_components->at(1), 2.5,3,3, 0,0,0, 1,0,0,0);
        Molecule m2(0, &plugin->_components->at(1), 3.5,3,3, 0,0,0, 1,0,0,0);
        m1.buildOwnSoA();
        m2.buildOwnSoA();
        container2->addParticle(m1);
        container2->addParticle(m2);

        double distanceVector[3];
        double dd = m2.dist2(m1, distanceVector);
        ASSERT_EQUAL(dd, 1.);
        forceAdapter->init();
        forceAdapter->processPair(m1, m2, distanceVector, MOLECULE_MOLECULE, dd, true);
        forceAdapter->finish();
        ASSERT_DOUBLES_EQUAL(-24., m1.ljcenter_F(0)[0]+m1.ljcenter_F(1)[0], 0);
        ASSERT_DOUBLES_EQUAL(24., m2.ljcenter_F(0)[0]+m2.ljcenter_F(1)[0], 0);
        container2->clear();
        m1.releaseOwnSoA();
        m2.releaseOwnSoA();
    }

    //check H-FP
    {
        Molecule m1(0, &plugin->_components->at(1), 3.5,5,5, 0,0,0, 1,0,0,0);
        Molecule m2(0, &plugin->_components->at(0), 4.5,5,5, 0,0,0, 1,0,0,0);
        m1.buildOwnSoA();
        m2.buildOwnSoA();
        container2->addParticle(m1);
        container2->addParticle(m2);

        double distanceVector[3];
        double dd = m2.dist2(m1, distanceVector);
        ASSERT_EQUAL(dd, 1.);
        forceAdapter->init();
        forceAdapter->processPair(m1, m2, distanceVector, MOLECULE_MOLECULE, dd, true);
        forceAdapter->finish();
        double weight = plugin->weight(m1.r_arr(), region);
        ASSERT_DOUBLES_EQUAL(-24. * weight, m1.ljcenter_F(0)[0]+m1.ljcenter_F(1)[0], 0);
        ASSERT_DOUBLES_EQUAL(24. * weight, m2.ljcenter_F(0)[0], 0);
        container2->clear();
        m1.releaseOwnSoA();
        m2.releaseOwnSoA();
    }

    //check FP-FP
    {
        Molecule m1(0, &plugin->_components->at(0), 4.5,5,5, 0,0,0, 1,0,0,0);
        Molecule m2(0, &plugin->_components->at(0), 5.5,5,5, 0,0,0, 1,0,0,0);
        m1.buildOwnSoA();
        m2.buildOwnSoA();
        container2->addParticle(m1);
        container2->addParticle(m2);

        double distanceVector[3];
        double dd = m2.dist2(m1, distanceVector);
        ASSERT_EQUAL(dd, 1.);
        forceAdapter->init();
        forceAdapter->processPair(m1, m2, distanceVector, MOLECULE_MOLECULE, dd, true);
        forceAdapter->finish();
        ASSERT_DOUBLES_EQUAL(-24., m1.ljcenter_F(0)[0], 0);
        ASSERT_DOUBLES_EQUAL(24., m2.ljcenter_F(0)[0], 0);
        container2->clear();
        m1.releaseOwnSoA();
        m2.releaseOwnSoA();
    }
}
