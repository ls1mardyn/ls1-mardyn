//
// Created by alex on 31.05.23.
//

#include "AdResSForceAdapterTest.h"
#include "particleContainer/ParticleContainer.h"
#include "plugins/AdResS/AdResS.h"

TEST_SUITE_REGISTRATION(AdResSForceAdapterTest);

AdResSForceAdapterTest::AdResSForceAdapterTest() = default;

AdResSForceAdapterTest::~AdResSForceAdapterTest() = default;

void AdResSForceAdapterTest::processPairTest() {
    /*const char* filename = "AdResS-empty-10x10x10.inp";
    ParticleContainer* container(initializeFromFile(ParticleContainerFactory::LinkedCell, filename, 2));

    std::unique_ptr<AdResS> plugin = std::make_unique<AdResS>();
    plugin->_particleContainer = container;
    FPRegion region({4,4,4}, {6,6,6}, {2,2,2});
    region.init();
    plugin->_fpRegions.emplace_back(region);
    plugin->_components = _simulation.getEnsemble()->getComponents();
    for(int i = 0; i < 6; i++) plugin->_comp_to_res[i] = static_cast<Resolution>(i % 3);
    plugin->_domain = _simulation.getDomain();
    plugin->weight = plugin->weightNearest;

    AdResSForceAdapter* forceAdapter = &plugin->_forceAdapter;
    ParticlePairs2PotForceAdapter pairsHandler(*_simulation.getDomain());

    //check CG-H
    {
        Molecule m1(0, &plugin->_components->at(2), 1.5,3,3, 0,0,0, 1,0,0,0);
        Molecule m2(0, &plugin->_components->at(1), 2.5,3,3, 0,0,0, 1,0,0,0);
        m1.buildOwnSoA();
        m2.buildOwnSoA();
        container->addParticle(m1);
        container->addParticle(m2);

        std::array<double,3> distanceVector = {};
        double dd = m2.dist2(m1, distanceVector.data());
        ASSERT_EQUAL(dd, 1.);

        pairsHandler.init();
        pairsHandler.processPair(m1,m2,distanceVector.data(), MOLECULE_MOLECULE, dd, true);
        pairsHandler.finish();
        forceAdapter->init(_simulation.getDomain());
        forceAdapter->processPair(m1, m2, distanceVector, MOLECULE_MOLECULE, dd, true, true, plugin->_comp_to_res, region);
        forceAdapter->finish();
        ASSERT_DOUBLES_EQUAL(0, m1.ljcenter_F(0)[0], 0);
        ASSERT_DOUBLES_EQUAL(0, m2.ljcenter_F(0)[0] + m2.ljcenter_F(1)[0], 0);

        forceAdapter->init(_simulation.getDomain());
        forceAdapter->processPair(m1, m2, distanceVector, MOLECULE_MOLECULE, dd, true, false, plugin->_comp_to_res, region);
        forceAdapter->finish();
        ASSERT_DOUBLES_EQUAL(-24., m1.ljcenter_F(0)[0], 0);
        ASSERT_DOUBLES_EQUAL(24., m2.ljcenter_F(0)[0] + m2.ljcenter_F(1)[0], 0);
        container->clear();
        m1.releaseOwnSoA();
        m2.releaseOwnSoA();
    }

    //check H-H
    {
        Molecule m1(0, &plugin->_components->at(1), 2.5,3,3, 0,0,0, 1,0,0,0);
        Molecule m2(0, &plugin->_components->at(1), 3.5,3,3, 0,0,0, 1,0,0,0);
        m1.buildOwnSoA();
        m2.buildOwnSoA();
        container->addParticle(m1);
        container->addParticle(m2);

        std::array<double,3> distanceVector = {};
        double dd = m2.dist2(m1, distanceVector.data());
        ASSERT_EQUAL(dd, 1.);

        pairsHandler.init();
        pairsHandler.processPair(m1,m2,distanceVector.data(), MOLECULE_MOLECULE, dd, true);
        pairsHandler.finish();
        forceAdapter->init(_simulation.getDomain());
        forceAdapter->processPair(m1, m2, distanceVector, MOLECULE_MOLECULE, dd, true, true, plugin->_comp_to_res, region);
        forceAdapter->finish();
        ASSERT_DOUBLES_EQUAL(0, m1.ljcenter_F(0)[0]+m1.ljcenter_F(1)[0], 0);
        ASSERT_DOUBLES_EQUAL(0, m2.ljcenter_F(0)[0] + m2.ljcenter_F(1)[0], 0);

        forceAdapter->init(_simulation.getDomain());
        forceAdapter->processPair(m1, m2, distanceVector, MOLECULE_MOLECULE, dd, true, false, plugin->_comp_to_res, region);
        forceAdapter->finish();
        ASSERT_DOUBLES_EQUAL(-24., m1.ljcenter_F(0)[0]+m1.ljcenter_F(1)[0], 0);
        ASSERT_DOUBLES_EQUAL(24., m2.ljcenter_F(0)[0]+m2.ljcenter_F(1)[0], 0);
        container->clear();
        m1.releaseOwnSoA();
        m2.releaseOwnSoA();
    }

    //check H-FP
    {
        Molecule m1(0, &plugin->_components->at(1), 3.5,5,5, 0,0,0, 1,0,0,0);
        Molecule m2(0, &plugin->_components->at(0), 4.5,5,5, 0,0,0, 1,0,0,0);
        m1.buildOwnSoA();
        m2.buildOwnSoA();
        container->addParticle(m1);
        container->addParticle(m2);

        std::array<double,3> distanceVector = {};
        double dd = m2.dist2(m1, distanceVector.data());
        ASSERT_EQUAL(dd, 1.);

        pairsHandler.init();
        pairsHandler.processPair(m1,m2,distanceVector.data(), MOLECULE_MOLECULE, dd, true);
        pairsHandler.finish();
        forceAdapter->init(_simulation.getDomain());
        forceAdapter->processPair(m1, m2, distanceVector, MOLECULE_MOLECULE, dd, true, true, plugin->_comp_to_res, region);
        forceAdapter->finish();
        ASSERT_DOUBLES_EQUAL(0, m1.ljcenter_F(0)[0]+m1.ljcenter_F(1)[0], 0);
        ASSERT_DOUBLES_EQUAL(0, m2.ljcenter_F(0)[0], 0);

        forceAdapter->init(_simulation.getDomain());
        forceAdapter->processPair(m1, m2, distanceVector, MOLECULE_MOLECULE, dd, true, false, plugin->_comp_to_res, region);
        forceAdapter->finish();
        double weight = plugin->weight(m1.r_arr(), region);
        ASSERT_DOUBLES_EQUAL(-24. * weight, m1.ljcenter_F(0)[0]+m1.ljcenter_F(1)[0], 0);
        ASSERT_DOUBLES_EQUAL(24. * weight, m2.ljcenter_F(0)[0], 0);
        container->clear();
        m1.releaseOwnSoA();
        m2.releaseOwnSoA();
    }*/
}