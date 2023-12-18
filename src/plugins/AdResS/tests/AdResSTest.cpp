//
// Created by alex on 31.05.23.
//

#include "AdResSTest.h"
#include "plugins/AdResS/AdResS.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"

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
    ParticlePairsHandler* pairHandler = new AdResSForceAdapter(*plugin);
    CellProcessor* cellProcessor = new LegacyCellProcessor(2, 2, pairHandler);

    /*
     * Cutoff is set to 2 and the domain is 10x10x10 with a FPRegion from 4,4,4 to 6,6,6 with hybrid dims 2,2,2
     * */

    //check every single cell -> place 2 molecules in one cell
    for(int cz = 1; cz <= 9; cz+=2) {
        for (int cy = 1; cy <= 9; cy+=2) {
            for (int cx = 1; cx <= 9; cx+=2) {
                //add 2 molecules, set correct component through plugin and let them interact
                Molecule m1(0, &plugin->_components->at(2), cx-0.5,cy,cz, 0,0,0, 1,0,0,0);
                Molecule m2(1, &plugin->_components->at(2), cx+0.5,cy,cz, 0,0,0, 1,0,0,0);
                container->addParticle(m1);
                container->addParticle(m2);
                plugin->beforeForces(container, &_simulation.domainDecomposition(), 0);

                //get ref back to mols and check the computed forces
                Molecule* mol1;
                Molecule* mol2;
                for(auto it = container->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it){
                    if(it->getID() == 0) mol1 = it.operator->();
                    if(it->getID() == 1) mol2 = it.operator->();
                }

                mol1->buildOwnSoA();
                mol2->buildOwnSoA();
                //ignore non-hybrid interactions
                if(mol1->componentid() == 0 && 0 == mol2->componentid()) goto clear0;
                if(mol1->componentid() == 2 && 2 == mol2->componentid()) goto clear0;

                container->traverseCells(*cellProcessor);
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

    delete container;
    delete pairHandler;
    delete cellProcessor;
}

void AdResSTest::checkGradient() {
    std::unique_ptr<AdResS> plugin = std::make_unique<AdResS>();
    // test function f(x) 2x² - x - 10
    // sample points 0..10..1
    std::vector<double> fun_val;
    std::vector<double> symb_grad;
    fun_val.resize(10);
    symb_grad.resize(10);

    for(int i = 0; i < 10; i++) {
        double x = i + 1;
        double f_x = 2 * std::pow(x, 2) - x - 10;
        fun_val[i] = f_x;

        double df_x = 4 * x - 1;
        symb_grad[i] = df_x;
    }

    std::vector<double> num_grad;
    plugin->computeGradient(fun_val, num_grad);

    // we do not care for the border values
    for(int i = 1; i < 9; i++) {
        double rel_error = std::abs(num_grad[i] - symb_grad[i] - 1e-8) / std::abs(symb_grad[i] + 1e-8);
        ASSERT_TRUE(rel_error < 0.05);
    }
}

void AdResSTest::checkMatrixSolver() {
    std::unique_ptr<AdResS> plugin = std::make_unique<AdResS>();
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> x;
    std::vector<double> d {-21, -12, 0, 12, 21};
    std::vector<double> solution {-4.8, -1.8, 0.0, 1.8, 4.8};

    a.resize(5, 1);
    b.resize(5, 4);
    c.resize(5, 1);

    a[0] = 0;
    c[4] = 0;

    plugin->solveTriDiagonalMatrix(a, b, c, x, d);
    for(int i = 0; i < 5; i++) {
        ASSERT_DOUBLES_EQUAL(solution[i], x[i], 1e-6);
    }
}
