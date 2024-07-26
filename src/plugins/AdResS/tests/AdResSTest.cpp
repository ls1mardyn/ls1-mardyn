//
// Created by alex on 31.05.23.
//

#include "AdResSTest.h"
#include "plugins/AdResS/AdResS.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"
#include "plugins/AdResS/Interpolation.h"
#include "plugins/AdResS/util/Region.h"
#include "plugins/AdResS/util/WeightFunction.h"
#include "plugins/AdResS/features/Resolution.h"
#include "plugins/AdResS/features/FTH.h"

TEST_SUITE_REGISTRATION(AdResSTest);

AdResSTest::AdResSTest() = default;

AdResSTest::~AdResSTest() = default;

void AdResSTest::computeForcesTest() {
    // std::unique_ptr<AdResS> plugin;
    // ParticleContainer *container = nullptr;
	// Resolution::FPRegion region({4, 4, 4}, {6, 6, 6}, {2, 2, 2});
	// std::vector<Component>* components = nullptr;
    // init
    // {
        // const char *filename = "AdResS-empty-10x10x10.inp";
        // container = initializeFromFile(ParticleContainerFactory::LinkedCell, filename, 2);
// 
        // plugin = std::make_unique<AdResS>();
		// components = _simulation.getEnsemble()->getComponents();
		// Resolution::Config r_conf {
			// components, _simulation.getDomain(), {}, {}
		// };
		// FTH::Config f_conf { false, false, false, false };
        // region.init();
		// r_conf.fpRegions.emplace_back(region);
		// r_conf.comp_to_res.resize(6);
        // for (int i = 0; i < 6; i++) r_conf.comp_to_res[i] = static_cast<Resolution::ResolutionType>(i % 3);
        // plugin->weight = Weight::nearest;
// 
		// plugin->_resolutionHandler.init(r_conf);
    // }
    // ParticlePairsHandler *pairHandler = new AdResSForceAdapter(plugin->_resolutionHandler);
    // CellProcessor *cellProcessor = new LegacyCellProcessor(2, 2, pairHandler);
// 
    // /*
    //  * Cutoff is set to 2 and the domain is 10x10x10 with a FPRegion from 4,4,4 to 6,6,6 with hybrid dims 2,2,2
    //  * */
// 
    // check every single cell -> place 2 molecules in one cell
    // for (int cz = 1; cz <= 9; cz += 2) {
        // for (int cy = 1; cy <= 9; cy += 2) {
            // for (int cx = 1; cx <= 9; cx += 2) {
                // add 2 molecules, set correct component through plugin and let them interact
                // Molecule m1(0, &components->at(2), cx - 0.5, cy, cz, 0, 0, 0, 1, 0, 0, 0);
                // Molecule m2(1, &components->at(2), cx + 0.5, cy, cz, 0, 0, 0, 1, 0, 0, 0);
                // container->addParticle(m1);
                // container->addParticle(m2);
                // plugin->beforeForces(container, &_simulation.domainDecomposition(), 0);
// 
                // get ref back to mols and check the computed forces
                // Molecule *mol1;
                // Molecule *mol2;
                // for (auto it = container->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it) {
                    // if (it->getID() == 0) mol1 = it.operator->();
                    // if (it->getID() == 1) mol2 = it.operator->();
                // }
// 
                // mol1->buildOwnSoA();
                // mol2->buildOwnSoA();
                // ignore non-hybrid interactions
                // if (mol1->componentid() == 0 && 0 == mol2->componentid()) goto clear0;
                // if (mol1->componentid() == 2 && 2 == mol2->componentid()) goto clear0;
// 
                // container->traverseCells(*cellProcessor);
                // if (mol1->componentid() == 0 && mol2->componentid() == 1) {
                    // double weight = plugin->weight(mol2->r_arr(), region);
                    // ASSERT_DOUBLES_EQUAL(-24. * weight, mol1->ljcenter_F(0)[0], 0);
                    // ASSERT_DOUBLES_EQUAL(24. * weight, mol2->ljcenter_F(0)[1], 0);
                // } else if (mol1->componentid() == 1 && mol2->componentid() == 0) {
                    // double weight = plugin->weight(mol1->r_arr(), region);
                    // ASSERT_DOUBLES_EQUAL(-24. * weight, mol1->ljcenter_F(0)[1], 0);
                    // ASSERT_DOUBLES_EQUAL(24. * weight, mol2->ljcenter_F(0)[0], 0);
                // } else {
                    // ASSERT_DOUBLES_EQUAL(-24., mol1->ljcenter_F(0)[0] +
                                            //    ((mol1->componentid() == 1) ? mol1->ljcenter_F(1)[0] : 0), 0);
                    // ASSERT_DOUBLES_EQUAL(24., mol2->ljcenter_F(0)[0] +
                                            //   ((mol2->componentid() == 1) ? mol2->ljcenter_F(1)[0] : 0), 0);
                // }
// 
                // clear0:
                // mol1->releaseOwnSoA();
                // mol2->releaseOwnSoA();
                // container->clear();
            // }
        // }
    // }
// 
    // delete container;
    // delete pairHandler;
    // delete cellProcessor;
}

void AdResSTest::checkGradient() {
    // test function f(x) 2x² - x - 10
    // sample points 0..10..1
    std::vector<double> fun_val;
    std::vector<double> symb_grad;
    fun_val.resize(10);
    symb_grad.resize(10);

    for (int i = 0; i < 10; i++) {
        double x = i + 1;
        double f_x = 2 * std::pow(x, 2) - x - 10;
        fun_val[i] = f_x;

        double df_x = 4 * x - 1;
        symb_grad[i] = df_x;
    }

    std::vector<double> num_grad;
    Interpolation::computeGradient(fun_val, num_grad);

    // we do not care for the border values
    for (int i = 1; i < 9; i++) {
        double rel_error = std::abs(num_grad[i] - symb_grad[i] - 1e-8) / std::abs(symb_grad[i] + 1e-8);
        ASSERT_TRUE(rel_error < 0.05);
    }
}

void AdResSTest::checkMatrixSolver() {
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> x;
    std::vector<double> d{-21, -12, 0, 12, 21};
    std::vector<double> solution{-4.8, -1.8, 0.0, 1.8, 4.8};

    a.resize(5, 1);
    b.resize(5, 4);
    c.resize(5, 1);

    a[0] = 0;
    c[4] = 0;

    Interpolation::solveTriDiagonalMatrix(a, b, c, x, d);
    for (int i = 0; i < 5; i++) {
        ASSERT_DOUBLES_EQUAL(solution[i], x[i], 1e-6);
    }
}

void AdResSTest::checkHermite() {
    double begin = 1.0;
    double end = 11.0; // exclusive
    double step = 0.1;

    std::vector<double> fun_vals;
    fun_vals.resize(100, 0.0);
    for (int i = 0; i < 100; i++) {
        double x = begin + i * step;
        fun_vals[i] = 4 * std::pow(x, 4)
                      - 63 * std::pow(x, 3)
                      + 298 * std::pow(x, 2)
                      - 411 * std::pow(x, 1);
    }

    std::vector<double> x_interpolate{2.0, 2.10126582, 2.20253165, 2.30379747, 2.40506329, 2.50632911,
                                      2.60759494, 2.70886076, 2.81012658, 2.91139241, 3.01265823, 3.11392405,
                                      3.21518987, 3.3164557, 3.41772152, 3.51898734, 3.62025316, 3.72151899,
                                      3.82278481, 3.92405063, 4.02531646, 4.12658228, 4.2278481, 4.32911392,
                                      4.43037975, 4.53164557, 4.63291139, 4.73417722, 4.83544304, 4.93670886,
                                      5.03797468, 5.13924051, 5.24050633, 5.34177215, 5.44303797, 5.5443038,
                                      5.64556962, 5.74683544, 5.84810127, 5.94936709, 6.05063291, 6.15189873,
                                      6.25316456, 6.35443038, 6.4556962, 6.55696203, 6.65822785, 6.75949367,
                                      6.86075949, 6.96202532, 7.06329114, 7.16455696, 7.26582278, 7.36708861,
                                      7.46835443, 7.56962025, 7.67088608, 7.7721519, 7.87341772, 7.97468354,
                                      8.07594937, 8.17721519, 8.27848101, 8.37974684, 8.48101266, 8.58227848,
                                      8.6835443, 8.78481013, 8.88607595, 8.98734177, 9.08860759, 9.18987342,
                                      9.29113924, 9.39240506, 9.49367089, 9.59493671, 9.69620253, 9.79746835,
                                      9.89873418, 10.0};

    Interpolation::Function function;
    std::vector<double> steps;
    steps.resize(99, step);
    Interpolation::computeHermite(begin, fun_vals, steps, 100, function);

    for(int i = 0; i < 80; i++) {
        double x = x_interpolate[i];
        double y_eval = Interpolation::computeHermiteAt(x, function);
        double y_test = 4 * std::pow(x, 4)
                        - 63 * std::pow(x, 3)
                        + 298 * std::pow(x, 2)
                        - 411 * std::pow(x, 1);
        ASSERT_DOUBLES_EQUAL(y_test, y_eval, 1e-4);
    }
}
