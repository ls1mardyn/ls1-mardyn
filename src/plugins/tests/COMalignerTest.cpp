//
// Created by kruegener on 5/31/2018.
//

#include "COMalignerTest.h"

TEST_SUITE_REGISTRATION(COMalignerTest);

COMalignerTest::COMalignerTest() {}

COMalignerTest::~COMalignerTest() {}

void COMalignerTest::testCOMalign() {

    if (_domainDecomposition->getNumProcs() >= 10){
        test_log -> info() << "COMalignerTest::testCOMalign: SKIPPED (required fewer than 10 processes but was run with " << _domainDecomposition->getNumProcs() << " => bounding box of test setup is too small to support decomposition)" << std::endl;
        return;
    }

    const char* filename = "1clj-regular-2x2x2-offset.inp";
    double cutoff = .5;
    ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, filename, cutoff);
    ParticleContainer* oldContainer = initializeFromFile(ParticleContainerFactory::LinkedCell, filename, cutoff);

    COMaligner* plugin = new COMaligner();

    plugin->init(container, _domainDecomposition, _domain);
    plugin->beforeForces(container, _domainDecomposition, 1);

    // TEST MASS
    // TODO: Mass check doesnt work with multiple Particles with "World Record" engaged
    double m = plugin->_mass;
    if (_domainDecomposition->getNumProcs() != 1) {
        test_log->info() << "COMalignerTest::testCOMalign: Mass Check SKIPPED (required exactly 1 process but was run with " <<  _domainDecomposition->getNumProcs() << " processes)" << std::endl;
    }
    else{
		double expectedMass; // 8.0 in Normal mode, 16.0 in RMM mode, cause of different behaviour of Molecule::mass() !
#ifdef ENABLE_REDUCED_MEMORY_MODE
		expectedMass = 16.0;
#else
		expectedMass = 8.0;
#endif
		ASSERT_EQUAL_MSG("Mass does not match number of particles", expectedMass, m);
    }
    // Hard Coded 8.0 instead
    //ASSERT_EQUAL_MSG("Mass does not match number of particles", 8.0, m);

    // TEST MOTION
    ASSERT_EQUAL_MSG("x motion is wrong", -.25, plugin->_motion[0]);
    ASSERT_EQUAL_MSG("y motion is wrong", -.25, plugin->_motion[1]);
    ASSERT_EQUAL_MSG("z motion is wrong", -.25, plugin->_motion[2]);

    // TEST IF MOTION WAS APPLIED
    ParticleIterator newPos = container->iterator();
    ParticleIterator oldPos = oldContainer->iterator();
    while(newPos.hasNext()){
        for(int d = 0; d < 3; d++){
            ASSERT_EQUAL_MSG("Motion has not been properly applied" ,oldPos->r(d) - .25, newPos->r(d));
        }
        newPos.next();
        oldPos.next();
    }

}
