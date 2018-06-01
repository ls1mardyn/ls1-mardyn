//
// Created by kruegener on 5/31/2018.
//

#include "COMalignerTest.h"

TEST_SUITE_REGISTRATION(COMalignerTest);

COMalignerTest::COMalignerTest() {}

COMalignerTest::~COMalignerTest() {}

void COMalignerTest::testCOMalign() {

    const char* filename = "1clj-regular-2x2x2-offset.inp";
    double cutoff = 1.5;
    ParticleContainer* container = initializeFromFile(ParticleContainerFactory::LinkedCell, filename, cutoff);
    ParticleContainer* oldContainer = initializeFromFile(ParticleContainerFactory::LinkedCell, filename, cutoff);

    COMaligner* plugin = new COMaligner();

    plugin->init(container, _domainDecomposition, _domain);
    plugin->beforeForces(container, _domainDecomposition, 1);

    // TEST MASS
    double m = plugin->_mass;
    ASSERT_EQUAL_MSG("Mass does not match number of particles", double(container->getNumberOfParticles()), m);

    // TEST MOTION
    ASSERT_EQUAL_MSG("x motion is wrong", -.5, plugin->_motion[0]);
    ASSERT_EQUAL_MSG("y motion is wrong", -.5, plugin->_motion[1]);
    ASSERT_EQUAL_MSG("z motion is wrong", -.5, plugin->_motion[2]);

    // TEST IF MOTION WAS APPLIED
    ParticleIterator newPos = container->iterator();
    ParticleIterator oldPos = oldContainer->iterator();
    while(newPos.hasNext()){
        for(int d = 0; d < 3; d++){
            ASSERT_EQUAL_MSG("Motion has not been properly applied" ,oldPos->r(d) - .5, newPos->r(d));
        }
        newPos.next();
        oldPos.next();
    }

}