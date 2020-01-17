/**
 * @file AutoPasContainerTest.cpp
 * @author seckler
 * @date 19.09.18
 */

#include "AutoPasContainerTest.h"
#include "particleContainer/AutoPasContainer.h"

TEST_SUITE_REGISTRATION(AutoPasContainerTest);

void AutoPasContainerTest::testConstructor() { AutoPasContainer container(1.); }

void AutoPasContainerTest::testUpdate() {
	AutoPasContainer autoPasContainer(1.);
	double min[3] = {0., 0., 0.};
	double max[3] = {10., 10., 10.};
	autoPasContainer.rebuild(min, max);
	int id = 1;
	for (double x : {0., 5., 9.999}) {
		for (double y : {0., 5., 9.999}) {
			for (double z : {0., 5., 9.999}) {
				Molecule p(id++, nullptr, x, y, z, 0., 0., 0.);

				autoPasContainer.addParticle(p);       // inside, therefore ok!
			}
		}
	}
	std::set<unsigned long> movedIDs;
	// we move particles that are close to the boundary to outside of the container and remember the id's we moved
	for(auto iter = autoPasContainer.iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); iter.isValid(); ++iter){
		for(unsigned short dim = 0; dim < 3; ++dim){
			if(iter->getR()[dim] < 0.5){
				auto r = iter->getR();
				// smallest double smaller than 0
				r[dim] = std::nexttoward(0., -1.);
				iter->setR(r);
				movedIDs.insert(iter->getID());
			}
			if(iter->getR()[dim] > 9.5){
				auto r = iter->getR();
				r[dim] = 10.;
				iter->setR(r);
				movedIDs.insert(iter->getID());
			}
		}
	}

	// now update the container!
	autoPasContainer.update();

	// the particles should no longer be in the inner cells!
	for(auto iter = autoPasContainer.iterator(ParticleIterator::Type::ONLY_INNER_AND_BOUNDARY); iter.isValid(); ++iter){
		ASSERT_EQUAL(movedIDs.count(iter->getID()), 0ul);
	}
}
