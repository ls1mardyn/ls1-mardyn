/*
 * main.cpp
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */


#include "autopas.h"
#include <iostream>

using namespace std;
using namespace autopas;

class MyMolecule : public Particle {
public:
	MyMolecule(double r = 0.0, double f = 0.0, unsigned long i = 0ul) : Particle(r, f), _id(i) {}
	void print() { cout << "Molecule (pos = " << _x << ", force = " << _f << ", id =" << _id << ")" << endl; }
	unsigned long _id;
};

// how do we recommend to use this thing?
// the user should create a wrapper around our ParticleContainer? That way they can compile whatever they actually need w.r.t. templates once.
// with a bit of documentation, they can

void testFullParticleCell();
void testRMMParticleCell();
void testParticleIteratorFull();
void testParticleIteratorRMM();
void testParticleContainerFull();
void testParticleContainerRMM();

template<class ParticleCell>
void addAFewParticles(ParticleCell& pc) {
	static int i = 0;
	int iEnd = i + 4;
	for ( ; i < iEnd; ++i) {
		MyMolecule m(static_cast<double>(i), static_cast<double>(i), static_cast<unsigned long>(i));
		pc.addParticle(m);
	}
}

int main(void) {
	// test FullParticleCell:
	testFullParticleCell();
	testRMMParticleCell();
	testParticleIteratorFull();
	testParticleIteratorRMM();
	testParticleContainerFull();
	testParticleContainerRMM();

	cout << "Hodor" << endl;
	return EXIT_SUCCESS;
}

void testFullParticleCell() {
	FullParticleCell<MyMolecule> fpc;
	addAFewParticles(fpc);


	SingleCellIterator<MyMolecule, FullParticleCell<MyMolecule>> it(&fpc);
	for (; it.isValid(); ++it) {
		it->print();
	}
}

void testRMMParticleCell() {
	RMMParticleCell<MyMolecule> rmmc;
	addAFewParticles(rmmc);

	cout << "===============" << endl;
	cout << "well, ids are not saved yet: " << endl;
	SingleCellIterator<MyMolecule, RMMParticleCell<MyMolecule>> it(&rmmc);
	for (; it.isValid(); ++it) {
		it->print();
	}
	cout << "===============" << endl;
}

void testParticleIteratorFull() {
	cout << " ============ testing FullCell ============ " << endl;
	std::vector<FullParticleCell<MyMolecule>> vec;
	vec.resize(1);

	cout << "one set" << endl;
	addAFewParticles(vec[0]);
	for (ParticleIterator<MyMolecule, FullParticleCell<MyMolecule>> pi(&vec);
			pi.isValid(); ++pi) {
		pi->print();
	}

	cout << "two sets" << endl;
	vec.resize(2);
	addAFewParticles(vec[1]);

	for (ParticleIterator<MyMolecule, FullParticleCell<MyMolecule>> pi(&vec);
			pi.isValid(); ++pi) {
		pi->print();
	}

	cout << "empty " << endl;
	vec.clear();
	for (ParticleIterator<MyMolecule, FullParticleCell<MyMolecule>> pi(&vec);
			pi.isValid(); ++pi) {
		pi->print();
	}

	cout << "empty full " << endl;
	vec.resize(2);
	addAFewParticles(vec[1]);
	for (ParticleIterator<MyMolecule, FullParticleCell<MyMolecule>> pi(&vec);
			pi.isValid(); ++pi) {
		pi->print();
	}

	cout << "empty full empty" << endl;
	vec.resize(3);
	for (ParticleIterator<MyMolecule, FullParticleCell<MyMolecule>> pi(&vec);
			pi.isValid(); ++pi) {
		pi->print();
	}

	cout << "empty full empty full" << endl;
	vec.resize(4);
	addAFewParticles(vec[3]);
	for (ParticleIterator<MyMolecule, FullParticleCell<MyMolecule>> pi(&vec);
			pi.isValid(); ++pi) {
		pi->print();
	}

	cout << "full empty " << endl;
	vec.clear();
	vec.resize(2);
	addAFewParticles(vec[0]);
	for (ParticleIterator<MyMolecule, FullParticleCell<MyMolecule>> pi(&vec);
			pi.isValid(); ++pi) {
		pi->print();
	}
}

void testParticleIteratorRMM() {
	cout << " ============ testing RMMCell ============ " << endl;
	std::vector<RMMParticleCell<MyMolecule>> vec;
	vec.resize(1);

	cout << "one set" << endl;
	addAFewParticles(vec[0]);
	for (ParticleIterator<MyMolecule, RMMParticleCell<MyMolecule>> pi(&vec);
			pi.isValid(); ++pi) {
		pi->print();
	}

	cout << "two sets" << endl;
	vec.resize(2);
	addAFewParticles(vec[1]);

	for (ParticleIterator<MyMolecule, RMMParticleCell<MyMolecule>> pi(&vec);
			pi.isValid(); ++pi) {
		pi->print();
	}

	cout << "empty " << endl;
	vec.clear();
	for (ParticleIterator<MyMolecule, RMMParticleCell<MyMolecule>> pi(&vec);
			pi.isValid(); ++pi) {
		pi->print();
	}

	cout << "empty full " << endl;
	vec.resize(2);
	addAFewParticles(vec[1]);
	for (ParticleIterator<MyMolecule, RMMParticleCell<MyMolecule>> pi(&vec);
			pi.isValid(); ++pi) {
		pi->print();
	}

	cout << "empty full empty" << endl;
	vec.resize(3);
	for (ParticleIterator<MyMolecule, RMMParticleCell<MyMolecule>> pi(&vec);
			pi.isValid(); ++pi) {
		pi->print();
	}

	cout << "empty full empty full" << endl;
	vec.resize(4);
	addAFewParticles(vec[3]);
	for (ParticleIterator<MyMolecule, RMMParticleCell<MyMolecule>> pi(&vec);
			pi.isValid(); ++pi) {
		pi->print();
	}

	cout << "full empty " << endl;
	vec.clear();
	vec.resize(2);
	addAFewParticles(vec[0]);
	for (ParticleIterator<MyMolecule, RMMParticleCell<MyMolecule>> pi(&vec);
			pi.isValid(); ++pi) {
		pi->print();
	}
}

void testParticleContainerFull() {
	cout << " =========== testing ParticleContainerFull ============" << endl;
	ParticleContainer<MyMolecule, FullParticleCell<MyMolecule>> pc;
	pc.init(5); // empty full empty full empty

	// add a few particles to the second cell
	for (int i = 0; i < 4; ++i) {
		MyMolecule m(static_cast<double>(i), static_cast<double>(i), static_cast<unsigned long>(i));
		pc.addParticle(m, 1);
	}
	// add a few particles to the fourth cell
	for (int i = 4; i < 8; ++i) {
		MyMolecule m(static_cast<double>(i), static_cast<double>(i), static_cast<unsigned long>(i));
		pc.addParticle(m, 3);
	}

	for (auto pi = pc.begin(); pi.isValid(); ++pi) {
		pi->print();
	}

	cout << " =========== done testing ParticleContainerFull ============" << endl;
}

void testParticleContainerRMM() {
	cout << " =========== testing ParticleContainerRMM ============" << endl;
	ParticleContainer<MyMolecule, RMMParticleCell<MyMolecule>> pc;
	pc.init(5); // empty full empty full empty

	// add a few particles to the second cell
	for (int i = 0; i < 4; ++i) {
		MyMolecule m(static_cast<double>(i), static_cast<double>(i), static_cast<unsigned long>(i));
		pc.addParticle(m, 1);
	}
	// add a few particles to the fourth cell
	for (int i = 4; i < 8; ++i) {
		MyMolecule m(static_cast<double>(i), static_cast<double>(i), static_cast<unsigned long>(i));
		pc.addParticle(m, 3);
	}

	for (auto pi = pc.begin(); pi.isValid(); ++pi) {
		pi->print();
	}

	cout << " =========== done testing ParticleContainerRMM ============" << endl;
}



