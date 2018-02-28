/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "NeighbourCommunicationSchemeTest.h"
#define PUSH_PULL_NEIGHBOURS 0


using namespace std;

TEST_SUITE_REGISTRATION(NeighbourCommunicationSchemeTest);


NeighbourCommunicationSchemeTest::NeighbourCommunicationSchemeTest() {
	_fullShell = new FullShell();
	_directScheme = new DirectNeighbourCommunicationScheme(_fullShell);
};

NeighbourCommunicationSchemeTest::~NeighbourCommunicationSchemeTest() {
	delete _fullShell;
	delete _directScheme;
}

#if PUSH_PULL_NEIGHBOURS

void NeighbourCommunicationSchemeTest::testShiftIfNecessary() {
	HaloRegion region; // rmin, rmax, offset, width
	double domainLength[3] = {10.0, 10.0, 10.0};
	double shift[3] = {0.0};
	
	// region does not need to be shifted
	for(int i = 0; i < 3; i++) region.rmax[i] = 5.0;
	for(int i = 0; i < 3; i++) region.rmin[i] = 3.0;
	
	_directScheme->shiftIfNeccessary(domainLength, &region, shift); // region is within domain box
	
	for(int i = 0; i < 3; i++) ASSERT_EQUAL(shift[i], 0.0);
	
	for(int i = 0; i < 3; i++) region.rmax[i] = 12.0;
	for(int i = 0; i < 3; i++) region.rmin[i] = 11.0;
	
	_directScheme->shiftIfNeccessary(domainLength, &region, shift);
	
	for(int i = 0; i < 3; i++) { ASSERT_EQUAL(shift[i], -10.0); shift[i] = 0.0; }
	
	for(int i = 0; i < 3; i++) region.rmax[i] = 0.0;
	for(int i = 0; i < 3; i++) region.rmin[i] = -1.0;
	
	_directScheme->shiftIfNeccessary(domainLength, &region, shift);
	
	for(int i = 0; i < 3; i++) { ASSERT_EQUAL(shift[i], 10.0); shift[i] = 0.0; }
	
	
}

void NeighbourCommunicationSchemeTest::testOverlap() { // assume this one works for now, because you thought about it long and hard.
	HaloRegion region01;
	HaloRegion region02;
	
}

void NeighbourCommunicationSchemeTest::testIOwnThis() { // i own a part of this
	HaloRegion region01;
	HaloRegion region02;
	
}
#endif

