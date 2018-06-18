/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "NeighbourCommunicationSchemeTest.h"


using namespace std;

TEST_SUITE_REGISTRATION(NeighbourCommunicationSchemeTest);


NeighbourCommunicationSchemeTest::NeighbourCommunicationSchemeTest() {
	_fullShell = new FullShell();
	_directScheme = new DirectNeighbourCommunicationScheme(_fullShell, true);
}

NeighbourCommunicationSchemeTest::~NeighbourCommunicationSchemeTest() {
	//delete _fullShell;
	delete _directScheme;
}


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
	
	for(int i = 0; i < 3; i++) { 
		region01.rmax[i] = 4.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 6.0;
		region02.rmin[i] = 3.0;
	}
	
	_directScheme->overlap(&region01, &region02);
	
	for(int i = 0; i < 3; i++) {
		ASSERT_EQUAL(region02.rmax[i], 4.0);
		ASSERT_EQUAL(region02.rmin[i], 3.0);
	}
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 6.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 5.0;
		region02.rmin[i] = 3.0;
	}
	
	_directScheme->overlap(&region01, &region02);
	
	for(int i = 0; i < 3; i++) {
		ASSERT_EQUAL(region02.rmax[i], 5.0);
		ASSERT_EQUAL(region02.rmin[i], 3.0);
	}
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 4.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 3.0;
		region02.rmin[i] = 1.0;
	}
	
	_directScheme->overlap(&region01, &region02);
	
	for(int i = 0; i < 3; i++) {
		ASSERT_EQUAL(region02.rmax[i], 3.0);
		ASSERT_EQUAL(region02.rmin[i], 2.0);
	}
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 4.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 6.0;
		region02.rmin[i] = 1.0;
	}
	
	_directScheme->overlap(&region01, &region02);
	
	for(int i = 0; i < 3; i++) {
		ASSERT_EQUAL(region02.rmax[i], 4.0);
		ASSERT_EQUAL(region02.rmin[i], 2.0);
	}
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 6.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 6.0;
		region02.rmin[i] = 2.0;
	}
	
	_directScheme->overlap(&region01, &region02);
	
	for(int i = 0; i < 3; i++) {
		ASSERT_EQUAL(region02.rmax[i], 6.0);
		ASSERT_EQUAL(region02.rmin[i], 2.0);
	}
}

void NeighbourCommunicationSchemeTest::testIOwnThis() { // i own a part of this
	HaloRegion region01;
	HaloRegion region02;
	
	for(int i = 0; i < 3; i++) { 
		region01.rmax[i] = 4.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 6.0;
		region02.rmin[i] = 3.0;
	}
	
	ASSERT_EQUAL(_directScheme->iOwnThis(&region01, &region02), true);
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 6.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 5.0;
		region02.rmin[i] = 3.0;
	}
	
	ASSERT_EQUAL(_directScheme->iOwnThis(&region01, &region02), true);
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 4.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 3.0;
		region02.rmin[i] = 1.0;
	}
	
	ASSERT_EQUAL(_directScheme->iOwnThis(&region01, &region02), true);
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 4.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 6.0;
		region02.rmin[i] = 1.0;
	}
	
	ASSERT_EQUAL(_directScheme->iOwnThis(&region01, &region02), true); 
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 6.0;
		region01.rmin[i] = 2.0;
		region02.rmax[i] = 6.0;
		region02.rmin[i] = 2.0;
	}
	
	ASSERT_EQUAL(_directScheme->iOwnThis(&region01, &region02), true); 
	
	for(int i = 0; i < 3; i++) {
		region01.rmax[i] = 5.0;
		region01.rmin[i] = 1.0;
	}
	
	region02.rmax[0] = 6.0;
	region02.rmax[1] = 7.0;
	region02.rmax[2] = 5.0;
	
	region02.rmin[0] = 5.0;
	region02.rmin[1] = 6.0;
	region02.rmin[2] = 5.0;
	
	ASSERT_EQUAL(_directScheme->iOwnThis(&region01, &region02), false);
	
	
}
