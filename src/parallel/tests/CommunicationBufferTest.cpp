/*
 * CommunicationBufferTest.cpp
 *
 *  Created on: 16 Oct 2017
 *      Author: tchipevn
 */

#include "CommunicationBufferTest.h"

#include "parallel/CommunicationBuffer.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "parallel/DomainDecomposition.h"
#include "Simulation.h"
#include "ensemble/EnsembleBase.h"
#include "ensemble/CanonicalEnsemble.h"

#include <iostream>

TEST_SUITE_REGISTRATION(CommunicationBufferTest);

CommunicationBufferTest::CommunicationBufferTest() {
}

CommunicationBufferTest::~CommunicationBufferTest() {
	delete global_simulation;
}

void CommunicationBufferTest::testEmplaceRead() {
	CommunicationBuffer buf;
	buf.resizeForRawBytes(1000);

	float f1 = 22. / 7., f2 = 13. / 7.6;
	double d1 = 17. / 187.;
	int i1 = -129387;
	unsigned long ul1 = 203482020348;
	size_t ret = 0;
	ret = buf.emplaceValue<float>(ret, f1);
	ASSERT_EQUAL(ret, 4ul);

	ret = buf.emplaceValue<float>(ret, f2);
	ASSERT_EQUAL(ret, 8ul);

	ret = buf.emplaceValue<double>(ret, d1);
	ASSERT_EQUAL(ret, 16ul);

	ret = buf.emplaceValue<int>(ret, i1);
	ASSERT_EQUAL(ret, 20ul);

	ret = buf.emplaceValue<unsigned long>(ret, ul1);
	ASSERT_EQUAL(ret, 28ul);

	ret = buf.emplaceValue<float>(ret, f1);
	ASSERT_EQUAL(ret, 32ul);

	float fread1, fread2;
	double dread1;
	int iread1;
	unsigned long ulread1;
	ret = 0;

	ret = buf.readValue<float>(ret, fread1);
	ASSERT_EQUAL(ret, 4ul);
	ASSERT_DOUBLES_EQUAL(f1, fread1, 1e-16);

	ret = buf.readValue<float>(ret, fread2);
	ASSERT_EQUAL(ret, 8ul);
	ASSERT_DOUBLES_EQUAL(f2, fread2, 1e-16);

	ret = buf.readValue<double>(ret, dread1);
	ASSERT_EQUAL(ret, 16ul);
	ASSERT_DOUBLES_EQUAL(d1, dread1, 1e-16);

	ret = buf.readValue<int>(ret, iread1);
	ASSERT_EQUAL(ret, 20ul);
	ASSERT_EQUAL(i1, iread1);

	ret = buf.readValue<unsigned long>(ret, ulread1);
	ASSERT_EQUAL(ret, 28ul);
	ASSERT_EQUAL(ul1, ulread1);

	ret = buf.readValue<float>(ret, fread1);
	ASSERT_EQUAL(ret, 32ul);
	ASSERT_DOUBLES_EQUAL(f1, fread1, 1e-16);
}

void CommunicationBufferTest::testHalo() {
	Component dummyComponent(0);
	dummyComponent.addLJcenter(0, 0, 0, 1, 1, 1, 0, false);
	global_simulation->getEnsemble()->addComponent(dummyComponent);

	CommunicationBuffer buf;
	buf.resizeForAppendingLeavingMolecules(0);
	buf.resizeForAppendingHaloMolecules(2);

	unsigned long uid = 0;
	Molecule m[2];
	m[0] = Molecule(0, global_simulation->getEnsemble()->getComponent(0), 1.,
			2., 3., -1., -2., -3.);
	m[1] = Molecule(1, global_simulation->getEnsemble()->getComponent(0), 11.,
			12., 13., -11., -12., -13.);
	buf.addHaloMolecule(0, m[0]);
	buf.addHaloMolecule(1, m[1]);

	Molecule mread[2];
	buf.readHaloMolecule(0, mread[0]);
	buf.readHaloMolecule(1, mread[1]);

	for (int i = 0; i < 2; ++i) {
#ifdef LS1_SEND_UNIQUE_ID_FOR_HALO_COPIES
		ASSERT_EQUAL(m[i].getID(), mread[i].getID());
#else
		ASSERT_EQUAL(UINT64_MAX, mread[i].getID());
#endif
		for (int d = 0; d < 3; ++d) {
			ASSERT_DOUBLES_EQUAL(m[i].r(d), mread[i].r(d), 1e-16);
			ASSERT_DOUBLES_EQUAL(0.0, mread[i].v(d), 1e-16);
		}
	}
}

void CommunicationBufferTest::testLeaving() {
	Component dummyComponent(0);
	dummyComponent.addLJcenter(0, 0, 0, 1, 1, 1, 0, false);
	global_simulation->getEnsemble()->addComponent(dummyComponent);

	CommunicationBuffer buf;
	buf.resizeForAppendingLeavingMolecules(2);

	unsigned long uid = 0;
	Molecule m[2];
	m[0] = Molecule(0, global_simulation->getEnsemble()->getComponent(0), 1.,
			2., 3., -1., -2., -3.);
	m[1] = Molecule(1, global_simulation->getEnsemble()->getComponent(0), 11.,
			12., 13., -11., -12., -13.);
	buf.addLeavingMolecule(0, m[0]);
	buf.addLeavingMolecule(1, m[1]);

	Molecule mread[2];
	buf.readLeavingMolecule(0, mread[0]);
	buf.readLeavingMolecule(1, mread[1]);

	for (int i = 0; i < 2; ++i) {
		ASSERT_EQUAL(m[i].getID(), mread[i].getID());
		for (int d = 0; d < 3; ++d) {
			ASSERT_DOUBLES_EQUAL(m[i].r(d), mread[i].r(d), 1e-16);
			ASSERT_DOUBLES_EQUAL(m[i].v(d), mread[i].v(d), 1e-16);
		}
	}
}

void CommunicationBufferTest::testLeavingAndHalo() {
	Component dummyComponent(0);
	dummyComponent.addLJcenter(0, 0, 0, 1, 1, 1, 0, false);
	global_simulation->getEnsemble()->addComponent(dummyComponent);

	Molecule m[5];
	m[0] = Molecule(0, global_simulation->getEnsemble()->getComponent(0), 1.,
			2., 3., -1., -2., -3.);
	m[1] = Molecule(1, global_simulation->getEnsemble()->getComponent(0), 11.,
			12., 13., -11., -12., -13.);
	m[2] = Molecule(2, global_simulation->getEnsemble()->getComponent(0), 21.,
			22., 23., -21., -22., -23.);
	m[3] = Molecule(3, global_simulation->getEnsemble()->getComponent(0), 31.,
			32., 33., -31., -32., -33.);
	m[4] = Molecule(4, global_simulation->getEnsemble()->getComponent(0), 41.,
			42., 43., -41., -42., -43.);

	CommunicationBuffer buf;
	buf.resizeForAppendingLeavingMolecules(2);
	buf.addLeavingMolecule(0, m[0]);
	buf.addLeavingMolecule(1, m[1]);

	buf.resizeForAppendingHaloMolecules(3);
	buf.addHaloMolecule(0, m[2]);
	buf.addHaloMolecule(1, m[3]);
	buf.addHaloMolecule(2, m[4]);

	Molecule mread[5];
	buf.readLeavingMolecule(0, mread[0]);
	buf.readLeavingMolecule(1, mread[1]);
	buf.readHaloMolecule(0, mread[2]);
	buf.readHaloMolecule(1, mread[3]);
	buf.readHaloMolecule(2, mread[4]);

	for (int i = 0; i < 5; ++i) {
		if (i < 2) {  // leaving
			ASSERT_EQUAL(m[i].getID(), mread[i].getID());
		} else {  // halo
#ifdef LS1_SEND_UNIQUE_ID_FOR_HALO_COPIES
		ASSERT_EQUAL(m[i].getID(), mread[i].getID());
#else
			ASSERT_EQUAL(UINT64_MAX, mread[i].getID());
#endif
		}
		for (int d = 0; d < 3; ++d) {
			ASSERT_DOUBLES_EQUAL(m[i].r(d), mread[i].r(d), 1e-16);
			if (i < 2) {  // leaving
				ASSERT_DOUBLES_EQUAL(m[i].v(d), mread[i].v(d), 1e-16);
			} else {  // halo
				ASSERT_DOUBLES_EQUAL(0.0, mread[i].v(d), 1e-16);
			}
		}
	}
}

void CommunicationBufferTest::testPackSendRecvUnpack() {
	if (_domainDecomposition->getNumProcs() < 2) {
		test_log->info() << "CommunicationBufferTest::testPackSendRecvUnpack"
				<< " not executed (rerun with more than 2 Processes!)"
				<< std::endl;
		std::cout << "numProcs:" << _domainDecomposition->getNumProcs()
				<< std::endl;
		return;
	}
	Component dummyComponent(0);
	dummyComponent.addLJcenter(0, 0, 0, 1, 1, 1, 0, false);
	global_simulation->getEnsemble()->addComponent(dummyComponent);

	Molecule m[5];
	m[0] = Molecule(0, global_simulation->getEnsemble()->getComponent(0), 1.,
			2., 3., -1., -2., -3.);
	m[1] = Molecule(1, global_simulation->getEnsemble()->getComponent(0), 11.,
			12., 13., -11., -12., -13.);
	m[2] = Molecule(2, global_simulation->getEnsemble()->getComponent(0), 21.,
			22., 23., -21., -22., -23.);
	m[3] = Molecule(3, global_simulation->getEnsemble()->getComponent(0), 31.,
			32., 33., -31., -32., -33.);
	m[4] = Molecule(4, global_simulation->getEnsemble()->getComponent(0), 41.,
			42., 43., -41., -42., -43.);

	int rank = _domainDecomposition->getRank();
	int Tag = 17;

	// rank 0 packs and sends
	if (rank == 0) {
		CommunicationBuffer buf;
		buf.resizeForAppendingLeavingMolecules(2);

		buf.addLeavingMolecule(0, m[0]);
		buf.addLeavingMolecule(1, m[1]);

		buf.resizeForAppendingHaloMolecules(3);
		buf.addHaloMolecule(0, m[2]);
		buf.addHaloMolecule(1, m[3]);
		buf.addHaloMolecule(2, m[4]);

		// send of course
		MPI_Send(buf.getDataForSending(), buf.getNumElementsForSending(),
				buf.getMPIDataType(), 1, Tag, MPI_COMM_WORLD);
	} else if (rank == 1) {
		// probe
		MPI_Status stat;
		MPI_Probe(0, Tag, MPI_COMM_WORLD, &stat);

		// get the message length
		CommunicationBuffer buf;
		int messageLengthInBytes;
		MPI_Get_count(&stat, buf.getMPIDataType(), &messageLengthInBytes);

		// recv
		buf.resizeForRawBytes(messageLengthInBytes);
		MPI_Recv(buf.getDataForSending(), messageLengthInBytes,
				buf.getMPIDataType(), 0, Tag, MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

		// start reading
		unsigned long numLeaving, numHalo;
		buf.resizeForReceivingMolecules(numLeaving, numHalo);

		Molecule mread[5];

		// check contents
		for (unsigned long i = 0; i < numLeaving; ++i) {
			buf.readLeavingMolecule(i, mread[i]);
		}
		for (unsigned long i = 0; i < numHalo; ++i) {
			buf.readHaloMolecule(i, mread[numLeaving + i]);
		}

		for (int i = 0; i < 5; ++i) {
			if (i < 2) {  // leaving
				ASSERT_EQUAL(m[i].getID(), mread[i].getID());
			} else {  // halo
#ifdef LS1_SEND_UNIQUE_ID_FOR_HALO_COPIES
				ASSERT_EQUAL(m[i].getID(), mread[i].getID());
#else
				ASSERT_EQUAL(UINT64_MAX, mread[i].getID());
#endif
			}
			for (int d = 0; d < 3; ++d) {
				ASSERT_DOUBLES_EQUAL(m[i].r(d), mread[i].r(d), 1e-16);
				if (i < 2) {
					ASSERT_DOUBLES_EQUAL(m[i].v(d), mread[i].v(d), 1e-16);
				} else {
					ASSERT_DOUBLES_EQUAL(0.0, mread[i].v(d), 1e-16);
				}
			}
		}
	}

	//MPI_Barrier(MPI_COMM_WORLD);
}
