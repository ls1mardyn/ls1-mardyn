/*
 * Adios2IOTest.cpp
 *
 * Check whether an ADIOS2 checkpoint can be successfully written and read again.
 *
 *  Created on: July 2021
 *      Author: Heinen, Gralka, Rau
 */
#ifdef ENABLE_ADIOS2

#include "Domain.h"
#include "molecules/Quaternion.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/LinkedCells.h"
#include "parallel/DomainDecompBase.h"
#include <iostream>
#include <random>
#include <filesystem>
#include <numeric>

#include "io/tests/Adios2IOTest.h"
#include "io/Adios2Writer.h"
#include "io/Adios2Reader.h"

#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#endif

#define NUM_PARTICLES 5



#if !defined(ENABLE_REDUCED_MEMORY_MODE)
TEST_SUITE_REGISTRATION(Adios2IOTest);
#else
#pragma message "Compilation Info: Adios2IOTest disabled in reduced memory mode."
#endif


/*
 * Initialize test simulation
 */
void Adios2IOTest::initParticles() {

	_ids.resize(NUM_PARTICLES);
	_positions.resize(NUM_PARTICLES);
	_velocities.resize(NUM_PARTICLES);
	_Ds.resize(NUM_PARTICLES);
	_quaternions.resize(NUM_PARTICLES);

	// init bbox
	_box_lower = {0,0,0};
	_box_upper = {10,10,10};

	// init ids
	std::iota(std::begin(_ids), std::end(_ids), 0);
	
	// init positions
    std::mt19937 rnd;
    rnd.seed(666);
    std::uniform_real_distribution<double> dist_pos(_box_lower[0], _box_upper[0]);
    for (auto& pos : _positions) {
		for (int i = 0; i < 3; ++i) {
        	pos[i] = dist_pos(rnd);
		}
    }

	// init component
	_comp.resize(1);
	_comp[0].setName("1CLJ");
	_comp[0].addLJcenter(0,0,0,1,1,1);

	// init velocities
	std::uniform_real_distribution<double> dist_v(-1, 1);
    for (auto& v : _velocities) {
		for (int i = 0; i < 3; ++i) {
        	v[i] = dist_v(rnd);
		}
    }

	// init angular momentum
    for (auto& D : _Ds) {
		for (int i = 0; i < 3; ++i) {
        	D[i] = dist_v(rnd);
		}
    }

	// init quaternions
	for (auto& q : _quaternions) {
			q[0] = 1.0;
			q[1] = 1.0;
			q[2] = 0.0;
			q[3] = 0.0;
	}
}


/*
 * Actual test if a checkpoint can successfully be written.
 */
void Adios2IOTest::testWriteCheckpoint() {
	initParticles();

	auto particleContainer = std::make_shared<LinkedCells>(_box_lower.data(), _box_upper.data(), _cutoff);

	std::vector<Molecule> particles(NUM_PARTICLES);
	for (int i = 0; i < particles.size(); ++i) {
		particles[i].setid(_ids[i]);
		particles[i].setComponent(&_comp[0]);
		particles[i].setr(0,_positions[i][0]);
		particles[i].setr(1,_positions[i][1]);
		particles[i].setr(2,_positions[i][2]);
		particles[i].setv(0,_velocities[i][0]);
		particles[i].setv(1,_velocities[i][1]);
		particles[i].setv(2,_velocities[i][2]);
		Quaternion q(_quaternions[i][0], _quaternions[i][1], _quaternions[i][2], _quaternions[i][3]);
		particles[i].setq(q);
		particles[i].setD(0,_Ds[i][0]);
		particles[i].setD(1,_Ds[i][1]);
		particles[i].setD(2,_Ds[i][2]);
	}
	particleContainer->addParticles(particles);

	
	auto adios2writer = std::make_shared<Adios2Writer>();
	adios2writer->init(nullptr, nullptr, nullptr);
	adios2writer->testInit(_comp ,_filename);
	int ownrank = 0;
	std::shared_ptr<DomainDecompBase> domaindecomp;
#ifdef ENABLE_MPI
	MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &ownrank) );
	global_log->info() << "[Adios2IOTest] Creating standard domain decomposition ... " << endl;
	domaindecomp = std::make_shared<DomainDecomposition>();
#else
	global_log->info() << "[Adios2IOTest] Creating alibi domain decomposition ... " << endl;
	domaindecomp = std::make_shared<DomainDecompBase>();
#endif

	global_log->info() << "[Adios2IOTest] ParticleContainer num particles:"
					   << particleContainer->getNumberOfParticles() << std::endl;
	
	ASSERT_EQUAL(static_cast<unsigned long>(NUM_PARTICLES), particleContainer->getNumberOfParticles());
	
	auto domain = std::make_shared<Domain>(ownrank);
	domain->setGlobalLength(0, _box_upper[0]);
	domain->setGlobalLength(1, _box_upper[1]);
	domain->setGlobalLength(2, _box_upper[2]);
	domain->updateglobalNumMolecules(particleContainer.get(), domaindecomp.get());
	adios2writer->endStep(particleContainer.get(), domaindecomp.get(), domain.get(), 0);
	adios2writer->finish(particleContainer.get(), domaindecomp.get(), domain.get());
	
	ASSERT_EQUAL(true, std::filesystem::is_directory(_filename));
}


/*
 * Actual test if a written checkpoint can successfully be read again.
 */
void Adios2IOTest::testReadCheckpoint() {
	initParticles();
	
	_inputPatricleContainer =
		std::make_shared<LinkedCells>(_box_lower.data(), _box_upper.data(), _cutoff);
	_inputDomain = std::make_shared<Domain>(0);
	_inputDomain->setGlobalLength(0, _box_upper[0]);
	_inputDomain->setGlobalLength(1, _box_upper[1]);
	_inputDomain->setGlobalLength(2, _box_upper[2]);
	_inputDomainDecomp = std::make_shared<DomainDecompBase>();

	auto adios2read = std::make_shared<Adios2Reader>();
	adios2read->testInit(_filename);
	unsigned long pcount = 0;
	try {
		pcount =
			adios2read->readPhaseSpace(_inputPatricleContainer.get(), _inputDomain.get(), _inputDomainDecomp.get());
	} catch (const std::exception& e) {
		global_log->error() << "[Adios2IOTest] exception: " << e.what() << std::endl;
	}
	
	ASSERT_EQUAL(static_cast<unsigned long>(NUM_PARTICLES), pcount);
	
	for (auto it = _inputPatricleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it) {
		auto i  = it->getID();
		for (int j = 0; j < 3; ++j) {
			//global_log->error() << "[Adios2IOTest] position container: " << it->r(j) << " vector: " << _positions[i][j]
			//					<< std::endl;
			
			ASSERT_EQUAL(it->r(j), _positions[i][j]);
			ASSERT_EQUAL(it->D(j), _Ds[i][j]);
			ASSERT_EQUAL(it->v(j), _velocities[i][j]);
			ASSERT_EQUAL(it->getID(), _ids[i]);
		}
	}
}

#endif