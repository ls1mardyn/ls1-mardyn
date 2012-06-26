/*
 * ParticleInsertionTest.cpp
 *
 *  Created on: Jun 18, 2012
 *      Author: tijana
 */
#include "ParticleInsertionTest.h"
#include "ensemble/PressureGradient.h"
#include <sstream>
TEST_SUITE_REGISTRATION(ParticleInsertionTest);
ParticleInsertionTest::ParticleInsertionTest() {
	// TODO Auto-generated constructor stub

}

ParticleInsertionTest::~ParticleInsertionTest() {
	// TODO Auto-generated destructor stub
}
//FIXME: Correct cutoff for linkedCells
void ParticleInsertionTest::testRotation() {

	// domain bounding box
	double boundings_min[] = { 0, 0, 0 };
	double boundings_max[] = { 33.4278, 33.4278, 33.4278 };

	LinkedCells linkedCells(boundings_min, boundings_max, 100, 100, 3, 1);
	std::vector<Component> components;

	// component - 2 LJ centers
	Component dummyComponent(0);
	dummyComponent.addLJcenter(-1, -1, 0, 0.015, 0.000429325357, 6.61404407, 0,
			false);
	dummyComponent.addLJcenter(1, 1, 0, 0.015, 0.000429325357, 6.61404407, 0,
			false);
	components.push_back(dummyComponent);

	// domain
	int ownrank = 0;
	PressureGradient* pressureGradient = new PressureGradient(ownrank);
	_domain = new Domain(ownrank, pressureGradient);
	_domain->addComponent(dummyComponent);
	_domain->initParameterStreams(100, 100);

	// particle handler for the linked cells
	ParticlePairsHandler* _particlePairsHandler =
			new ParticlePairs2PotForceAdapter(*_domain);
	linkedCells.setPairsHandler(_particlePairsHandler);

	// molecules 1 and 2
	Molecule dummyMolecule1(0, 0, 8, 8, 8, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
			&components);
	dummyMolecule1.upd_cache();
	linkedCells.addParticle(dummyMolecule1);
	linkedCells.update();

	Molecule dummyMolecule2(1, 0, 18, 8, 8, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
			&components);
	dummyMolecule2.upd_cache();
	linkedCells.addParticle(dummyMolecule2);
	linkedCells.update();

	//  ParticleInsertion instance
	moleculardynamics::coupling::ParticleInsertion<Molecule, LinkedCells, 3>
			insertion;

	// molecule to be inserted
	Molecule newMolecule(2, 0, 13, 8, 8, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
			&components);

	int num_tests = 3;
	vector<Quaternion> qs;
	qs.push_back(Quaternion(9.453842e-01, 1.701856e-01, 1.283354e-02,
			2.777065e-01));
	double angle = -PI / 4;
	double torque[3] = { 0, 0, 1 };
	qs.push_back(Quaternion(cos(angle / 2), torque[0] * sin(angle / 2),
			torque[1] * sin(angle / 2), torque[2] * sin(angle / 2)));
	qs.push_back(Quaternion(6.711773e-01, -6.365373e-01, 3.561545e-01,
			1.322697e-01));

	double equals[3] = { 0.944337, 5.01764, 0.944164 };

	double targets[3] = { 0.943416768341188, 5.015456e+00, 0.943416768341188 };
	Simulation* sim = new Simulation();
	moleculardynamics::coupling::interface::MardynMoleculeWrapper<Molecule, 3>
			wrapper(&newMolecule);
	for (int i = 0; i < num_tests; i++) {
		Quaternion q = qs[i];
		newMolecule.setq(q);
		newMolecule.upd_cache();

		// calculate starting energy & force
		double force[3] = { 0, 0, 0 };
		newMolecule.clearFM();

		double energy = linkedCells.getForceAndEnergy(&newMolecule, force);
		cout << "starting energy " << energy << endl;
		tarch::la::Vector<3, double> force_vec(0.0);
		for (int j = 0; j < 3; j++)
			force_vec(j) = force[j];
		double absForce = std::sqrt(force_vec * force_vec);
		wrapper.setForce(force_vec);

		// allocate vectors for logging energy, angle, lj position and center position
		vector<double> vec_energy;
		vector<double> vec_angle;
		vector<double*> vec_lj;
		vector<double*> vec_center;
		string name_energy = "energy.txt";
		string name_angle = "angle.txt";
		string name_lj = "lj.txt";
		string name_center = "center.txt";
		double target_e_2m_xy = 0.943416768341188;
		int timestep = 0;
		double energy_old = energy;

		// perform insertion
		int doInsertion = insertion.rotateMolecule(&newMolecule, wrapper,
				&linkedCells, sim, 1000, 45 * 3.14 / 180, 3.14, 0.1 * PI / 180,
				targets[i], 10000 * newMolecule.getEps(), &timestep, 0.001,
				&energy, &energy_old, &force_vec, &absForce, &q, 0, 1,
				&vec_energy, &vec_angle, &vec_lj, &vec_center, name_energy,
				name_angle, name_lj, name_center);

		ASSERT_DOUBLES_EQUAL(equals[i], energy, 1e-5);
	}

}

void ParticleInsertionTest::testTranslationAndRotation() {
	// domain bounding box
	double boundings_min[] = { 0, 0, 0 };
	double boundings_max[] = { 33.4278, 33.4278, 33.4278 };

	LinkedCells linkedCells(boundings_min, boundings_max, 100, 100, 3, 1);
	std::vector<Component> components;

	// component - 2 LJ centers
	Component dummyComponent(0);
	dummyComponent.addLJcenter(-1, -1, 0, 0.015, 0.000429325357, 6.61404407, 0,
			false);
	dummyComponent.addLJcenter(1, 1, 0, 0.015, 0.000429325357, 6.61404407, 0,
			false);
	components.push_back(dummyComponent);

	// domain
	int ownrank = 0;
	PressureGradient* pressureGradient = new PressureGradient(ownrank);
	_domain = new Domain(ownrank, pressureGradient);
	_domain->addComponent(dummyComponent);
	_domain->initParameterStreams(100, 100);

	// particle handler for the linked cells
	ParticlePairsHandler* _particlePairsHandler =
			new ParticlePairs2PotForceAdapter(*_domain);
	linkedCells.setPairsHandler(_particlePairsHandler);

	// molecules 1 and 2
	Molecule dummyMolecule1(0, 0, 8, 8, 8, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
			&components);
	dummyMolecule1.upd_cache();
	linkedCells.addParticle(dummyMolecule1);
	linkedCells.update();

	Molecule dummyMolecule2(1, 0, 18, 8, 8, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
			&components);
	dummyMolecule2.upd_cache();
	linkedCells.addParticle(dummyMolecule2);
	linkedCells.update();

	//  ParticleInsertion instance
	moleculardynamics::coupling::ParticleInsertion<Molecule, LinkedCells, 3>
			insertion;

	// molecule to be inserted
	Molecule newMolecule(2, 0, 13, 8, 8, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
			&components);

	int seeds[2] = { 4, 5 };
	double equals[2] = { 0.943985, 0.943758 };
	double target = 0.943416768341188;
	int num_tests = 2;
	Simulation* sim = new Simulation();
	ParticleCell cell = linkedCells.getCell(linkedCells.getCellIndexOfMolecule(
			linkedCells.begin()));
	moleculardynamics::coupling::interface::MardynMoleculeWrapper<Molecule, 3>
			wrapper(&newMolecule);

	for (int i = 0; i < num_tests; i++) {
		// calculate starting energy & force
		double force[3] = { 0, 0, 0 };
		newMolecule.clearFM();
		Simulation* sim = new Simulation();
		double energy = linkedCells.getForceAndEnergy(&newMolecule, force);
		cout << "starting energy " << energy << endl;
		tarch::la::Vector<3, double> force_vec(0.0);
		for (int j = 0; j < 3; j++)
			force_vec(j) = force[j];
		double absForce = std::sqrt(force_vec * force_vec);
		wrapper.setForce(force_vec);

		// allocate for logging
		vector<double> vec_energy;
		vector<double> vec_angle;
		vector<double*> vec_lj;
		vector<double*> vec_center;
		string name_energy = "energy.txt";
		string name_angle = "angle.txt";
		string name_lj = "lj.txt";
		string name_center = "center.txt";
		double target = 0.943416768341188;
		int timestep = 0;
		double energy_old = energy;

insertion		.findParticlePosition(
				&linkedCells,
				cell,
				tarch::la::Vector<3,unsigned int>(0,0,0),
				tarch::la::Vector<3,unsigned int>(0,0,0),
				tarch::la::Vector<3,unsigned int>(0,0,0),
				&newMolecule, sim, target, &energy, &energy_old, 0, 1, seeds[i],
				100, 100, 100, 45 * PI / 180, 3.14, 0.1 * PI / 180, 0.001,
				&vec_energy, &vec_angle, &vec_lj, &vec_center, name_energy,
				name_angle, name_lj, name_center);
		ASSERT_DOUBLES_EQUAL(equals[i], energy, 1e-5);

	}
}
void ParticleInsertionTest::testParameterSetup() {

	LinkedCells* linkedCells = (LinkedCells*) initializeFromFile(
			ParticleContainerFactory::LinkedCell, "Ethan_equilibrated.inp",
			32.1254);

	int ownrank = 0;
	PressureGradient* pressureGradient = new PressureGradient(ownrank);
	std::vector<Component> components;

	// component - 2 LJ centers
	Component dummyComponent(0);
	dummyComponent.addLJcenter(0, 0, -2.2157048, 0.0150347,	0.00042932536, 6.6140441, 0, 0);
	dummyComponent.addLJcenter(0, 0, 2.2157048,	0.0150347,	0.00042932536, 6.6140441, 0, 0);
	components.push_back(dummyComponent);
	_domain = new Domain(ownrank, pressureGradient);
	_domain->addComponent(dummyComponent);
	_domain->initParameterStreams(32.1254, 32.1254);
	moleculardynamics::coupling::ParticleInsertion<Molecule, LinkedCells, 3>
			insertion;

	// molecule to be inserted
	Molecule newMolecule(linkedCells->getNumberOfParticles() + 1, 0, 13, 8, 8, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
			&components);

	// particle handler for the linked cells
	ParticlePairsHandler* _particlePairsHandler =
			new ParticlePairs2PotForceAdapter(*_domain);
	linkedCells->setPairsHandler(_particlePairsHandler);
	int num_tests = 17;
	Simulation* sim = new Simulation();
	ParticleCell cell = linkedCells->getCell(
			linkedCells->getCellIndexOfMolecule(linkedCells->begin()));
	moleculardynamics::coupling::interface::MardynMoleculeWrapper<Molecule, 3>
			wrapper(&newMolecule);
	double u_avg = 0;
	int count = 0;
	double* force = new double[3];
	for (Molecule* molec = linkedCells->begin(); molec != linkedCells->end(); molec
			= linkedCells->next()) {
		count++;
		u_avg += linkedCells->getForceAndEnergy(molec, force);
	}

	u_avg /= count;
	cout << "target energy: " << u_avg << endl;
	string file_name_base = "../test_input/param_config";
	for (int i = 0; i < num_tests; i++) {
		stringstream file_stream;
		file_stream << file_name_base << "_" << i + 1 << ".txt"<<endl;
		string file_name;
		file_stream >> file_name;
		// calculate starting energy & force
		double force[3] = { 0, 0, 0 };
		newMolecule.clearFM();
		Simulation* sim = new Simulation();
		double energy = linkedCells->getForceAndEnergy(&newMolecule, force);

		tarch::la::Vector<3, double> force_vec(0.0);
		for (int j = 0; j < 3; j++)
			force_vec(j) = force[j];
		double absForce = std::sqrt(force_vec * force_vec);
		wrapper.setForce(force_vec);

		// allocate for logging
		vector<double> vec_energy;
		vector<double> vec_angle;
		vector<double*> vec_lj;
		vector<double*> vec_center;
		string name_energy = "energy.txt";
		string name_angle = "angle.txt";
		string name_lj = "lj.txt";
		string name_center = "center.txt";
		double target = 0.943416768341188;
		int timestep = 0;
		double energy_old = energy;
		cout << "starting insertion" << endl;
		int* maxIter = new int[1];
		int* maxRestarts = new int[1];
		int* maxRotations = new int[1];
		double* tolerance = new double[1];
		double* maxAngle = new double[1];
		double* maxAllowedAngle = new double[1];
		double* minAngle = new double[1];
		bool* largeStepsizeOnOverlap = new bool[1];
		bool* restartIfIncreases = new bool[1];
		readParamFile(file_name, maxIter, maxRestarts, maxRotations, tolerance,
				maxAngle, maxAllowedAngle, minAngle, largeStepsizeOnOverlap,
				restartIfIncreases);
insertion	.findParticlePosition(
			linkedCells,
			cell,
			tarch::la::Vector<3,unsigned int>(0,0,0),
			tarch::la::Vector<3,unsigned int>(0,0,0),
			tarch::la::Vector<3,unsigned int>(0,0,0),
			&newMolecule, sim, u_avg, &energy, &energy_old, *largeStepsizeOnOverlap, *restartIfIncreases, 5,
			*maxIter, *maxRestarts, *maxRotations, *maxAngle, *maxAllowedAngle, *minAngle, *tolerance,
			&vec_energy, &vec_angle, &vec_lj, &vec_center, name_energy,
			name_angle, name_lj, name_center);
		delete maxIter;
		delete maxRestarts;
		delete maxRotations;
		delete maxAngle;
		delete maxAllowedAngle;
		delete minAngle;
		delete tolerance;
		delete largeStepsizeOnOverlap;
		delete restartIfIncreases;
	//ASSERT_DOUBLES_EQUAL(equals[i], energy, 1e-5);
	}
}
void ParticleInsertionTest::readParamFile(string file_name, int* maxIter, int* maxRestarts, int* maxRotations, double* tolerance,
		double* maxAngle, double* maxAllowedAngle, double* minAngle, bool* largeStepsizeOnOverlap,
		bool* restartIfIncreases) {
	FILE* file = NULL;
	file = fopen(file_name.c_str(), "r");
	if (file == NULL) cout<<"problem with the file!!"<<endl;

	int* restartIfIncreasesInt = new int[1];
	int* largeStepsizeOnOverlapInt = new int[1];

	fscanf(file, "%d\n%d\n%d\n%lf\n%lf\n%lf\n%lf\n%d\n%d", maxIter, maxRestarts, maxRotations, tolerance,
			 maxAngle, maxAllowedAngle, minAngle, largeStepsizeOnOverlapInt, restartIfIncreasesInt);
	if (*restartIfIncreasesInt == 1) {
		*restartIfIncreases = true;
	} else if (*restartIfIncreasesInt == 0) {
		*restartIfIncreases = false;
	} else {
		cout<<"parameter not specified correctly"<<endl;
	}
	if (*largeStepsizeOnOverlapInt == 1) {
		*largeStepsizeOnOverlap = true;
	} else if (*largeStepsizeOnOverlapInt == 0) {
		*largeStepsizeOnOverlap= false;
	} else {
		cout<<"parameter not specified correctly"<<endl;
	}


	cout<<"Chosen configuration for filename: "<< file_name << endl;
	printf( "%d\n%d\n%d\n%g\n%g\n%g\n%g\n%d\n%d\n", *maxIter, *maxRestarts, *maxRotations,
			*tolerance, *maxAngle, *maxAllowedAngle, *minAngle, *largeStepsizeOnOverlapInt, *restartIfIncreasesInt);
	delete largeStepsizeOnOverlapInt;
	delete restartIfIncreasesInt;
	fclose(file);
}


