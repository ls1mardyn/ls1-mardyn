/*
 * RDFDummyDecomposition.cpp
 *
 * Domain decomposition to be used if RDF boundaries are used.
 * Makes sure no halo layers are created
 *
 *  Created on: Jul 5, 2012
 *      Author: tijana
 */

#include "RDFDummyDecomposition.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "utils/Logger.h"
#include "particleContainer/LinkedCells.h"
#include "ParticleInsertion.h"
#include "Simulation.h"

using namespace std;
using Log::global_log;

bool RDFDummyDecomposition::have_avg_energy = false;
int RDFDummyDecomposition::num_calles = 0;
bool RDFDummyDecomposition::first_unif = true;
double* RDFDummyDecomposition::unif_rand = new double[2];

RDFDummyDecomposition::RDFDummyDecomposition(ParticlePairsHandler* ph,
		int boundary, int insertion_type, Simulation* sim) :
	DomainDecompDummy(), _particlePairsHandler(ph), rdfBoundary(boundary),
			particle_insertion_type(insertion_type), simulation(sim) {
	if (particle_insertion_type == 1)
		std::cout << "chosen insertion type is usehr";
	else if (particle_insertion_type == 0)
		std::cout << "chosen insertion type is reflective wall";
	else
		std::cout << "unknown insertion type" << std::endl;
	// TODO Auto-generated constructor stub

}

double RDFDummyDecomposition::getUniformRandomNumber() {

	return ((double) rand()) / ((double) RAND_MAX);
}

double RDFDummyDecomposition::getGaussianRandomNumber() {

	if (first_unif) {

		// in this case: generate new numbers
		double s = 2.0;
		double v[2] = { 0, 0 };

		while (s >= 1.0) {
			unif_rand[0] = ((double) rand()) / ((double) RAND_MAX);
			unif_rand[1] = ((double) rand()) / ((double) RAND_MAX);
			v[0] = (2.0 * unif_rand[0] - 1.0);
			v[1] = (2.0 * unif_rand[1] - 1.0);
			s = v[0] * v[0] + v[1] * v[1];
		}
		unif_rand[0] = v[0] * sqrt((-2.0 * log(s)) / s);
		unif_rand[1] = v[1] * sqrt((-2.0 * log(s)) / s);

		// change to the other variable in the next call
		first_unif = false;
		return unif_rand[0];

		// otherwise: change to the other random number
	} else {
		first_unif = true;
		return unif_rand[1];
	}
}

double RDFDummyDecomposition::getAverageEnergy(LinkedCells* linkedCells,
		double* low_limit, double* high_limit) {
	double u_avg = 0;
	double* force = new double[3]; // placeholder
	int num = 0;
	Molecule* currentMolecule;
	std::cout << "calling for: " << low_limit[0] << " " << low_limit[1] << " "
			<< low_limit[2] << std::endl;
	std::cout << "high: " << high_limit[0] << " " << high_limit[1] << " "
			<< high_limit[2] << std::endl;
	for (currentMolecule = linkedCells->begin(); currentMolecule
			!= linkedCells->end(); currentMolecule = linkedCells->next()) {
		if (currentMolecule->r(0) < high_limit[0] && currentMolecule->r(0)
				> low_limit[0] && currentMolecule->r(1) < high_limit[1]
				&& currentMolecule->r(1) > low_limit[1]
				&& currentMolecule->r(2) < high_limit[2] && currentMolecule->r(
				2) > low_limit[2]) {
			u_avg += linkedCells->getEnergy(currentMolecule, force);
			num++;
		}
	}
	std::cout << "num: " << num << std::endl;
	if (num > 0)
		u_avg /= num;

	return u_avg;
}

void RDFDummyDecomposition::generateRandomVelocity(double temperature,
		double mass, double* v) {

	// Velocity mardyn code
	//	for (int dim = 0; dim < 3; dim++) {
	//		v[dim] = randdouble(-0.5, 0.5);
	//	}
	//	double dotprod_v = 0;
	//	for (unsigned int i = 0; i < 3; i++) {
	//		dotprod_v += v[i] * v[i];
	//	}
	//	// Velocity Correction
	//	double vCorr = sqrt(3.0 * temperature / dotprod_v);
	//	for (unsigned int i = 0; i < 3; i++) {
	//		v[i] *= vCorr;
	//	}

	double randomNumbers[2] = { 0, 0 };
	double kB = 1;
	// standard deviation for fluctuation of velocity around meanVelocity
	double stdDeviation = std::sqrt(3 * kB * temperature / mass);
	double meanVelocity = 0;//std::sqrt(2 * kB * temperature / mass);

	// random number (unit mean
	randomNumbers[0] = getGaussianRandomNumber();
	// then, get D-1 uniform random numbers over 2PI
	for (unsigned int d = 1; d < 3; d++) {
		randomNumbers[d] = 2.0 * 3.14 * getUniformRandomNumber();
	}

	// put initial velocity together. For 2D/ 3D, we use the polar/ spherical coordinates to generate the fluctuation around the mean velocity

	v[0] = meanVelocity + stdDeviation * (randomNumbers[0] * std::sin(
			randomNumbers[1]) * std::cos(randomNumbers[2]));
	v[1] = meanVelocity + stdDeviation * (randomNumbers[0] * std::sin(
			randomNumbers[1]) * std::sin(randomNumbers[2]));
	v[2] = meanVelocity + stdDeviation * (randomNumbers[0] * std::cos(
			randomNumbers[1]));

}

void RDFDummyDecomposition::generateRandomAngularVelocity(double temperature,
		double* w, Domain* domain, Molecule* currentMolecule) {

	std::vector<Component> components = domain->getComponents();
	double I[3] = { 0., 0., 0. };
	I[0] = components[currentMolecule->componentid()].I11();
	I[1] = components[currentMolecule->componentid()].I22();
	I[2] = components[currentMolecule->componentid()].I33();
	for (int d = 0; d < 3; d++) {
		w[d] = (I[d] == 0) ? 0.0 : ((randdouble(0, 1) > 0.5) ? 1 : -1) * sqrt(
				2.0 * randdouble(0, 1) * temperature / I[d]);
		w[d] = w[d] * 0.030619994; // fs_2_mardyn
	}
}
void RDFDummyDecomposition::balanceAndExchange(bool balance,
		ParticleContainer* moleculeContainer,
		const std::vector<Component>& components, Domain* domain) {
	exchangeMolecules(moleculeContainer, components, domain);

}

RDFDummyDecomposition::~RDFDummyDecomposition() {
	// TODO Auto-generated destructor stub
}

void RDFDummyDecomposition::addPeriodicCopies(
		ParticleContainer* moleculeContainer, double* rmin, double* rmax,
		double* phaseSpaceSize, double* halo_L,
		const std::vector<Component>& components) {
	Molecule* currentMolecule = moleculeContainer->begin();
	double low_limit[3] = { 0, 0, 0 };
	double high_limit[3] = { 0, 0, 0 };
	double new_position[3] = { 0, 0, 0 };
	// now traverse everything that's not the rdf boundary
	for (unsigned short d = 0; d < 3; ++d) {
		if (d == rdfBoundary)
			continue;

		phaseSpaceSize[d] = rmax[d] - rmin[d];
		// set limits (outside "inner" region)
		low_limit[d] = rmin[d] + halo_L[d];
		high_limit[d] = rmax[d] - halo_L[d];
		currentMolecule = moleculeContainer->begin();

		while (currentMolecule != moleculeContainer->end()) {
			//if wall, enable reflections
			if (!particle_insertion_type)
				currentMolecule->enableBouncingBack(rmin[rdfBoundary],
						rmax[rdfBoundary], 0, 0);
			const double& rd = currentMolecule->r(d);
			if (currentMolecule->r(rdfBoundary) < rmin[rdfBoundary]
					|| currentMolecule->r(rdfBoundary) > rmax[rdfBoundary]) {
				currentMolecule = moleculeContainer->next();
				continue;
			}
			if (rd < low_limit[d]) {
				// determine the position for the copy of the molecule
				for (unsigned short d2 = 0; d2 < 3; d2++) {
					// when moving parallel to the coordinate d2 to another process, the
					// local coordinates in d2 change
					if (d2 == d)
						new_position[d2] = rd + phaseSpaceSize[d2];
					else
						new_position[d2] = currentMolecule->r(d2);
				}

				Molecule m1 = Molecule(currentMolecule->id(),
						currentMolecule->componentid(), new_position[0],
						new_position[1], new_position[2],
						currentMolecule->v(0), currentMolecule->v(1),
						currentMolecule->v(2), currentMolecule->q().qw(),
						currentMolecule->q().qx(), currentMolecule->q().qy(),
						currentMolecule->q().qz(), currentMolecule->D(0),
						currentMolecule->D(1), currentMolecule->D(2),
						&components);
				// if wall, enable reflections
				if (!particle_insertion_type)
					m1.enableBouncingBack(rmin[rdfBoundary], rmax[rdfBoundary],
							0, 0);
				moleculeContainer->addParticle(m1);
				currentMolecule = moleculeContainer->next();
			} else if (rd >= high_limit[d]) {
				// determine the position for the copy of the molecule
				for (unsigned short d2 = 0; d2 < 3; d2++) {
					// when moving parallel to the coordinate d2 to another process, the
					// local coordinates in d2 change
					if (d2 == d)
						new_position[d2] = rd - phaseSpaceSize[d2];
					else
						new_position[d2] = currentMolecule->r(d2);
				}
				Molecule m1 = Molecule(currentMolecule->id(),
						currentMolecule->componentid(), new_position[0],
						new_position[1], new_position[2],
						currentMolecule->v(0), currentMolecule->v(1),
						currentMolecule->v(2), currentMolecule->q().qw(),
						currentMolecule->q().qx(), currentMolecule->q().qy(),
						currentMolecule->q().qz(), currentMolecule->D(0),
						currentMolecule->D(1), currentMolecule->D(2),
						&components);
				moleculeContainer->addParticle(m1);

				// if wall enable reflections
				if (!particle_insertion_type)
					m1.enableBouncingBack(rmin[rdfBoundary], rmax[rdfBoundary],
							0, 0);
				currentMolecule = moleculeContainer->next();
			} else
				currentMolecule = moleculeContainer->next();
		}
	}
}

void RDFDummyDecomposition::exchangeMolecules(
		ParticleContainer* moleculeContainer,
		const std::vector<Component>& components, Domain* domain) {

	double rmin[3]; // lower corner of the process-specific domain //ENABLE_MPI
	double rmax[3];
	double halo_L[3]; // width of the halo strip //ENABLE_MPI
	for (int i = 0; i < 3; i++) {
		rmin[i] = moleculeContainer->getBoundingBoxMin(i);
		rmax[i] = moleculeContainer->getBoundingBoxMax(i);
		halo_L[i] = moleculeContainer->get_halo_L(i);
	}

	double phaseSpaceSize[3];

	double low_limit[3]; // particles below this limit have to be copied or moved to the lower process
	double high_limit[3]; // particles above(or equal) this limit have to be copied or moved to the higher process

	std::vector<Molecule*> molecules_to_delete;
	vector<double*> starting_positions;
	vector<int> deleted_ids;
	vector<int> deleted_comp_ids;

	for (unsigned short d = 0; d < 3; ++d) {

		phaseSpaceSize[d] = rmax[d] - rmin[d];
		// set limits (outside "inner" region)
		low_limit[d] = rmin[d] + halo_L[d];
		high_limit[d] = rmax[d] - halo_L[d];
	}

	LinkedCells* linkedCells = (LinkedCells*) moleculeContainer;
	linkedCells->setPairsHandler(_particlePairsHandler);

	// average energy for usher

	addPeriodicCopies(moleculeContainer, rmin, rmax, phaseSpaceSize, halo_L,
			components);
	moleculeContainer->update();

	num_calles++;

}
void RDFDummyDecomposition::insertUsher(ParticleContainer* moleculeContainer,
		const std::vector<Component>& components, Domain* domain) {
	if (num_calles <= 1)
		return;
	double rmin[3]; // lower corner of the process-specific domain //ENABLE_MPI
	double rmax[3];
	double phaseSpaceSize[3];
	double halo_L[3]; // width of the halo strip //ENABLE_MPI
	for (int i = 0; i < 3; i++) {
		rmin[i] = moleculeContainer->getBoundingBoxMin(i);
		rmax[i] = moleculeContainer->getBoundingBoxMax(i);
		halo_L[i] = moleculeContainer->get_halo_L(i);
		phaseSpaceSize[i] = rmax[i] - rmin[i];
	}

	Molecule* molecule;
	int deleted_id, deleted_comp_id, i = 1;
	double v_old[3];
	double w_old[3];
	double r_old[3];
	molecule = moleculeContainer->begin();
	double rc = moleculeContainer->getCutoff();

	LinkedCells* linkedCells = (LinkedCells*) moleculeContainer;
	Molecule* m;

	if (num_calles > 1) {
		double rmin_inner[3] = { rc, rc, rc };
		double rmax_inner[3] = { rmax[0] - rc, rmax[1] - rc, rmax[2] - rc };
		u_avg = this->getAverageEnergy(linkedCells, rmin_inner, rmax_inner);
		have_avg_energy = true;
	}

	double temperature = domain->getCurrentTemperature(0);
	double v[3] = { 0, 0, 0 };
	double w[3] = { 0, 0, 0 };
	this->generateRandomVelocity(temperature, molecule->mass(), v);
	this->generateRandomAngularVelocity(temperature, w, domain, molecule);

	while (molecule != moleculeContainer->end()) {

		if (molecule->r(rdfBoundary) < rmin[rdfBoundary] || molecule->r(
				rdfBoundary) > rmax[rdfBoundary]) {
			deleted_id = molecule->id();
			deleted_comp_id = molecule->componentid();

			double *energy = new double;
			double *old_energy = new double;
			double force[3] = { 0, 0, 0 };
			*energy = 1;//linkedCells->getForceAndEnergy(newMolecule, force);
			*old_energy = 1;//*energy;
			for (int j = 0; j < 3; j++) {
				// invert velocity
				v_old[j] = -molecule->v(j);

				// keep angular velocity
				w_old[j] = molecule->D(j);

			}
			Quaternion q_old = molecule->q();
			double allowed_low[3] = { 0, 0, 0 };
			double allowed_high[3] = { rmax[0], rmax[1], rmax[2] };
			//molecule = moleculeContainer->deleteCurrent();
			int iterations = -1;
			int seed = 1;
			vector<double> vec_energy;
			vector<double> vec_angle;
			vector<double*> vec_lj;
			vector<double*> vec_center;
			string name_energy = "energy.txt";
			string name_angle = "angle.txt";
			string name_lj = "lj.txt";
			string name_center = "center.txt";
			std::cout << "about to call the wrapper" << std::endl;
			moleculardynamics::coupling::interface::MardynMoleculeWrapper<
					Molecule, 3> wrapper(molecule);
			moleculardynamics::coupling::ParticleInsertion<Molecule,
					LinkedCells, 3> usher;
			int boundary;

			// define allowed region to be the thin strip of width rc
			if (molecule->r(rdfBoundary) < rmin[rdfBoundary]) {
				boundary = -1;
				allowed_high[rdfBoundary] = rmin[rdfBoundary] + rc;
			}
			if (molecule->r(rdfBoundary) > rmax[rdfBoundary]) {
				boundary = 1;
				allowed_low[rdfBoundary] = rmax[rdfBoundary] - rc;
			}

			// perform usher. after usher is done molecule will have the
			// position usher found
			while (iterations == -1 && seed < 10) {
				iterations = usher.findParticlePosition(linkedCells, molecule,
						u_avg, energy, old_energy, true, false, seed, 100,
						10000, 100, 45 * 3.14 / 180, 3.14, 3.14 / 180, 0.02,
						&vec_energy, &vec_angle, &vec_lj, &vec_center,
						name_energy, name_angle, name_lj, name_center,
						allowed_low, allowed_high);
				seed++;
			}
			molecule->setv(v_old);

			molecule->setD(w_old);

			molecule = moleculeContainer->next();

		} else {
			molecule = moleculeContainer->next();
		}
		i++;
	}
	moleculeContainer->update();
}
