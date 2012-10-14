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
	if (particle_insertion_type == 1) std::cout<<"chosen insertion type is usehr";
	else if (particle_insertion_type == 0) std::cout<<"chosen insertion type is reflective wall";
	else std::cout<<"unknown insertion type"<<std::endl;
	// TODO Auto-generated constructor stub

}

double RDFDummyDecomposition::getUniformRandomNumber() {
	srand((unsigned) time(NULL));
	return ((double) rand()) / ((double) RAND_MAX);
}

double RDFDummyDecomposition::getGaussianRandomNumber() {
	srand((unsigned) time(NULL));
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
	for (currentMolecule = linkedCells->begin(); currentMolecule
			!= linkedCells->end(); currentMolecule = linkedCells->next()) {
		if (currentMolecule->r(0) < high_limit[0] && currentMolecule->r(0)
				> low_limit[0] && currentMolecule->r(1) < high_limit[1]
				&& currentMolecule->r(1) > low_limit[1]
				&& currentMolecule->r(2) < high_limit[2] && currentMolecule->r(
				2) > low_limit[2]) {
			u_avg += linkedCells->getForceAndEnergy(currentMolecule, force);
			num++;
		}
	}
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

	srand((unsigned) time(NULL));
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
	for (unsigned short d = 1; d < 3; ++d) {
		phaseSpaceSize[d] = rmax[d] - rmin[d];
		// set limits (outside "inner" region)
		low_limit[d] = rmin[d] + halo_L[d];
		high_limit[d] = rmax[d] - halo_L[d];
		currentMolecule = moleculeContainer->begin();

		//cout << "low_limit: " << low_limit << " / high_limit: " << high_limit << endl;
		//cout << "halo_L: " << halo_L[0] << " / " << halo_L[1] << " / " << halo_L[2] << endl;
		//cout << "proc_domain_L: " << proc_domain_L[0] << " / " << proc_domain_L[1] << " / " << proc_domain_L[2] << endl;
		while (currentMolecule != moleculeContainer->end()) {
			if (!particle_insertion_type)
				currentMolecule->enableBouncingBack(rmin[rdfBoundary], rmax[rdfBoundary], 0, 0);
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
				if (!particle_insertion_type)
					m1.enableBouncingBack(rmin[rdfBoundary], rmax[rdfBoundary], 0, 0);
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
				if (!particle_insertion_type)
					m1.enableBouncingBack(rmin[rdfBoundary], rmax[rdfBoundary], 0, 0);
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
	bool have_avg_energy = false;
	double u_avg;
	Molecule* currentMolecule = moleculeContainer->begin();
	// molecules that have to be copied (because of halo), get a new position
	double new_position[3];

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

	if (!have_avg_energy && num_calles > 1) {
		u_avg = domain->getAverageGlobalUpot();//this->getAverageEnergy(linkedCells, rmin, rmax);
		std::cout << "average energy: " << u_avg << std::endl;
		have_avg_energy = true;
	} else if (!have_avg_energy && num_calles <= 1) {
		std::cout << "calculating average energy" << std::endl;
		u_avg = this->getAverageEnergy(linkedCells, rmin, rmax);
	}
	num_calles++;
	int curr_id = 0;

	addPeriodicCopies(moleculeContainer, rmin, rmax, phaseSpaceSize, halo_L,
			components);
	std::cout << "finished inserting periodic copies" << std::endl;
	moleculeContainer->update();
	std::cout << "num particles: " << moleculeContainer->getNumberOfParticles()
			<< std::endl;
	int num = 0;

	/*
	 * If particle insertion type is usher, do usher
	 */
	currentMolecule = moleculeContainer->begin();
	while (particle_insertion_type == 1 && currentMolecule != moleculeContainer->end()) {

		num++;
		curr_id = currentMolecule->id();

		// if this is not a boundary molecule, nothing needs to be done
		if (currentMolecule->r(0) > low_limit[0] && currentMolecule->r(0)
				< high_limit[0] && currentMolecule->r(1) > low_limit[1]
				&& currentMolecule->r(1) < high_limit[1] && currentMolecule->r(
				2) > low_limit[2] && currentMolecule->r(2) < high_limit[2]) {

			currentMolecule = moleculeContainer->next();
			continue;
		}

		// if leaving the domain
		if (currentMolecule->r(rdfBoundary) < rmin[rdfBoundary]
				|| currentMolecule->r(rdfBoundary) > rmax[rdfBoundary]) {

			// save relevant values to delete and add with usher


			std::cout << "to delete: " << currentMolecule->id() << " "
					<< currentMolecule->r(0) << " " << currentMolecule->r(1)
					<< " " << currentMolecule->r(2) << std::endl;

			double* starting_position = new double[3];
			for (unsigned short d = 0; d < 3; ++d) {
				const double& rd = currentMolecule->r(d);
				if (rd < rmin[d]) {
					starting_position[d] = 2 * rmin[d] - rd; //rmax[d] + rd;//
				} else if (rd > rmax[d]) {
					starting_position[d] = 2 * rmax[d] - rd; //rd - rmax[d]; //
				} else {
					starting_position[d] = rd;
				}
			}
			std::cout << "starting position: " << starting_position[0] << " "
					<< starting_position[1] << " " << starting_position[2]
					<< std::endl;
			starting_positions.push_back(starting_position);
			//molecules_to_delete.push_back(currentMolecule);

			deleted_ids.push_back(currentMolecule->id());
			deleted_comp_ids.push_back(currentMolecule->componentid());

			//			double axis[3] = { 0, 0, 0 };
			//			if (currentMolecule->r(rdfBoundary) < rmin[0])
			//				axis[0] = 1;
			//			else
			//				axis[0] = -1;
			//			currentMolecule->bounceBack(rdfBoundary, axis);
			//			double temperature = domain->getCurrentTemperature(0);
			//			double w[3] = { 0, 0, 0 };
			//			this->generateRandomAngularVelocity(temperature, w, domain,
			//					currentMolecule);
			//			currentMolecule->setD(w);
			//			std::cout<<"w: "<<w[0]<<std::endl;
			//			// get torque
			//			currentMolecule->calcFM();
			//			double absTorque = 0;
			//			for (int i = 0; i < 3; i++)
			//				absTorque += currentMolecule->M(i) * currentMolecule-> M(i);
			//
			//			absTorque = std::sqrt(absTorque);
			//			double torque[3];
			//			for (int i = 0; i < 3; i++)
			//				torque[i] = (1 / absTorque) * currentMolecule->M(i);
			//
			//
			//			double energy = 500;
			//			int i = 1;
			//			Quaternion q = currentMolecule->q();
			//			Quaternion old_q = q;
			//			srand ( (unsigned)time(NULL));
			//			while (energy >= 500 && i < 1000) {
			//				double angle = -30 * 3.14 / 180 + ((double) rand() / (double) RAND_MAX) * 60 * 3.14
			//									/ 180;
			//				std::cout<<angle<<std::endl;
			////				q.multiply_left(Quaternion(cos(angle / 2), torque[0] * sin(angle / 2),
			////				torque[1] * sin(angle / 2), torque[2] * sin(angle / 2)));
			//				q = Quaternion(rand(), rand(), rand(), rand());
			//				q.normalize();
			//				currentMolecule->setq(q);
			//				currentMolecule->upd_cache();
			//				std::cout<<"try number: "<<i<<" energy: "<<energy<<std::endl;
			//				i++;
			//				double tempf[3] = {0, 0, 0};
			//				energy = linkedCells->getForceAndEnergy(currentMolecule, tempf);
			//				q = old_q;
			//			}
			currentMolecule = moleculeContainer->deleteCurrent();//moleculeContainer->next();
			currentMolecule = moleculeContainer->next();
		} else {
			currentMolecule = moleculeContainer->next();
		}

	}
	//	std::cout << "num " << num << std::endl;
	std::cout << "finished searching for usher candidates" << std::endl;
	//
	//	std::vector<double*> velocities;
	//	// delete the ones that need to be deleted
	//
	linkedCells->update();
	std::cout << "updated" << std::endl;

	// insert with usher

	for (unsigned int i = 0; i < starting_positions.size(); i++) {

		double cell_coordinates[3];
		std::cout << "processing " << i << " out of "
				<< starting_positions.size() << std::endl;
		std::cout << starting_positions[i][0] << " "
				<< starting_positions[i][1] << " " << starting_positions[i][2]
				<< std::endl;

		linkedCells->getCellCoordinates(starting_positions[i], cell_coordinates);

		std::cout << "coordinates " << cell_coordinates[0] << " "
				<< cell_coordinates[1] << " " << cell_coordinates[2]
				<< std::endl;
		double* density_bins = simulation->getDensityBins();
		//double* density_bins = new double[200];
		//		for (int bin = 0; bin < 200; bin++) {
		//			density_bins[bin] = sim_density_bins[2 * bin] + sim_density_bins[2 * bin + 1];
		//		}
		int min_bin_idx = 0;
		int min_bin = 10000;
		double rnd = ((double) rand()) / ((double) RAND_MAX);
		if (rnd > 0.5) {
			for (int bin = 0; bin < 400; bin++) {
				if (density_bins[bin] < min_bin) {
					min_bin = density_bins[bin];
					min_bin_idx = bin;
				}
			}
		} else {
			for (int bin = 399; bin >= 0; bin--) {
				if (density_bins[bin] < min_bin) {
					min_bin = density_bins[bin];
					min_bin_idx = bin;
				}
			}
		}
		double length = rmax[0] / 400;
		std::cout << "bin " << min_bin_idx << " with " << min_bin << std::endl;
		//		double allowed_low[3] = {std::max(
		//										starting_positions[i][0] - 0.2 * linkedCells->cellLength()[0], rmin[0]),
		//				 0, 0 };
		//double allowed_low[3] = {std::max(cell_coordinates[0] - 1 * linkedCells->cellLength()[0], rmin[0]), std::max(cell_coordinates[1] - 1 * linkedCells->cellLength()[1], rmin[1]), std::max(cell_coordinates[2] - 1 * linkedCells->cellLength()[2], rmin[2])};
		//		double allowed_high[3] = {
		//			std::min(
		//					cell_coordinates[0] + 0 * linkedCells->cellLength()[0], rmax[0]), std::min(
		//							cell_coordinates[1] + 0 * linkedCells->cellLength()[1],
		//							rmax[1]), std::min(
		//							cell_coordinates[2] + 0 * linkedCells->cellLength()[2],
		//							rmax[2])};
		//		double allowed_high[3] = {std::min(
		//										starting_positions[i][0] + 0.2 * linkedCells->cellLength()[0], rmax[0]), rmax[1], rmax[2] };

		double allowed_low[3] = { 0, 0, 0 };
		double allowed_high[3] = { rmax[0], rmax[1], rmax[2] };
		std::cout << "allowed low and high: " << allowed_low[0] << " "
				<< allowed_low[1] << " " << allowed_low[2] << " high: "
				<< allowed_high[0] << " " << allowed_high[1] << " "
				<< allowed_high[2] << std::endl;
		Molecule* newMolecule = new Molecule(deleted_ids[i],
				deleted_comp_ids[i], starting_positions[i][0],
				starting_positions[i][1], starting_positions[i][2], 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, &components);

		// needed for usher
		double *energy = new double;
		double *old_energy = new double;
		double force[3] = { 0, 0, 0 };
		*energy = 1;//linkedCells->getForceAndEnergy(newMolecule, force);
		*old_energy = 1;//*energy;

		// stuff for debugging
		vector<double> vec_energy;
		vector<double> vec_angle;
		vector<double*> vec_lj;
		vector<double*> vec_center;
		string name_energy = "energy.txt";
		string name_angle = "angle.txt";
		string name_lj = "lj.txt";
		string name_center = "center.txt";
		std::cout << "about to call the wrapper" << std::endl;
		moleculardynamics::coupling::interface::MardynMoleculeWrapper<Molecule,
				3> wrapper(newMolecule);
		moleculardynamics::coupling::ParticleInsertion<Molecule, LinkedCells, 3>
				usher;
		std::cout << "started usher, target energy: " << u_avg << std::endl;
		int iterations = -1;
		int seed = 1;
		while (iterations == -1 && seed < 10) {
			iterations = usher.findParticlePosition(linkedCells, newMolecule,
					u_avg, energy, old_energy, true, false, seed, 100, 10000,
					100, 45 * 3.14 / 180, 3.14, 3.14 / 180, 0.01, &vec_energy,
					&vec_angle, &vec_lj, &vec_center, name_energy, name_angle,
					name_lj, name_center, allowed_low, allowed_high);
			seed++;
		}
		if (seed == 10) {
			allowed_low[0] = 0;
			allowed_high[0] = rmax[0];
			iterations = usher.findParticlePosition(linkedCells, newMolecule,
					u_avg, energy, old_energy, true, false, seed, 100, 10000,
					100, 45 * 3.14 / 180, 3.14, 3.14 / 180, 0.01, &vec_energy,
					&vec_angle, &vec_lj, &vec_center, name_energy, name_angle,
					name_lj, name_center, allowed_low, allowed_high);
		}

		std::cout << "usher did " << iterations << " iterations " << std::endl;
		std::cout << std::endl;

		// random velocity
		double temperature = domain->getCurrentTemperature(0);
		double v[3] = { 0, 0, 0 };
		double w[3] = { 0, 0, 0 };
		this->generateRandomVelocity(temperature, newMolecule->mass(), v);
		this->generateRandomAngularVelocity(temperature, w, domain, newMolecule);
		std::cout << "velocity before " << newMolecule->v(0) << " "
				<< newMolecule->v(1) << " " << newMolecule->v(2) << std::endl;
		newMolecule->setv(v);
		std::cout << "velocity after " << newMolecule->v(0) << " "
				<< newMolecule->v(1) << " " << newMolecule->v(2) << std::endl;
		std::cout << "angular velocity before " << newMolecule->D(0) << " "
				<< newMolecule->D(1) << " " << newMolecule->D(2) << std::endl;
		newMolecule->setD(w);
		std::cout << "angular velocity after " << newMolecule->D(0) << " "
				<< newMolecule->D(1) << " " << newMolecule->D(2) << std::endl;

		std::cout << "new position " << newMolecule->r(0) << " "
				<< newMolecule->r(1) << " " << newMolecule->r(2) << std::endl;

		if (iterations != -1) {
			// if usher succeeded
			moleculeContainer->addParticle(*newMolecule);
			moleculeContainer->update();
			std::cout << "added " << std::endl;

			// add periodic copies of the usher's particle, periodic in one coord
			double new_position_diag_periodic[3] = { newMolecule->r(0),
					newMolecule->r(1), newMolecule->r(2) };
			int diag = 0;
			for (int d = 0; d < 3; d++) {
				double new_position_periodic[3] = { newMolecule->r(0),
						newMolecule->r(1), newMolecule->r(2) };
				bool add_periodic = false;
				if (d != rdfBoundary && newMolecule->r(d) < low_limit[d]) {
					new_position_periodic[d] += phaseSpaceSize[d];
					new_position_diag_periodic[d] += phaseSpaceSize[d];
					add_periodic = true;
					diag++;
				} else if (d != rdfBoundary && newMolecule->r(d)
						> high_limit[d]) {
					new_position_periodic[d] -= phaseSpaceSize[d];
					new_position_diag_periodic[d] -= phaseSpaceSize[d];
					add_periodic = true;
					diag++;
				}
				if (add_periodic) {
					Molecule m1 = Molecule(newMolecule->id(),
							newMolecule->componentid(),
							new_position_periodic[0], new_position_periodic[1],
							new_position_periodic[2], newMolecule->v(0),
							newMolecule->v(1), newMolecule->v(2),
							newMolecule->q().qw(), newMolecule->q().qx(),
							newMolecule->q().qy(), newMolecule->q().qz(),
							newMolecule->D(0), newMolecule->D(1),
							newMolecule->D(2), &components);
					moleculeContainer->addParticle(m1);
					moleculeContainer->update();
					std::cout
							<< "ATTENTION: added the periodic copy of what usher added "
							<< m1.r(0) << " " << m1.r(1) << " " << m1.r(2)
							<< std::endl;
				}

			}

			// diagonal periodic copy if usher added in a corner
			if (diag == 2) {
				Molecule m1 = Molecule(newMolecule->id(),
						newMolecule->componentid(),
						new_position_diag_periodic[0],
						new_position_diag_periodic[1],
						new_position_diag_periodic[2], newMolecule->v(0),
						newMolecule->v(1), newMolecule->v(2),
						newMolecule->q().qw(), newMolecule->q().qx(),
						newMolecule->q().qy(), newMolecule->q().qz(),
						newMolecule->D(0), newMolecule->D(1),
						newMolecule->D(2), &components);
				moleculeContainer->addParticle(m1);
				moleculeContainer->update();
				std::cout
						<< "ATTENTION: added the periodic copy of what usher added diagonally"
						<< m1.r(0) << " " << m1.r(1) << " " << m1.r(2)
						<< std::endl;
			}

			// now add the periodic copy diagonally if necessary

		} else {
			std::cout << "IMPORTANT!!! usher failed, not inserting"
					<< std::endl;
		}

	}

}
