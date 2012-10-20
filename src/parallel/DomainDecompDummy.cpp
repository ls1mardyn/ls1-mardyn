#include <ctime>

#include "DomainDecompDummy.h"

#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "utils/Logger.h"
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

bool DomainDecompDummy::first_unif = true;
double* DomainDecompDummy::unif_rand = new double[2];
bool DomainDecompDummy::have_avg_energy = false;
int DomainDecompDummy::num_calles = 0;

DomainDecompDummy::DomainDecompDummy() {
}

DomainDecompDummy::~DomainDecompDummy() {
}

double DomainDecompDummy::getUniformRandomNumber() {
	//srand((unsigned) time(NULL));
	return ((double) rand()) / ((double) RAND_MAX);
}

double DomainDecompDummy::getGaussianRandomNumber() {
	//srand((unsigned) time(NULL));
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

double DomainDecompDummy::getAverageEnergy(LinkedCells* linkedCells,
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

void DomainDecompDummy::generateRandomVelocity(double temperature, double mass,
		double* v) {

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

void DomainDecompDummy::generateRandomAngularVelocity(double temperature,
		double* w, Domain* domain, Molecule* currentMolecule) {

	//srand((unsigned) time(NULL));
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

void DomainDecompDummy::exchangeMolecules(ParticleContainer* moleculeContainer,
		const vector<Component>& components, Domain* domain) {

	double rmin[3]; // lower corner of the process-specific domain //ENABLE_MPI
	double rmax[3];
	double halo_L[3]; // width of the halo strip //ENABLE_MPI
	for (int i = 0; i < 3; i++) {
		rmin[i] = moleculeContainer->getBoundingBoxMin(i);
		rmax[i] = moleculeContainer->getBoundingBoxMax(i);
		halo_L[i] = moleculeContainer->get_halo_L(i);
	}

	Molecule* currentMolecule;
	// molecules that have to be copied (because of halo), get a new position
	double new_position[3];

	//
	double phaseSpaceSize[3];

	double low_limit; // particles below this limit have to be copied or moved to the lower process
	double high_limit; // particles above(or equal) this limit have to be copied or moved to the higher process

	for (unsigned short d = 0; d < 3; ++d) {
		phaseSpaceSize[d] = rmax[d] - rmin[d];
		// set limits (outside "inner" region)
		low_limit = rmin[d] + halo_L[d];
		high_limit = rmax[d] - halo_L[d];
		currentMolecule = moleculeContainer->begin();

		//cout << "low_limit: " << low_limit << " / high_limit: " << high_limit << endl;
		//cout << "halo_L: " << halo_L[0] << " / " << halo_L[1] << " / " << halo_L[2] << endl;
		//cout << "proc_domain_L: " << proc_domain_L[0] << " / " << proc_domain_L[1] << " / " << proc_domain_L[2] << endl;
		while (currentMolecule != moleculeContainer->end()) {
			const double& rd = currentMolecule->r(d);
			if (rd < low_limit) {
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
				//				if (currentMolecule->r(0) < moleculeContainer->get_halo_L(0)
				//						&& currentMolecule->r(1)
				//								< moleculeContainer->get_halo_L(1)
				//						&& currentMolecule->r(2)
				//								< moleculeContainer->get_halo_L(2)) {
				//					std::cout << "Periodic copy1: " << m1.r(0) << " "
				//							<< m1.r(1) << " " << m1.r(2)
				//							<< " of molecule with " << currentMolecule->id()
				//							<< " and pos " << currentMolecule->r(0) << " "
				//							<< currentMolecule->r(1) << " "
				//							<< currentMolecule->r(2) << std::endl;
				//				}
				moleculeContainer->addParticle(m1);
				currentMolecule = moleculeContainer->next();
			} else if (rd >= high_limit) {
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
				//				if (currentMolecule->r(0) < moleculeContainer->get_halo_L(0)
				//						&& currentMolecule->r(1)
				//								< moleculeContainer->get_halo_L(1)
				//						&& currentMolecule->r(2)
				//								< moleculeContainer->get_halo_L(2)) {
				//					std::cout << "Periodic copy2: " << m1.r(0) << " "
				//							<< m1.r(1) << " " << m1.r(2)
				//							<< " of molecule with " << currentMolecule->id()
				//							<< " and pos " << currentMolecule->r(0) << " "
				//							<< currentMolecule->r(1) << " "
				//							<< currentMolecule->r(2) << std::endl;
				//				}
				moleculeContainer->addParticle(m1);
				currentMolecule = moleculeContainer->next();
			} else
				currentMolecule = moleculeContainer->next();
		}
	}
	//moleculeContainer->updateMoleculeCaches();
	//	if (num_calles > 1)
	//		validateUsher(moleculeContainer, components, domain);
	num_calles++;
}

void DomainDecompDummy::validateUsher(ParticleContainer* moleculeContainer,
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

	// choose a random molecule to delete and delete it
	double rnd = (double) rand() / (double) RAND_MAX;
	int int_rnd = 1 + 9826 * rnd;
	std::cout << int_rnd << std::endl;
	Molecule* molecule;
	int deleted_id, deleted_comp_id, i = 1;
	double v_old[3];
	double w_old[3];
	double r_old[3];
	molecule = moleculeContainer->begin();
	double rc = moleculeContainer->getCutoff();
	bool found = false;
	LinkedCells* linkedCells = (LinkedCells*) moleculeContainer;
	Molecule* m;

	if (num_calles > 1) {
		u_avg = this->getAverageEnergy(linkedCells, rmin, rmax);//domain->getAverageGlobalUpot();
		std::cout << "uavg1 " << u_avg << std::endl;
		//std::cout<<"uavg2 "<<this->getAverageEnergy(linkedCells, rmin, rmax);
		std::cout << "average energy: " << u_avg << std::endl;
//		if (!have_avg_energy) {
//			energy_file
//					= fopen(
//							"/home_local/kovacevt/Desktop/thesis_rep/masters-thesis-kovacevic-tijana/Ethan_10k_supercritical/results/Ethan_10k_supercritica_validate_usher_energy.txt",
//							"w");
//		}
		have_avg_energy = true;
	}

	double temperature = domain->getCurrentTemperature(0);
	double v[3] = { 0, 0, 0 };
	double w[3] = { 0, 0, 0 };
	this->generateRandomVelocity(temperature, molecule->mass(), v);
	this->generateRandomAngularVelocity(temperature, w, domain, molecule);

	while (molecule != moleculeContainer->end()) {
		if (molecule->id() == 2794) {
			m = molecule;
		}
		if (molecule->id() == int_rnd && molecule->r(0) > rmin[0]
				&& molecule->r(0) < rmax[0] && molecule->r(1) > rmin[1]
				&& molecule->r(1) < rmax[1] && molecule->r(2)
				> rmin[2] && molecule->r(2) < rmax[2]) {
			deleted_id = molecule->id();
			deleted_comp_id = molecule->componentid();
//			fprintf(energy_file, "%g \n", domain->getAverageGlobalUpot());
//			fflush(energy_file);
			std::cout << "ushering molecule " << molecule->id()<<std::endl;
			double *energy = new double;
			double *old_energy = new double;
			double force[3] = { 0, 0, 0 };
			*energy = 1;//linkedCells->getForceAndEnergy(newMolecule, force);
			*old_energy = 1;//*energy;
			for (int j = 0; j < 3; j++) {
				v_old[j] = molecule->v(j);
				w_old[j] = molecule->D(j);
				r_old[j] = molecule->r(j);
			}
			double allowed_low[3] = { 0, 0, 0 };
			std::cout << "halol: " << halo_L[0] << std::endl;
			double allowed_high[3] =
					{ rmax[0], rmax[1], rmax[2]};
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
			while (iterations == -1 && seed < 10) {
				std::cout<<"starting usher "<<std::endl;
				iterations = usher.findParticlePosition(linkedCells,
						molecule, u_avg, energy, old_energy, true, false,
						seed, 100, 10000, 100, 45 * 3.14 / 180, 3.14,
						3.14 / 180, 0.02, &vec_energy, &vec_angle, &vec_lj,
						&vec_center, name_energy, name_angle, name_lj,
						name_center, allowed_low, allowed_high);
				seed++;
			}
			molecule->setv(v);
			molecule->setD(w);
			found = true;
			//break;
			molecule = moleculeContainer->next();

		} else {
			molecule = moleculeContainer->next();
		}
		i++;
	}
	moleculeContainer->update();


}

void DomainDecompDummy::balanceAndExchange(bool balance,
		ParticleContainer* moleculeContainer,
		const vector<Component>& components, Domain* domain) {
	exchangeMolecules(moleculeContainer, components, domain);
}

unsigned long DomainDecompDummy::countMolecules(
		ParticleContainer* moleculeContainer, vector<unsigned long> &compCount) {
	for (unsigned i = 0; i < compCount.size(); i++) {
		compCount[i] = 0;
	}
	Molecule* tempMolecule;
	for (tempMolecule = moleculeContainer->begin(); tempMolecule
			!= moleculeContainer->end(); tempMolecule
			= moleculeContainer->next()) {
		compCount[tempMolecule->componentid()] += 1;
	}
	int numMolecules = 0;
	for (unsigned i = 0; i < compCount.size(); i++) {
		numMolecules += compCount[i];
	}
	return numMolecules;
}

double DomainDecompDummy::getBoundingBoxMin(int dimension, Domain* domain) {
	return 0.0;
}
double DomainDecompDummy::getBoundingBoxMax(int dimension, Domain* domain) {
	return domain->getGlobalLength(dimension);
}

double DomainDecompDummy::getTime() {
	return double(clock()) / CLOCKS_PER_SEC;
}

unsigned DomainDecompDummy::Ndistribution(unsigned localN, float* minrnd,
		float* maxrnd) {
	*minrnd = 0.0;
	*maxrnd = 1.0;
	return localN;
}

void DomainDecompDummy::assertIntIdentity(int IX) {
}

void DomainDecompDummy::assertDisjunctivity(TMoleculeContainer* mm) {
}

void DomainDecompDummy::printDecomp(std::string filename, Domain* domain) {
	global_log->warning() << "printDecomp useless in serial mode" << endl;
}
