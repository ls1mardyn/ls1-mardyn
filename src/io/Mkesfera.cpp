#include "Mkesfera.h"

/* the following macro has to be defined to use math constants in cmath */
#define _USE_MATH_DEFINES  1

#include <cmath>
#include <vector>

#include "Domain.h"
#include "ensemble/EnsembleBase.h"
#include "molecules/Component.h"
#include "particleContainer/ParticleContainer.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "utils/Random.h"

#define VARFRACTION 0.07
#define BOXOVERLOAD 1.3333

//might be necessary to increase for exascale
#define OVERLAPFACTOR 1.5


void MkesferaGenerator::readXML(XMLfileUnits& xmlconfig) {
#define MAX(a, b) (((a) >= (b)) (a) : (b));
	double boxLength = 0;
	for(int d = 0; d < 3; d++) {
		double length = _simulation.getEnsemble()->domain()->length(d);
		if(boxLength < length) {
			boxLength = length;
		}
	}
	Log::global_log->info() << "Box length: " << boxLength << std::endl;
	R_o = boxLength / 2;

	xmlconfig.getNodeValueReduced("outer-density", rho_o);
	Log::global_log->info() << "Outer density: " << rho_o << std::endl;

	xmlconfig.getNodeValueReduced("droplet/radius", R_i);
	Log::global_log->info() << "Droplet radius: " << R_i << std::endl;
	xmlconfig.getNodeValueReduced("droplet/density", rho_i);
	Log::global_log->info() << "Droplet density: " << rho_i << std::endl;
	for(int d = 0; d < 3; d++) {
		center[d] = R_o;
	}
	xmlconfig.getNodeValueReduced("droplet/center/x", center[0]);
	xmlconfig.getNodeValueReduced("droplet/center/y", center[1]);
	xmlconfig.getNodeValueReduced("droplet/center/z", center[2]);
	Log::global_log->info() << "Droplet center: " << center[0] << ", " << center[0] << ", " << center[0] << std::endl;
}

unsigned long
MkesferaGenerator::readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase*) {

	int fl_units;
	double rhomax = (rho_i > rho_o) ? rho_i : rho_o;
	double N_boxes = 8.0 * R_o * R_o * R_o * (BOXOVERLOAD * rhomax) / 3.0;
	fl_units = ceil(pow(N_boxes, 1.0 / 3.0));
	double fl_unit = 2.0 * R_o / (double) fl_units;


	int fl_units_local[3];
	double box_max[3];
	//borders of local subregion in parallel computation
	double box_min_local[3];
	double box_max_local[3];
	/* box min is assumed to be 0 */
	int startx[3];
	int endx[3];
	for(int d = 0; d < 3; d++) {
		box_max[d] = _simulation.getEnsemble()->domain()->length(d);
#ifdef ENABLE_MPI
		box_max_local[d] = particleContainer->getBoundingBoxMax(d);
		box_min_local[d] = particleContainer->getBoundingBoxMin(d);

#else
		box_min_local[d] = 0;
		fl_units_local[d] = fl_units;
		box_max_local[d] = box_max[d];
#endif
	}
	for(int d = 0; d < 3; d++) {
		startx[d] = floor(box_min_local[d] / fl_unit - 0.5);
		endx[d] = ceil(box_max_local[d] / fl_unit);

	}
	for(int d = 0; d < 3; d++) {

		endx[d] = std::min(fl_units - 1, endx[d] + 1);

		startx[d] = std::max(0, startx[d] - 1);


		fl_units_local[d] = endx[d] - startx[d] + 1;

	}

	double T = _simulation.getEnsemble()->T();
	Log::global_log->info() << "Temperature: " << T << std::endl;

	double cutoff = _simulation.getcutoffRadius();
	Random* rnd = new Random();
	rnd->init(
			(int) (10000.0 * R_o) - (int) (3162.3 * cutoff)
			+ (int) (1000.0 * T) - (int) (316.23 * rho_i)
			+ (int) (100.0 * R_i)
	);

	bool**** fill = new bool*** [fl_units_local[0]];
	for(int i = 0; i < fl_units_local[0]; i++) {
		fill[i] = new bool** [fl_units_local[1]];
		for(int j = 0; j < fl_units_local[1]; j++) {
			fill[i][j] = new bool* [fl_units_local[2]];
			for(int k = 0; k < fl_units_local[2]; k++) {
				fill[i][j][k] = new bool[3];
			}
		}
	}
	unsigned slots = 3.0 * fl_units * fl_units * fl_units;
	double boxdensity = (double) slots / (8.0 * R_o * R_o * R_o);
	Log::global_log->debug() << "Box density: " << boxdensity << " (unit cell: " << fl_unit << ")" << std::endl;
	double P_in = rho_i / boxdensity;
	double P_out = rho_o / boxdensity;
	Log::global_log->debug() << "Insertion probability: " << P_in << " inside, " << P_out << " outside" << std::endl;

	/* box min is assumed to be 0 (not in parallel!)*/
	for(int d = 0; d < 3; d++) {
		box_max[d] = _simulation.getEnsemble()->domain()->length(d);
	}

	double goffset[3][3]; /**< relative koordinates of face centers */
	goffset[0][0] = 0.0;
	goffset[1][0] = 0.5;
	goffset[2][0] = 0.5;
	goffset[0][1] = 0.5;
	goffset[1][1] = 0.0;
	goffset[2][1] = 0.5;
	goffset[0][2] = 0.5;
	goffset[1][2] = 0.5;
	goffset[2][2] = 0.0;

	unsigned N = 0;
	int idx[3];


	for(idx[0] = 0; idx[0] < fl_units; idx[0]++) {
		for(idx[1] = 0; idx[1] < fl_units; idx[1]++) {
			for(idx[2] = 0; idx[2] < fl_units; idx[2]++) {
				for(int p = 0; p < 3; p++) {
					double qq = 0.0;
					double q[3];
					bool notInBox = 0;
					for(int d = 0; d < 3; d++) {
						q[d] = (idx[d] + goffset[d][p]) * fl_unit;
						if(q[d] > box_max_local[d] or q[d] < box_min_local[d]) {
							notInBox = 1;

						}
						q[d] -= center[d];
						q[d] = q[d] - round(q[d] / box_max[d]) * box_max[d];
						qq += q[d] * q[d];
					}
					double tP = (qq > R_i * R_i) ? P_out : P_in;
					bool tfill = (tP >= rnd->rnd());
					if(notInBox) {
						if(idx[0] >= startx[0] and idx[0] <= endx[0] and idx[1] >= startx[1] and idx[1] <= endx[1] and
						   idx[2] >= startx[2] and idx[2] <= endx[2]) {
							fill[idx[0] - startx[0]][idx[1] - startx[1]][idx[2] - startx[2]][p] = 0;
						}
						continue;
					}


					if(idx[0] - startx[0] >= fl_units_local[0] or idx[1] - startx[1] >= fl_units_local[1] or
					   idx[2] - startx[2] >= fl_units_local[2] or startx[0] > idx[0] or startx[1] > idx[1] or
					   startx[2] > idx[2]) {
						Log::global_log->error() << "Error in calculation of start and end values! \n";
						Simulation::exit(0);
					}
					fill[idx[0] - startx[0]][idx[1] - startx[1]][idx[2] - startx[2]][p] = tfill;
					if(tfill) {
						N++;
					}
				}
			}
		}
	}

	int startID = 0;
#ifdef ENABLE_MPI
	MPI_Exscan(&N, &startID, 1 , MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	MPI_Allreduce(MPI_IN_PLACE, &N, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
	Log::global_log->debug() << "Filling " << N << " out of " << slots << " slots" << std::endl;
	Log::global_log->debug() << "Density: " << N / (8.0 * R_o * R_o * R_o) << std::endl;


	double v_avg = sqrt(3.0 * T);

	Component* component = _simulation.getEnsemble()->getComponent(0);
	unsigned ID = 1 + startID;
	unsigned int numberOfMolecules = 0;
	for(idx[0] = 0; idx[0] < fl_units; idx[0]++) {
		for(idx[1] = 0; idx[1] < fl_units; idx[1]++) {
			for(idx[2] = 0; idx[2] < fl_units; idx[2]++) {
				for(unsigned p = 0; p < 3; p++) {
					float random1 = rnd->rnd();
					float random2 = rnd->rnd();
					float random3 = rnd->rnd();//every MPI process needs to get the same random numbers for the same molecules
					if(idx[0] >= startx[0] and idx[0] <= endx[0] and idx[1] >= startx[1] and idx[1] <= endx[1] and
					   idx[2] >= startx[2] and idx[2] <= endx[2]) {
						if(fill[idx[0] - startx[0]][idx[1] - startx[1]][idx[2] - startx[2]][p]) {
							//global_log->debug() << "Inserting: " << idx[0] << "," << idx[1] << "," << idx[2] << "; " << p << std::endl;
							double q[3];
							bool notInBox = false;
							for(int d = 0; d < 3; d++) {
								q[d] = (idx[d] + VARFRACTION * (random1 - 0.5) + goffset[d][p]) * fl_unit;
								if(q[d] < box_min_local[d]) notInBox = true;
								else if(q[d] > box_max_local[d]) notInBox = true;
							}
							if(notInBox) continue;
							double phi = 2 * M_PI * random2;
							double omega = 2 * M_PI * random3;

							double v[3];
							v[0] = v_avg * cos(phi) * cos(omega);
							v[1] = v_avg * cos(phi) * sin(omega);
							v[2] = v_avg * sin(phi);
							Molecule molecule(ID, component, q[0], q[1], q[2], v[0], v[1], v[2], 1, 0, 0, 0, 0, 0, 0);

							// particle position is already checked for this container.
							particleContainer->addParticle(molecule, true, false);

							ID++;
							numberOfMolecules++;
						}
					}
				}
			}
		}
	}

	for(int i = 0; i < fl_units_local[0]; i++) {
		for(int j = 0; j < fl_units_local[1]; j++) {
			for(int k = 0; k < fl_units_local[2]; k++) {
				delete[] fill[i][j][k];
			}
			delete[] fill[i][j];
		}
		delete[] fill[i];
	}
	delete[] fill;
#ifdef ENABLE_MPI
	MPI_Allreduce(MPI_IN_PLACE, &numberOfMolecules, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
	domain->setGlobalTemperature(T);
	domain->setglobalNumMolecules(numberOfMolecules);
	domain->setglobalRho(numberOfMolecules / _simulation.getEnsemble()->V());

	Log::global_log->info() << "Inserted number of molecules: " << numberOfMolecules << std::endl;
	return ID;
}

