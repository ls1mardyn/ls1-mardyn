/* the following macro has to be defined to use math constants in cmath */
#define _USE_MATH_DEFINES  1

#include "TcTS.h"

#include <cmath>

#include "Domain.h"
#include "ensemble/EnsembleBase.h"
#include "molecules/Molecule.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "utils/Random.h"


#define PRECISION 5
#define VARFRACTION 0.125

void MkTcTSGenerator::readXML(XMLfileUnits& xmlconfig) {
	//this is the ratio of the length of the domain box that layer1 encompasses (has to be between 0 and 1)
	l1_ratio = 0.5;
	//this is the offset of the beginning of layer 1 from the box_boundaries given in realtive values to the box_length
	//if l1_offset is 0.5 layer1 begins in the middle of the box
	//has to be between 0 and 1
	l1_offset = 0.3;
	xmlconfig.getNodeValue("layer1/l1_ratio", l1_ratio);
	Log::global_log->info() << "Layer 1, ratio: " << l1_ratio << std::endl;
	xmlconfig.getNodeValue("layer1/l1_offset", l1_offset);
	xmlconfig.getNodeValue("layer1/density", rho1);
	Log::global_log->info() << "Layer 1, density: " << rho1 << std::endl;
	rho2 = rho1;
	xmlconfig.getNodeValue("layer2/density", rho2);
	Log::global_log->info() << "Layer 2, density: " << rho2 << std::endl;
}

void MkTcTSGenerator::readPhaseSpaceHeader(Domain* domain, double /*timestep*/) {
}

unsigned long
MkTcTSGenerator::readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase*) {
	// Mixing coefficients
	std::vector<double>& dmixcoeff = domain->getmixcoeff();
	dmixcoeff.clear();

	unsigned int numcomponents = _simulation.getEnsemble()->getComponents()->size();
	for(unsigned int i = 1; i < numcomponents; i++) {
		for(unsigned int j = i + 1; j <= numcomponents; j++) {
			double xi = 1., eta = 1.;

			dmixcoeff.push_back(xi);
			dmixcoeff.push_back(eta);
		}
	}

	double T = _simulation.getEnsemble()->T();
	double box[3];
	for(int d = 0; d < 3; d++) {
		box[d] = _simulation.getEnsemble()->domain()->length(d);
	}

	size_t fl_units[3][2];
	double fl_unit[3][2];
	double N_id[2];
	for(int i = 0; i < 2; i++) {
		double curr_ratio = (i == 0) ? l1_ratio : (1 - l1_ratio);
		N_id[i] = box[0] * (curr_ratio * box[1]) * box[2] * ((i == 0) ? rho1 : rho2);
		double N_boxes = N_id[i] / 3.0; /* 3 molecules per box */
		fl_units[1][i] = static_cast<size_t>(round(
				pow(
						(N_boxes * (curr_ratio * box[1]) * (curr_ratio * box[1]))
						/ (box[0] * box[2]), 1.0 / 3.0
				)
		));
		if(fl_units[1][i] == 0ul) fl_units[1][i] = 1;
		double bxbz_id = N_boxes / fl_units[1][i];
		fl_units[0][i] = static_cast<size_t>(round(sqrt(box[0] * bxbz_id / box[2])));
		if(fl_units[0][i] == 0ul) fl_units[0][i] = 1;
		fl_units[2][i] = static_cast<size_t>(ceil(bxbz_id / fl_units[0][i]));
		for(int d = 0; d < 3; d++) fl_unit[d][i] = ((d == 1) ? curr_ratio : 1.0) * box[d] / static_cast<double>(fl_units[d][i]);
		Log::global_log->debug() << "Elementary cell " << i << ": " << fl_unit[0][i] << " x " << fl_unit[1][i] << " x "
							<< fl_unit[2][i] << std::endl;
	}

	Random* rnd = new Random();
	rnd->init(
			(int) (10000.0 * box[0]) - (int) (3162.3 * _simulation.getcutoffRadius())
			+ (int) (1000.0 * T)
			+ (int) (100.0 * box[1])
	);

	//bool fill0[fl_units[0][0]][fl_units[1][0]][fl_units[2][0]][3];
	std::vector<bool> fill0(fl_units[0][0] * fl_units[1][0] * fl_units[2][0] * 3);
	//bool fill1[fl_units[0][1]][fl_units[1][1]][fl_units[2][1]][3];
	std::vector<bool> fill1(fl_units[0][1] * fl_units[1][1] * fl_units[2][1] * 3);
	unsigned N[2];
	unsigned slots[2];
	for(unsigned l = 0; l < 2; l++) {
		for(unsigned i = 0; i < fl_units[0][l]; i++) {
			for(unsigned j = 0; j < fl_units[1][l]; j++) {
				for(unsigned k = 0; k < fl_units[2][l]; k++) {
					for(unsigned d = 0; d < 3; d++) {
						if(l == 0)
							fill0[((fl_units[1][0] * i + j) * fl_units[2][0] + k) * 3 + d] = true;
						else
							fill1[((fl_units[1][1] * i + j) * fl_units[2][1] + k) * 3 + d] = true;
					}
				}
			}
		}
		slots[l] = 3 * fl_units[0][l] * fl_units[1][l] * fl_units[2][l];
		N[l] = slots[l];
	}
	bool tswap;
	double pswap;
	for(unsigned l = 0; l < 2; l++) {
		for(unsigned m = 0; m < PRECISION; m++) {
			tswap = (N[l] < N_id[l]);
			pswap = (N_id[l] - (double) N[l]) / ((tswap ? slots[l] : 0) - (double) N[l]);
			// cout << "N = " << N[l] << ", N_id = " << N_id[l] << " => tswap = " << tswap << ", pswap = " << pswap << "\n";
			for(unsigned i = 0; i < fl_units[0][l]; i++) {
				for(unsigned j = 0; j < fl_units[1][l]; j++) {
					for(unsigned k = 0; k < fl_units[2][l]; k++) {
						for(unsigned d = 0; d < 3; d++) {
							if(pswap >= rnd->rnd()) {
								if(((l == 0)
									&& fill0[((fl_units[1][0] * i + j) * fl_units[2][0] + k) * 3 + d])
								   || ((l == 1) && fill1[((fl_units[1][1] * i + j) * fl_units[2][1] + k) * 3 + d]))
									N[l]--;
								if(l == 0)
									fill0[((fl_units[1][0] * i + j) * fl_units[2][0] + k) * 3 + d] = tswap;
								else
									fill1[((fl_units[1][1] * i + j) * fl_units[2][1] + k) * 3 + d] = tswap;
								if(tswap)
									N[l]++;
							}
						}
					}
				}
			}
		}
		Log::global_log->debug() << "Filling " << N[l] << " of 3*"
							<< fl_units[0][l] << "*" << fl_units[1][l] << "*" << fl_units[2][l]
							<< " = " << slots[l] << " slots (ideally " << N_id[l] << ")" << std::endl;
	}

	double loffset[3][2];
	loffset[0][0] = 0.1;
	loffset[1][0] = l1_offset;
	loffset[2][0] = 0.1;
	double l2_offset = l1_offset + l1_ratio;
	if(l1_offset > 1) {
		--l1_offset;
	}
	loffset[0][1] = 0.1;
	loffset[1][1] = l2_offset;
	loffset[2][1] = 0.1;
	double goffset[3][3] = {
			{0.0, 0.5, 0.5},
			{0.5, 0.0, 0.5},
			{0.5, 0.5, 0.0}
	};

	double v_avg = sqrt(3.0 * T);

	std::vector<Component>* components = _simulation.getEnsemble()->getComponents();
	if(components->size() > 2) {
		Log::global_log->warning() << "MkTcTs only supports 2 components but " << components->size() << "where given!";
	}
	Component* component1 = &(*components)[0];
	Component* component2 = components->size() >= 2 ? &(*components)[1] : &(*components)[0];
	//Component* currComponent = _simulation.getEnsemble()->getComponent(0);
	unsigned ID = 1;
	for(unsigned l = 0; l < 2; l++) {
		Component* currComponent = l == 0 ? component1 : component2;
		for(unsigned i = 0; i < fl_units[0][l]; i++) {
			for(unsigned j = 0; j < fl_units[1][l]; j++) {
				for(unsigned k = 0; k < fl_units[2][l]; k++) {
					for(unsigned d = 0; d < 3; d++) {
						if(((l == 0) && fill0[((fl_units[1][0] * i + j) * fl_units[2][0] + k) * 3 + d])
						   || ((l == 1)
							   && fill1[((fl_units[1][1] * i + j) * fl_units[2][1] + k) * 3 + d])) {
							double q[3];
							q[0] = i * fl_unit[0][l];
							q[1] = j * fl_unit[1][l];
							q[2] = k * fl_unit[2][l];
							for(int m = 0; m < 3; m++) {
								q[m] += box[m] * loffset[m][l] + fl_unit[m][l] * goffset[m][d];
								q[m] += VARFRACTION * fl_unit[m][l] * (rnd->rnd() - 0.5);
								if(q[m] > box[m]) q[m] -= box[m];
							}
							double phi = 2 * M_PI * rnd->rnd();
							double omega = 2 * M_PI * rnd->rnd();

							double v[3];
							v[0] = v_avg * cos(phi) * cos(omega);
							v[1] = v_avg * cos(phi) * sin(omega);
							v[2] = v_avg * sin(phi);
							Molecule molecule(ID, currComponent, q[0], q[1], q[2], v[0], v[1], v[2], 1, 0, 0, 0, 0, 0,
											  0);
							// only add particle if it is inside of the own domain!
							if(particleContainer->isInBoundingBox(molecule.r_arr().data())) {
								particleContainer->addParticle(molecule, true, false);
							}
							ID++;
						}
					}
				}
			}
		}
	}
	domain->setGlobalTemperature(T);
	domain->setglobalNumMolecules(ID - 1);
	domain->setglobalRho(ID / _simulation.getEnsemble()->V());
	return ID;
}
