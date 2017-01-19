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

using namespace std;
using Log::global_log;

#define PRECISION 5
#define VARFRACTION 0.125

void MkTcTSGenerator::readXML(XMLfileUnits& xmlconfig) {
// 	TODO: Add option to controll the layer thicknesses
// 	xmlconfig.getNodeValue("layer1/heigth", heigth1);
// 	global_log->info() << "Layer 1, heigth: " << heigth1 << endl;
	xmlconfig.getNodeValue("layer1/density", rho1);
	global_log->info() << "Layer 1, density: " << rho1 << endl;
	rho2 = rho1;
	xmlconfig.getNodeValue("layer2/density", rho2);
	global_log->info() << "Layer 2, density: " << rho2 << endl;
}

long unsigned int MkTcTSGenerator::readPhaseSpace(ParticleContainer* particleContainer, list< ChemicalPotential >* lmu, Domain* domain, DomainDecompBase* domainDecomp) {

	double T = _simulation.getEnsemble()->T();
	double box[3];
	for( int d = 0; d < 3; d++ ) {
	   box[d] = _simulation.getEnsemble()->domain()->length(d);
	}

	unsigned fl_units[3][2];
	double fl_unit[3][2];
	double N_id[2];
	for(int i=0; i < 2; i++) {
		N_id[i] = box[0]*(0.5*box[1])*box[2] * ((i == 0)? rho1: rho2);
		double N_boxes = N_id[i] / 3.0; /* 3 molecules per box */
		fl_units[1][i] = (unsigned int) round(
							pow(
								(N_boxes * (0.5*box[1]) * (0.5*box[1]))
										/ (box[0] * box[2]), 1.0/3.0
							)
						);
		if(fl_units[1][i] == 0) fl_units[1][i] = 1;
		double bxbz_id = N_boxes / fl_units[1][i];
		fl_units[0][i] = (unsigned int) round(sqrt(box[0] * bxbz_id / box[2]));
		if(fl_units[0][i] == 0) fl_units[0][i] = 1;
		fl_units[2][i] = (unsigned int) ceil(bxbz_id / fl_units[0][i]);
		for(int d=0; d < 3; d++) fl_unit[d][i] = ((d == 1)? 0.5: 1.0) * box[d] / (double)fl_units[d][i];
		global_log->debug() << "Elementary cell " << i << ": " << fl_unit[0][i] << " x " << fl_unit[1][i] << " x " << fl_unit[2][i] << endl;
	}

	Random* rnd = new Random();
	rnd->init(
		(int)(10000.0*box[0]) - (int)(3162.3*_simulation.getcutoffRadius())
			+ (int)(1000.0*T)
			+ (int)(100.0*box[1])
	);

	bool fill0[fl_units[0][0]][fl_units[1][0]][fl_units[2][0]][3];
	bool fill1[fl_units[0][1]][fl_units[1][1]][fl_units[2][1]][3];
	unsigned N[2];
	unsigned slots[2];
	for(unsigned l=0; l < 2; l++) {
		for(unsigned i=0; i < fl_units[0][l]; i++) {
			for(unsigned j=0; j < fl_units[1][l]; j++) {
				for(unsigned k=0; k < fl_units[2][l]; k++) {
					for(unsigned d=0; d < 3; d++) {
						if(l == 0) fill0[i][j][k][d] = true;
						else fill1[i][j][k][d] = true;
					}
				}
			}
		}
		slots[l] = 3 * fl_units[0][l] * fl_units[1][l] * fl_units[2][l];
		N[l] = slots[l];
	}
	bool tswap;
	double pswap;
	for(unsigned l=0; l < 2; l++) {
		for(unsigned m=0; m < PRECISION; m++) {
			tswap = (N[l] < N_id[l]);
			pswap = (N_id[l] - (double)N[l]) / ((tswap? slots[l]: 0) - (double)N[l]);
			// cout << "N = " << N[l] << ", N_id = " << N_id[l] << " => tswap = " << tswap << ", pswap = " << pswap << "\n";
			for(unsigned i=0; i < fl_units[0][l]; i++) {
				for(unsigned j=0; j < fl_units[1][l]; j++) {
					for(unsigned k=0; k < fl_units[2][l]; k++) {
						for(unsigned d=0; d < 3; d++) {
							if(pswap >= rnd->rnd()) {
								if(((l == 0) && fill0[i][j][k][d]) || ((l == 1) && fill1[i][j][k][d])) N[l] --;
								if(l == 0) fill0[i][j][k][d] = tswap;
								else fill1[i][j][k][d] = tswap;
								if(tswap) N[l] ++;
							}
						}
					}
				}
			}
		}
		global_log->debug() << "Filling " << N[l] << " of 3*"
			<< fl_units[0][l] << "*" << fl_units[1][l] << "*" << fl_units[2][l]
			<< " = " << slots[l] << " slots (ideally " << N_id[l] << ")" << endl;
	}

	_simulation.initCanonical(10);
	_simulation.initStatistics(3003003);

	double loffset[3][2];
	loffset[0][0] = 0.1; loffset[1][0] = 0.3; loffset[2][0] = 0.1;
	loffset[0][1] = 0.1; loffset[1][1] = 0.8; loffset[2][1] = 0.1;
	double goffset[3][3] = {
		{0.0, 0.5, 0.5},
		{0.5, 0.0, 0.5},
		{0.5, 0.5, 0.0}
	};

	double v_avg = sqrt(3.0 * T);

	Component* component = _simulation.getEnsemble()->component(0);
	unsigned ID = 1;
	for(unsigned l=0; l < 2; l++) {
		for(unsigned i=0; i < fl_units[0][l]; i++) {
			for(unsigned j=0; j < fl_units[1][l]; j++) {
				for(unsigned k=0; k < fl_units[2][l]; k++) {
					for(unsigned d=0; d < 3; d++) {
						if(((l == 0) && fill0[i][j][k][d]) || ((l == 1) && fill1[i][j][k][d])) {
							double q[3];
							q[0] = i * fl_unit[0][l];
							q[1] = j * fl_unit[1][l];
							q[2] = k * fl_unit[2][l];
							for(int m=0; m < 3; m++) {
								q[m] += box[m]*loffset[m][l] + fl_unit[m][l]*goffset[m][d];
								q[m] += VARFRACTION * fl_unit[m][l] * (rnd->rnd() - 0.5);
								if(q[m] > box[m]) q[m] -= box[m];
							}
							double phi = 2*M_PI * rnd->rnd();
							double omega = 2*M_PI * rnd->rnd();

							double v[3];
							v[0] = v_avg*cos(phi)*cos(omega);
							v[1] = v_avg*cos(phi)*sin(omega);
							v[2] = v_avg*sin(phi);
							Molecule molecule(ID, component, q[0], q[1], q[2], v[0], v[1], v[2], 1, 0, 0, 0, 0, 0, 0);
							particleContainer->addParticle(molecule);
							ID++;
						}
					}
				}
			}
		}
	}
	domain->setGlobalTemperature(T);
	domain->setglobalNumMolecules(ID-1);
	domain->setglobalRho(ID / _simulation.getEnsemble()->V() );
	return ID;
}
