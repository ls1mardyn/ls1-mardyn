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

using namespace std;
using Log::global_log;

void MkesferaGenerator::readXML(XMLfileUnits& xmlconfig) {
#define MAX(a,b) (((a) >= (b)) (a) : (b));
	double boxLength = 0;
	for(int d = 0; d < 3; d++) {
		double length = _simulation.getEnsemble()->domain()->length(d);
		if(boxLength < length) {
			boxLength = length;
		}
	}
	global_log->info() << "Box length: " << boxLength << endl;
	R_o = boxLength / 2;

	xmlconfig.getNodeValueReduced("outer-density", rho_o);
	global_log->info() << "Outer density: " << rho_o << endl;

	xmlconfig.getNodeValueReduced("droplet/radius", R_i);
	global_log->info() << "Droplet radius: " << R_i << endl;
	xmlconfig.getNodeValueReduced("droplet/density", rho_i);
	global_log->info() << "Droplet density: " << rho_i << endl;
	for(int d = 0; d < 3; d++) {
		center[d] = R_o;
	}
	xmlconfig.getNodeValueReduced("droplet/center/x", center[0]);
	xmlconfig.getNodeValueReduced("droplet/center/y", center[1]);
	xmlconfig.getNodeValueReduced("droplet/center/z", center[2]);
	global_log->info() << "Droplet center: " << center[0] << ", " << center[0] << ", " << center[0] << endl;
}

long unsigned int MkesferaGenerator::readPhaseSpace(ParticleContainer* particleContainer, list< ChemicalPotential >* lmu, Domain* domain, DomainDecompBase* domainDecomp) {

	unsigned fl_units;
	double rhomax = (rho_i > rho_o)? rho_i: rho_o;
	double N_boxes = 8.0*R_o*R_o*R_o * (BOXOVERLOAD*rhomax) / 3.0;
	fl_units = ceil(pow(N_boxes, 1.0/3.0));
	double fl_unit = 2.0*R_o / (double)fl_units;

	double T = _simulation.getEnsemble()->T();
	global_log->info() << "Temperature: " << T << endl;

	double cutoff = _simulation.getcutoffRadius();
	Random* rnd = new Random();
	rnd->init(
		(int)(10000.0*R_o) - (int)(3162.3*cutoff)
		+ (int)(1000.0*T) - (int)(316.23*rho_i)
		+ (int)(100.0*R_i)
	);

   bool**** fill = new bool***[fl_units];
	for (unsigned int i = 0; i < fl_units; i++) {
		fill[i] = new bool**[fl_units];
		for (unsigned int j = 0; j < fl_units; j++) {
			fill[i][j] = new bool*[fl_units];
			for (unsigned int k = 0; k < fl_units; k++) {
					fill[i][j][k] = new bool[3];
				}
		}
	}
	unsigned slots = 3.0*fl_units*fl_units*fl_units;
	double boxdensity = (double)slots / (8.0*R_o*R_o*R_o);
	global_log->debug() << "Box density: " << boxdensity << " (unit cell: " << fl_unit << ")" << endl;
	double P_in = rho_i / boxdensity;
	double P_out = rho_o / boxdensity;
	global_log->debug() << "Insertion probability: " << P_in << " inside, " << P_out << " outside" << endl;

	double box_max[3];
	/* box min is assumed to be 0 */
	for(int d = 0; d < 3; d++) {
		box_max[d] = _simulation.getEnsemble()->domain()->length(d);
	}

	double goffset[3][3]; /**< relative koordinates of face centers */
	goffset[0][0] = 0.0; goffset[1][0] = 0.5; goffset[2][0] = 0.5;
	goffset[0][1] = 0.5; goffset[1][1] = 0.0; goffset[2][1] = 0.5;
	goffset[0][2] = 0.5; goffset[1][2] = 0.5; goffset[2][2] = 0.0;

	unsigned N = 0;
	unsigned idx[3];
	for(idx[0]=0; idx[0] < fl_units; idx[0]++) {
		for(idx[1]=0; idx[1] < fl_units; idx[1]++) {
			for(idx[2]=0; idx[2] < fl_units; idx[2]++) {
				for(int p=0; p < 3; p++) {
					double qq = 0.0;
					double q[3];
					for(int d = 0; d < 3; d++) {
						q[d]= (idx[d] + goffset[d][p])*fl_unit - center[d];
						q[d] = q[d] - round(q[d]/box_max[d])*box_max[d];
						qq += q[d]*q[d];
					}
					double tP = (qq > R_i*R_i)? P_out: P_in;
					bool tfill = (tP >= rnd->rnd());
					fill[idx[0]][idx[1]][idx[2]][p] = tfill;
					if(tfill) {
						N++;
					}
				}
			}
		}
	}
	global_log->debug() << "Filling " << N << " out of " << slots << " slots" << endl;
	global_log->debug() << "Density: " << N / (8.0*R_o*R_o*R_o) << endl;

	_simulation.initCanonical(10);
	_simulation.initStatistics(3003003);

	domain->setGlobalTemperature(T);
	domain->setglobalNumMolecules(N);
	domain->setglobalRho(N / _simulation.getEnsemble()->V() );

	double v_avg = sqrt(3.0 * T);

	Component* component = _simulation.getEnsemble()->component(0);
	unsigned ID = 1;
	for(idx[0]=0; idx[0] < fl_units; idx[0]++) {
		for(idx[1]=0; idx[1] < fl_units; idx[1]++) {
			for(idx[2]=0; idx[2] < fl_units; idx[2]++) {
				for(unsigned p=0; p < 3; p++) {
					if(fill[idx[0]][idx[1]][idx[2]][p]) {
						//global_log->debug() << "Inserting: " << idx[0] << "," << idx[1] << "," << idx[2] << "; " << p << endl;
						double q[3];
						for(int d=0; d < 3; d++)
						{
							q[d] = (idx[d] + VARFRACTION*(rnd->rnd() - 0.5) + goffset[d][p])*fl_unit;
							if(q[d] < 0.0) q[d] += 2.0*R_o;
							else if(q[d] > 2.0*R_o) q[d] -= 2.0*R_o;
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

	for (unsigned int i = 0; i < fl_units; i++) {
		for (unsigned int j = 0; j < fl_units; j++) {
			for (unsigned int k = 0; k < fl_units; k++) {
				delete[] fill[i][j][k];
			}
			delete[] fill[i][j];
		}
		delete[] fill[i];
	}
	delete[] fill;
	global_log->info() << "Inserted number of molecules: " << ID << endl;
	return ID;
}

