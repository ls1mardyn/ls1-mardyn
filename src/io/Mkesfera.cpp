/* the following macro has to be defined to use math constants in cmath */
#define _USE_MATH_DEFINES  1

#include "Mkesfera.h"

#include "particleContainer/ParticleContainer.h"
#include "particleContainer/LinkedCells.h"
#include "Domain.h"
#include "integrators/Integrator.h"
#include "integrators/Leapfrog.h"
#include "io/ResultWriter.h"
#include "molecules/Component.h"
#ifdef ENABLE_MPI
#include "parallel/DomainDecomposition.h"
#else
#include "parallel/DomainDecompDummy.h"
#endif
#include "Simulation.h"
#include "utils/Logger.h"
#include "utils/OptionParser.h"
#include "utils/Random.h"

#include <cmath>
#include <vector>

#define DT 0.002
#define TIME 20111102
#define VARFRACTION 0.07
#define BOXOVERLOAD 1.3333
#define AVGBIN 0.025

using namespace std;
using Log::global_log;

Mkesfera::Mkesfera(optparse::Values& options) {
	rho_i = options.get("density-1");  /* inner density */
	global_log->info() << "Inner density: " << rho_i << endl;
	R_i = options.get("droplet-radius-inner");
	global_log->info() << "Inner droplet radius: " << R_i << endl;
	rho_o = options.get("density-2"); /* outer density */
	global_log->info() << "Outer density: " << rho_o << endl;
	R_o = options.get("droplet-radius-outer");
	global_log->info() << "Outer droplet radius: " << R_o << endl;
	cutoff = options.get("cutoff-LJ");
	global_log->info() << "Cutoff radius: " << cutoff << endl;
	do_shift = options.get("shift_LJ");
	global_log->info() << "Shift LJ potential: " << do_shift << endl;
	T = options.get("temperature");
}

void Mkesfera::generate(Domain* domain, DomainDecompBase** domainDecomposition, Integrator** integrator, ParticleContainer** moleculeContainer, std::list< OutputBase* > &outputPlugins, Simulation* simulation) {

	unsigned fl_units;
	double rhomax = (rho_i > rho_o)? rho_i: rho_o;
	double N_boxes = 8.0*R_o*R_o*R_o * (BOXOVERLOAD*rhomax) / 3.0;
	fl_units = ceil(pow(N_boxes, 1.0/3.0));
	double fl_unit = 2.0*R_o / (double)fl_units;

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
					for(int d = 0; d < 3; d++) {
						double qrel = (idx[d] + goffset[d][p])*fl_unit - R_o;
						qq += qrel*qrel;
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

	simulation->initCanonical(10);
	simulation->initStatistics(3003003);
	simulation->setcutoffRadius(cutoff);
	simulation->setLJCutoff(cutoff);
	simulation->setTersoffCutoff(0.5);

	(*integrator) = new Leapfrog(DT);

	/*
	 * This part contains adapted code to enable most of the analysis features used in the original generator
	 * We leave this part out as this ist not tested in the current trunk version.
	 */
	/*
	unsigned xun = 1;
	unsigned yun = (unsigned)round(R_o / AVGBIN);
	unsigned zun = 1;

	simulation->_doRecordProfile = true;
	unsigned long profileRecordingTimesteps = 1;
	unsigned long profileOutputTimesteps = 500000;
	string profileOutputPrefix = simulation->getOutputPrefix();
	simulation->profileSettings(profileRecordingTimesteps, profileOutputTimesteps, profileOutputPrefix);
	domain->setupProfile(xun, yun, zun);

	esfera
	nprofiledComponent	1
	nomomentum 16384
	AlignCentre 333 0.003
	chemicalPotential 0 component 1 conduct " << (unsigned)round(N/3.0) << " tests every 1 steps
	Widom;
	*/

#ifdef ENABLE_MPI
	(*domainDecomposition) = (DomainDecompBase*) new DomainDecomposition();
#else
	(*domainDecomposition) = (DomainDecompBase*) new DomainDecompDummy();
#endif
	double box_max[3];
	box_max[0] = 2*R_o;
	box_max[1] = 2*R_o;
	box_max[2] = 2*R_o;
	domain->setCurrentTime(0.0);
	domain->setGlobalLength(0, box_max[0]);
	domain->setGlobalLength(1, box_max[1]);
	domain->setGlobalLength(2, box_max[2]);
	domain->setGlobalTemperature(T);
	domain->setglobalNumMolecules(N);

	double bBoxMin[3];
	double bBoxMax[3];
	for (int d = 0; d < 3; d++) {
	    bBoxMin[d] = (*domainDecomposition)->getBoundingBoxMin(d, domain);
	    bBoxMax[d] = (*domainDecomposition)->getBoundingBoxMax(d, domain);
	}
	(*moleculeContainer) = new LinkedCells(bBoxMin, bBoxMax, simulation->getcutoffRadius(), simulation->getLJCutoff(), 1);

	vector<Component> &components = domain->getComponents();
	components.resize(1);
	components[0].addLJcenter(0.0, 0.0, 0.0,  1.0, 1.0, 1.0, cutoff, (do_shift? 1: 0));
	components[0].setI11(0);
	components[0].setI22(0);
	components[0].setI33(0);
	domain->setepsilonRF(1e+10);

	double v_avg = sqrt(3.0 * T);

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
						Molecule molecule(ID, &components[0], q[0], q[1], q[2], v[0], v[1], v[2], 1, 0, 0, 0, 0, 0, 0);
						(*moleculeContainer)->addParticle(molecule);
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
}

