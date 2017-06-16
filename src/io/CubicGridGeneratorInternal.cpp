/*
 * CubicGridGeneratorInternal.cpp
 *
 *  Created on: Jun 9, 2017
 *      Author: seckler
 */

/* the following macro has to be defined to use math constants in cmath */
#define _USE_MATH_DEFINES  1
#include <cmath>

#include "CubicGridGeneratorInternal.h"
#include "utils/Logger.h"
#include "particleContainer/ParticleContainer.h"
#include "ensemble/GrandCanonical.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"

#include "utils/Random.h"



CubicGridGeneratorInternal::CubicGridGeneratorInternal() :
		_numMolecules(1000), _binaryMixture(false) {
}

void CubicGridGeneratorInternal::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("numMolecules", _numMolecules);
	global_log->info() << "numMolecules: " << _numMolecules << std::endl;
	xmlconfig.getNodeValue("binaryMixture", _binaryMixture);
	global_log->info() << "binaryMixture: " << _binaryMixture << std::endl;
}

unsigned long CubicGridGeneratorInternal::readPhaseSpace(ParticleContainer* particleContainer,
		std::list<ChemicalPotential>* lmu, Domain* domain, DomainDecompBase* domainDecomp) {
	global_simulation->startTimer("CUBIC_GRID_GENERATOR_INPUT");
	Log::global_log->info() << "Reading phase space file (CubicGridGenerator)." << std::endl;

	// create a body centered cubic layout, by creating by placing the molecules on the
	// vertices of a regular grid, then shifting that grid by spacing/2 in all dimensions.

	int numMoleculesPerDimension = ceil(pow((double) _numMolecules / 2.0, 1. / 3.));
	if (_binaryMixture) {
		global_simulation->getEnsemble()->getComponents()->at(1).updateMassInertia();
	}

	unsigned long int id = 0;

	double simBoxLength = domain->getGlobalLength(0);
	if(domain->getGlobalLength(0)!=simBoxLength || domain->getGlobalLength(0) != simBoxLength){
		Log::global_log->error() << "varying simBoxLength not yet supported for CubicGridGenerator" << std::endl;
		Simulation::exit(1);
	}
	double spacing = simBoxLength / numMoleculesPerDimension;
	double origin1 = spacing / 4.; // origin of the first DrawableMolecule
	double origin2 = spacing / 4. * 3.; // origin of the first DrawableMolecule

	int start_i = floor((domainDecomp->getBoundingBoxMin(0, domain) / simBoxLength) * numMoleculesPerDimension) - 1;
	int start_j = floor((domainDecomp->getBoundingBoxMin(1, domain) / simBoxLength) * numMoleculesPerDimension) - 1;
	int start_k = floor((domainDecomp->getBoundingBoxMin(2, domain) / simBoxLength) * numMoleculesPerDimension) - 1;

	int end_i = ceil((domainDecomp->getBoundingBoxMax(0, domain) / simBoxLength) * numMoleculesPerDimension) + 1;
	int end_j = ceil((domainDecomp->getBoundingBoxMax(1, domain) / simBoxLength) * numMoleculesPerDimension) + 1;
	int end_k = ceil((domainDecomp->getBoundingBoxMax(2, domain) / simBoxLength) * numMoleculesPerDimension) + 1;

	// only for console output
	int percentageRead = 0;
	double percentage = 1.0 / (end_i - start_i) * 100.0;

    const int blocksize = 4;

	for (int i = start_i; i < end_i; i+=blocksize) {
        for (int j = start_j; j < end_j; j+=blocksize) {
            for (int k = start_k; k < end_k; k+=blocksize) {
                for (int ii = i; ii < i+blocksize and ii < end_i; ii++) {
                    for (int jj = j; jj < j+blocksize and jj < end_j; jj++) {
                        for (int kk = k; kk < k+blocksize and kk < end_k; kk++) {

                            double x1 = origin1 + ii * spacing;
                            double y1 = origin1 + jj * spacing;
                            double z1 = origin1 + kk * spacing;
                            if (domainDecomp->procOwnsPos(x1, y1, z1, domain)) {
                                addMolecule(x1, y1, z1, id, particleContainer);
                                id++;
                            }

                            double x2 = origin2 + ii * spacing;
                            double y2 = origin2 + jj * spacing;
                            double z2 = origin2 + kk * spacing;
                            if (domainDecomp->procOwnsPos(x2, y2, z2, domain)) {
                                addMolecule(x2, y2, z2, id, particleContainer);
                                id++;
                            }
                        }
                    }
                }
            }
        }
        if ((int) (i * percentage) > percentageRead) {
            percentageRead = i * percentage;
            Log::global_log->info() << "Finished reading molecules: " << (percentageRead) << "%\r" << std::flush;
        }
    }

	domainDecomp->collCommInit(1);
	domainDecomp->collCommAppendUnsLong(id); //number of local molecules
	domainDecomp->collCommScanSum();
	unsigned long idOffset = domainDecomp->collCommGetUnsLong() - id;
	domainDecomp->collCommFinalize();
	// fix ID's to be unique:
	for (auto mol = particleContainer->iteratorBegin(); mol != particleContainer->iteratorEnd(); ++mol) {
		mol->setid(mol->id() + idOffset);
	}
	//std::cout << domainDecomp->getRank()<<": #num local molecules:" << id << std::endl;
	//std::cout << domainDecomp->getRank()<<": offset:" << idOffset << std::endl;

	removeMomentum(particleContainer, *(global_simulation->getEnsemble()->getComponents()));
	domain->evaluateRho(particleContainer->getNumberOfParticles(), domainDecomp);
	Log::global_log->info() << "Calculated Rho=" << domain->getglobalRho() << std::endl;
	global_simulation->stopTimer("CUBIC_GRID_GENERATOR_INPUT");
	global_simulation->setOutputString("CUBIC_GRID_GENERATOR_INPUT", "Initial IO took:                 ");
	Log::global_log->info() << "Initial IO took:                 "
			<< global_simulation->getTime("CUBIC_GRID_GENERATOR_INPUT") << " sec" << std::endl;

	return id + idOffset;
}

void CubicGridGeneratorInternal::addMolecule(double x, double y, double z, unsigned long id,
		ParticleContainer* particleContainer) {
	std::vector<double> velocity = getRandomVelocity(global_simulation->getEnsemble()->T());

	//double orientation[4] = {1, 0, 0, 0}; // default: in the xy plane
	// rotate by 30° along the vector (1/1/0), i.e. the angle bisector of x and y axis
	// o = cos 30° + (1 1 0) * sin 15°
	double orientation[4];
	getOrientation(15, 10, orientation);

	int componentType = 0;
	if (_binaryMixture) {
		componentType = randdouble(0, 1.999999);
	}

	double I[3] = { 0., 0., 0. };
	I[0] = global_simulation->getEnsemble()->getComponents()->at(0).I11();
	I[1] = global_simulation->getEnsemble()->getComponents()->at(0).I22();
	I[2] = global_simulation->getEnsemble()->getComponents()->at(0).I33();
	/*****  Copied from animake - initialize anular velocity *****/
	double w[3];
	for (int d = 0; d < 3; d++) {
		w[d] = (I[d] == 0) ?
				0.0 : ((randdouble(0, 1) > 0.5) ? 1 : -1) * sqrt(2.0 * randdouble(0, 1) * global_simulation->getEnsemble()->T() / I[d]);
		double fs_2_mardyn = 0.030619994;
		w[d] = w[d] * fs_2_mardyn;
	}
	/************************** End Copy **************************/

	Molecule m(id, &(global_simulation->getEnsemble()->getComponents()->at(componentType)), x, y, z, // position
			velocity[0], -velocity[1], velocity[2], // velocity
			orientation[0], orientation[1], orientation[2], orientation[3], w[0], w[1], w[2]);
	particleContainer->addParticle(m);
}

void CubicGridGeneratorInternal::removeMomentum(ParticleContainer* particleContainer,
		const std::vector<Component>& components) {
	double mass = 0.;
	double mass_sum = 0.;
	double momentum_sum[3] = { 0., 0., 0. };

	ParticleIterator molecule = particleContainer->iteratorBegin();
	while (molecule != particleContainer->iteratorEnd()) {
		mass = components[molecule->componentid()].m();
		mass_sum = mass_sum + mass;
		momentum_sum[0] = momentum_sum[0] + mass * molecule->v(0);
		momentum_sum[1] = momentum_sum[1] + mass * molecule->v(1);
		momentum_sum[2] = momentum_sum[2] + mass * molecule->v(2);
		++molecule;
	}

	double momentum_sub0 = momentum_sum[0] / mass_sum;
	double momentum_sub1 = momentum_sum[1] / mass_sum;
	double momentum_sub2 = momentum_sum[2] / mass_sum;

	molecule = particleContainer->iteratorBegin();
	while (molecule != particleContainer->iteratorEnd()) {
		molecule->vsub(momentum_sub0, momentum_sub1, momentum_sub2);
		++molecule;
	}

	//test
	momentum_sum[0] = 0.;
	momentum_sum[1] = 0.;
	momentum_sum[2] = 0.;

	molecule = particleContainer->iteratorBegin();
	while (molecule != particleContainer->iteratorEnd()) {
		mass = components[molecule->componentid()].m();
		mass_sum = mass_sum + mass;
		momentum_sum[0] = momentum_sum[0] + mass * molecule->v(0);
		momentum_sum[1] = momentum_sum[1] + mass * molecule->v(1);
		momentum_sum[2] = momentum_sum[2] + mass * molecule->v(2);
		++molecule;
	}

	//printf("momentum_sum[0] from removeMomentum is %lf\n", momentum_sum[0]);
	//printf("momentum_sum[1] from removeMomentum is %lf\n", momentum_sum[1]);
	//printf("momentum_sum[2] from removeMomentum is %lf\n", momentum_sum[2]);
}

void CubicGridGeneratorInternal::getOrientation(int base, int delta, double orientation[4]) {
	double offset = randdouble(-delta / 2., delta / 2.) / 180. * M_PI;
	double rad = base / 180. * M_PI;
	double angle = rad + offset;

	double cosinePart = cos(angle);
	double sinePart = sin(angle);

	double length = sqrt(cosinePart * cosinePart + 2 * (sinePart * sinePart));
	orientation[0] = cosinePart / length;
	orientation[1] = sinePart / length;
	orientation[2] = sinePart / length;
	orientation[3] = 0;
}

std::vector<double> CubicGridGeneratorInternal::getRandomVelocity(double temperature) const {
	std::vector<double> v_;
	v_.resize(3);

	// Velocity
	for (int dim = 0; dim < 3; dim++) {
		v_[dim] = randdouble(-0.5, 0.5);
	}
	double dotprod_v = 0;
	for (unsigned int i = 0; i < v_.size(); i++) {
		dotprod_v += v_[i] * v_[i];
	}
	// Velocity Correction
	double vCorr = sqrt(3.0 * temperature / dotprod_v);
	for (unsigned int i = 0; i < v_.size(); i++) {
		v_[i] *= vCorr;
	}

	return v_;
}
