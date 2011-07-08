#include "io/OneCLJGenerator.h"
#include "io/DropletPlacement.h"

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"

#include <cmath>
#include <iostream>

const double OneCLJGenerator::eps = 1.0;
const double OneCLJGenerator::sigma = 1.0;
const double OneCLJGenerator::mass = 1.0;

OneCLJGenerator::OneCLJGenerator(string mode, unsigned long N, double T) : _moleculeCountOffset(0) {
	_mode = mode;
	_temperature = T;
	_numberOfMolecules = N;
	_simBoxLength.resize(3);
}

void OneCLJGenerator::setHomogeneuosParameter(double rho) {
	_rho = rho;

	// calculate Size of the Simulation Box
	_simBoxLength[0] = pow(_numberOfMolecules / _rho, (1.0 / 3.0));
	_simBoxLength[1] = pow(_numberOfMolecules / _rho, (1.0 / 3.0));
	_simBoxLength[2] = pow(_numberOfMolecules / _rho, (1.0 / 3.0));
}

void OneCLJGenerator::setClusterParameters(double gasDensity, double fluidDensity, double fluidVolumePercent,
        double maxSphereVolume, double numSphereSizes) {
	_gasDensity = gasDensity;
	_fluidDensity = fluidDensity;

	_fluidVolume = fluidVolumePercent;
	_maxSphereVolume = maxSphereVolume;
	_numSphereSizes = numSphereSizes;

	_rho = fluidDensity * fluidVolumePercent / 100.0 + gasDensity * (1 - fluidVolumePercent / 100.0);
	_simBoxLength[0] = pow(_numberOfMolecules / _rho, (1.0 / 3.0));
	_simBoxLength[1] = pow(_numberOfMolecules / _rho, (1.0 / 3.0));
	_simBoxLength[2] = pow(_numberOfMolecules / _rho, (1.0 / 3.0));

}

void OneCLJGenerator::readPhaseSpaceHeader(Domain* domain, double timestepLength) {
	vector<Component>& dcomponents = domain->getComponents();
	domain->setCurrentTime(0);
	domain->disableComponentwiseThermostat();
	domain->setGlobalTemperature(_temperature);
	domain->setGlobalLength(0, _simBoxLength[0]);
	domain->setGlobalLength(1, _simBoxLength[1]);
	domain->setGlobalLength(2, _simBoxLength[2]);
	dcomponents.resize(1);
	dcomponents[0].setID(0);
	// add a LJ-center x/y/z = (0/0/0), m=1.0 eps=1.0 sigma=1.0, rc=1.0, truncated=false
	// I think, this "rc" parameter is not used as truncate is false.
	dcomponents[0].addLJcenter(0.0, 0.0, 0.0, mass, eps, sigma, 1.0, false);

	domain->setepsilonRF(1e+10);
}

//! @brief read the actual phase space information
unsigned long OneCLJGenerator::readPhaseSpace(ParticleContainer* particleContainer, std::list<ChemicalPotential>* lmu,
        Domain* domain, DomainDecompBase* domainDecomp) {

	vector<double> bBoxMin;
	vector<double> bBoxMax;

	// This is rather a hack to make sure that the ids of the molecules for each process are unique
	_moleculeCountOffset = (2 * domainDecomp->getRank()) * _numberOfMolecules / domainDecomp->getNumProcs();
	Log::global_log->info() << "MoleculeCountOffset=" << _moleculeCountOffset << endl;

	bBoxMin.resize(3);
	bBoxMax.resize(3);
	for (int i = 0; i < 3; i++) {
		bBoxMin[i] = domainDecomp->getBoundingBoxMin(i, domain);
		bBoxMax[i] = domainDecomp->getBoundingBoxMax(i, domain);
	}

	if (_mode == "Homogeneous") {
		Log::global_log->info() << "OneCLJGenerator  generating molecules. mode: " << _mode << " rho " << _rho << " T "
		        << _temperature << " #molecules " << _numberOfMolecules << endl;
		createHomogeneousDist(particleContainer, bBoxMin, bBoxMax, domain, domainDecomp);
	} else if (_mode == "Cluster") {
		Log::global_log->info() << "OneCLJGenerator  generating cluster distribution. " << " T " << _temperature
		        << " #molecules " << _numberOfMolecules << " rho_gas " << _gasDensity << " rho_fluid " << _fluidDensity
		        << endl;
		createClusters(particleContainer, bBoxMin, bBoxMax, domain, domainDecomp);
	} else {
		Log::global_log->error() << "Unknown mode \"" << _mode << "\"!" << endl;
		exit(1);
	}

	vector<Component>& dcomponents = domain->getComponents();
	vector<unsigned long> partsPerComp;
	partsPerComp.resize(1);
	particleContainer->update();
	particleContainer->deleteOuterParticles();
	domain->setglobalNumMolecules(domainDecomp->countMolecules(particleContainer, partsPerComp));

	for (unsigned int i = 0; i < partsPerComp.size(); i++) {
		dcomponents[i].setNumMolecules(partsPerComp[i]);
		domain->setglobalRotDOF(partsPerComp[i] * dcomponents[i].getRotationalDegreesOfFreedom());
	}
	domain->setglobalRho(domain->getglobalNumMolecules() / (_simBoxLength[0] * _simBoxLength[1] * _simBoxLength[2]));
	return domain->getglobalNumMolecules();
}


void OneCLJGenerator::createHomogeneousDist(ParticleContainer* particleContainer, vector<double> &bBoxMin, vector<
        double> &bBoxMax, Domain* domain, DomainDecompBase* domainDecomp) {

	vector<int> globalFccCells;
	vector<int> localFccCellsMin;
	vector<int> localFccCellsMax;
	vector<double> fccCellLength;

	globalFccCells.resize(3);
	localFccCellsMin.resize(3);
	localFccCellsMax.resize(3);
	fccCellLength.resize(3);

	//_numberOfMolecules=1;
	for (int dim = 0; dim < 3; dim++) {
		globalFccCells[dim] = (int) ceil(pow(_rho / 4.0, 1.0 / 3.0) * _simBoxLength[dim]);
		//_numberOfMolecules *= globalFccCells[dim];
	}
	//_numberOfMolecules *= 4;

	for (int dim = 0; dim < 3; dim++) {
		fccCellLength[dim] = _simBoxLength[dim] / globalFccCells[dim];
		localFccCellsMin[dim] = (int) floor(bBoxMin[dim] / fccCellLength[dim]);
		localFccCellsMax[dim] = (int) floor(bBoxMax[dim] / fccCellLength[dim]);
		//_numberOfMolecules *= globalFccCells[dim];
	}

	vector<vector<double> > fccOffsets;
	fccOffsets.resize(4);
	for (int i = 0; i < 4; i++) {
		fccOffsets[i].resize(3);
		for (int j = 0; j < 3; j++) {
			fccOffsets[i][j] = 0.0;
		}
	}
	fccOffsets[1][0] = 0.5 * fccCellLength[0];
	fccOffsets[1][1] = 0.5 * fccCellLength[1];
	fccOffsets[2][0] = 0.5 * fccCellLength[0];
	fccOffsets[2][2] = 0.5 * fccCellLength[2];
	fccOffsets[3][1] = 0.5 * fccCellLength[1];
	fccOffsets[3][2] = 0.5 * fccCellLength[2];

	vector<double> r_;
	r_.resize(3);

	unsigned long int molCount = _moleculeCountOffset;
	for (int fcc = 0; fcc < 4; fcc++) {
		for (int iz = localFccCellsMin[2]; iz <= localFccCellsMax[2]; iz++) {
			for (int iy = localFccCellsMin[1]; iy <= localFccCellsMax[1]; iy++) {
				for (int ix = localFccCellsMin[0]; ix <= localFccCellsMax[0]; ix++) {

					// Position
					r_[0] = ix * (fccCellLength[0]) + fccOffsets[fcc][0] + 0.000000001;
					r_[1] = iy * (fccCellLength[1]) + fccOffsets[fcc][1] + 0.000000001;
					r_[2] = iz * (fccCellLength[2]) + fccOffsets[fcc][2] + 0.000000001;

					// if the position is not in the domain of this proc,
					// the molecule must not be created
					if (domainDecomp->procOwnsPos(r_[0], r_[1], r_[2], domain)) {
						addParticle(molCount, r_[0], r_[1], r_[2], particleContainer, domain, domainDecomp);
						molCount++;
					}
				}
			}
		}
	}
}

void OneCLJGenerator::createClusters(ParticleContainer* particleContainer, vector<double> &bBoxMin,
        vector<double> &bBoxMax, Domain* domain, DomainDecompBase* domainDecomp) {

	readLocalClusters(domain, domainDecomp);
	vector<int> globalFccCells;
	vector<int> clusterFccCellsMin;
	vector<int> clusterFccCellsMax;
	vector<double> fccCellLength;

	globalFccCells.resize(3);
	clusterFccCellsMin.resize(3);
	clusterFccCellsMax.resize(3);
	fccCellLength.resize(3);

	// fluid properties
	for (int dim = 0; dim < 3; dim++) {
		globalFccCells[dim] = (int) ceil(pow(_fluidDensity / 4., 1. / 3.) * _simBoxLength[dim]);
		fccCellLength[dim] = _simBoxLength[dim] / globalFccCells[dim];
	}
	vector<vector<double> > fccOffsets;
	fccOffsets.resize(4);
	for (int i = 0; i < 4; i++) {
		fccOffsets[i].resize(3);
		for (int j = 0; j < 3; j++) {
			fccOffsets[i][j] = 0.0;
		}
	}
	fccOffsets[1][0] = 0.5 * fccCellLength[0];
	fccOffsets[1][1] = 0.5 * fccCellLength[1];
	fccOffsets[2][0] = 0.5 * fccCellLength[0];
	fccOffsets[2][2] = 0.5 * fccCellLength[2];
	fccOffsets[3][1] = 0.5 * fccCellLength[1];
	fccOffsets[3][2] = 0.5 * fccCellLength[2];

	double securityOffset = 0.5 * fccCellLength[0];

	vector<double> clusterPos;
	clusterPos.resize(3);
	double radius;

	vector<double> r_;
	r_.resize(3);

	unsigned long int molCount = _moleculeCountOffset;
	for (unsigned int cluster = 0; cluster < _localClusters.size(); cluster++) {
		radius = _localClusters[cluster][3];
		for (int dim = 0; dim < 3; dim++) {
			clusterPos[dim] = _localClusters[cluster][dim];
			clusterFccCellsMin[dim] = (int) floor((clusterPos[dim] - radius) / fccCellLength[dim]);
			clusterFccCellsMax[dim] = (int) floor((clusterPos[dim] + radius) / fccCellLength[dim]);
		}

		for (int fcc = 0; fcc < 4; fcc++) {
			for (int iz = clusterFccCellsMin[2]; iz <= clusterFccCellsMax[2]; iz++) {
				for (int iy = clusterFccCellsMin[1]; iy <= clusterFccCellsMax[1]; iy++) {
					for (int ix = clusterFccCellsMin[0]; ix <= clusterFccCellsMax[0]; ix++) {
						// Position
						r_[0] = ix * (fccCellLength[0]) + fccOffsets[fcc][0] + 0.000000001;
						r_[1] = iy * (fccCellLength[1]) + fccOffsets[fcc][1] + 0.000000001;
						r_[2] = iz * (fccCellLength[2]) + fccOffsets[fcc][2] + 0.000000001;

						if (sqrt(pow(r_[0] - clusterPos[0], 2.0) + pow(r_[1] - clusterPos[1], 2.0) + pow(r_[2]
						        - clusterPos[2], 2.0)) > radius) {
							// molecule is too far away from the center of the cluster
							continue;
						}

						if (belongsToPreviousCluster(r_[0], r_[1], r_[2], cluster)) {
							// some other cluster already created this molecule
							continue;
						}

						// if the position is not in the domain of this proc,
						// the molecule must not be created
						if (domainDecomp->procOwnsPos(r_[0], r_[1], r_[2], domain)) {
							addParticle(molCount, r_[0], r_[1], r_[2], particleContainer, domain, domainDecomp);
							molCount++;
						}
					}
				}
			}
		}
	}

	// gas properties
	vector<int> fccCellsMin;
	vector<int> fccCellsMax;
	fccCellsMin.resize(3);
	fccCellsMax.resize(3);

	for (int dim = 0; dim < 3; dim++) {
		globalFccCells[dim] = (int) ceil(pow(_gasDensity / 4., 1. / 3.) * _simBoxLength[dim]);
		fccCellLength[dim] = _simBoxLength[dim] / globalFccCells[dim];
		fccCellsMin[dim] = (int) floor(bBoxMin[dim] / fccCellLength[dim]);
		fccCellsMax[dim] = (int) floor(bBoxMax[dim] / fccCellLength[dim]);
	}

	fccOffsets[1][0] = 0.5 * fccCellLength[0];
	fccOffsets[1][1] = 0.5 * fccCellLength[1];
	fccOffsets[2][0] = 0.5 * fccCellLength[0];
	fccOffsets[2][2] = 0.5 * fccCellLength[2];
	fccOffsets[3][1] = 0.5 * fccCellLength[1];
	fccOffsets[3][2] = 0.5 * fccCellLength[2];

	// GAS
	for (int fcc = 0; fcc < 4; fcc++) {
		for (int iz = fccCellsMin[2]; iz <= fccCellsMax[2]; iz++) {
			for (int iy = fccCellsMin[1]; iy <= fccCellsMax[1]; iy++) {
				for (int ix = fccCellsMin[0]; ix <= fccCellsMax[0]; ix++) {
					// Position
					r_[0] = ix * (fccCellLength[0]) + fccOffsets[fcc][0] + 0.000000001;
					r_[1] = iy * (fccCellLength[1]) + fccOffsets[fcc][1] + 0.000000001;
					r_[2] = iz * (fccCellLength[2]) + fccOffsets[fcc][2] + 0.000000001;

					if (closeToAnyCluster(r_[0], r_[1], r_[2], securityOffset)) {
						// some other cluster already created this molecule
						continue;
					}

					// if the position is not in the domain of this proc,
					// the molecule must not be created
					if (domainDecomp->procOwnsPos(r_[0], r_[1], r_[2], domain)) {
						addParticle(molCount, r_[0], r_[1], r_[2], particleContainer, domain, domainDecomp);
						molCount++;
					}
				}
			}
		}
	}
}


void OneCLJGenerator::addParticle(int id, double x, double y, double z, ParticleContainer* particleContainer,
        Domain* domain, DomainDecompBase* domainDecomp) {
	vector<Component>& dcomponents = domain->getComponents();
	vector<double> v_;
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
	double vCorr = sqrt(3.0 * _temperature / dotprod_v);
	for (unsigned int i = 0; i < v_.size(); i++) {
		v_[i] *= vCorr;
	}

	Molecule m1 = Molecule(id, 0, x, y, z, v_[0], v_[1], v_[2], 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &dcomponents);
	particleContainer->addParticle(m1);
	//dcomponents[0].incrnumMolecules();
	//domain->setglobalRotDOF(dcomponents[0].rot_dof()+domain->getglobalRotDOF());
}

void OneCLJGenerator::readLocalClusters(Domain* domain, DomainDecompBase* domainDecomp){

	DropletPlacement dropletPlacement(_fluidVolume, _maxSphereVolume, _numSphereSizes);
	vector<DropletPlacement::Droplet> droplets = dropletPlacement.generateDroplets();

	vector<double> shiftedSphere;
	vector<double> sphere;
	double distanceToDomain;
	sphere.resize(4); // x y z r
	shiftedSphere.resize(4); // x y z r
	for (unsigned int count = 0; count < droplets.size(); count++) {
		for (int i = 0; i < 3; i++) {
			sphere[i] = droplets[count]._center[i] * _simBoxLength[i];
		}
		sphere[3] = droplets[count]._radius * _simBoxLength[0];
		shiftedSphere[3] = sphere[3];

		for (int ix = -1; ix <= 1; ix++) {
			for (int iy = -1; iy <= 1; iy++) {
				for (int iz = -1; iz <= 1; iz++) {
					shiftedSphere[0] = sphere[0] + ix * _simBoxLength[0];
					shiftedSphere[1] = sphere[1] + iy * _simBoxLength[1];
					shiftedSphere[2] = sphere[2] + iz * _simBoxLength[2];
					distanceToDomain = domainDecomp->guaranteedDistance(shiftedSphere[0], shiftedSphere[1],
					        shiftedSphere[2], domain);
					if (distanceToDomain <= shiftedSphere[3]) { // sphere[3] = radius
						// reduce number of spheres for domains with periodic boundary
						bool tooFar = false;
						for (int i = 0; i < 3; i++) {
							if (shiftedSphere[i] < -sphere[3] || shiftedSphere[i] > _simBoxLength[i] + sphere[3])
								tooFar = true;
						}
						if (not (tooFar)) {
							_localClusters.push_back(shiftedSphere);
						}
					}
				}
			}
		}
	}
	//clusterStream.close();
}

bool OneCLJGenerator::belongsToPreviousCluster(double x, double y, double z, int clusterid) {
	bool belongsToOther = false;
	for (int i = 0; i < clusterid; i++) {
		if (sqrt(pow(x - _localClusters[i][0], 2.0) + pow(y - _localClusters[i][1], 2.0) + pow(
		        z - _localClusters[i][2], 2.0)) <= _localClusters[i][3]) {
			belongsToOther = true;
			break;
		}
	}
	return belongsToOther;
}

bool OneCLJGenerator::closeToAnyCluster(double x, double y, double z, double offset) {
	bool belongsToOther = false;
	for (unsigned int i = 0; i < _localClusters.size(); i++) {
		if (sqrt(pow(x - _localClusters[i][0], 2.0) + pow(y - _localClusters[i][1], 2.0) + pow(
		        z - _localClusters[i][2], 2.0)) <= _localClusters[i][3] + offset) {
			belongsToOther = true;
			break;
		}
	}
	return belongsToOther;
}

double OneCLJGenerator::randdouble(double a, double b) {
	return a + rand() * (b - a) / (RAND_MAX);
}

