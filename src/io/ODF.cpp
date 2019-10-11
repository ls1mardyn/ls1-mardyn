// Created by Joshua Marx 08/2019

#include "ODF.h"

void ODF::readXML(XMLfileUnits& xmlconfig) {
	global_log->debug() << "[ODF] enabled. Dipole orientations must be set to [0 0 1]!" << std::endl;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "[ODF] Write frequency: " << _writeFrequency << endl;
	xmlconfig.getNodeValue("initstatistics", _initStatistics);
	global_log->info() << "[ODF] Init Statistics: " << _initStatistics << endl;
	xmlconfig.getNodeValue("recordingtimesteps", _recordingTimesteps);
	global_log->info() << "[ODF] Recording Timesteps: " << _recordingTimesteps << endl;
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "[ODF] Output prefix: " << _outputPrefix << endl;
	xmlconfig.getNodeValue("phi1increments", _phi1Increments);
	global_log->info() << "[ODF] Phi1 increments: " << _phi1Increments << endl;
	xmlconfig.getNodeValue("phi2increments", _phi2Increments);
	global_log->info() << "[ODF] Phi2 increments: " << _phi2Increments << endl;
	xmlconfig.getNodeValue("gammaincrements", _gammaIncrements);
	global_log->info() << "[ODF] Gamma increments: " << _gammaIncrements << endl;
	xmlconfig.getNodeValue("shellcutoff1", _shellCutOff[0]);
	global_log->info() << "[ODF] Shell cutoff component one: " << _shellCutOff[0] << endl;
	xmlconfig.getNodeValue("shellcutoff2", _shellCutOff[1]);
	global_log->info() << "[ODF] Shell cutoff component two: " << _shellCutOff[1] << endl;
	xmlconfig.getNodeValue(
		"applyshellmixingrule",
		_mixingRule);  // if mixingrule = 1 then the shell cutoff for heterogenous pairings is calculated as half the
					   // diameter of the central atom plus a full diameter of the surrounding atom
	global_log->info() << "[ODF] Shell mixing rule: " << _mixingRule << endl;
}

void ODF::init(ParticleContainer* /*particleContainer*/, DomainDecompBase* /*domainDecomp*/, Domain* /*domain*/) {
	std::vector<Component>* components = global_simulation->getEnsemble()->getComponents();
	this->_numComponents = components->size();
	std::vector<unsigned int> isDipole(this->_numComponents);
	unsigned int numPairs = 0;

	for (unsigned int i = 0; i < this->_numComponents; ++i) {
		Component& ci = (*components)[i];
		isDipole[i] = ci.numDipoles();

		if (isDipole[i] == 1) {
			numPairs++;
		}
	}

	this->_numPairs = numPairs * numPairs;
	this->_numElements = this->_phi1Increments * this->_phi2Increments * this->_gammaIncrements + 1;
	global_log->info() << "ODF arrays contains " << this->_numElements << " elements each for " << this->_numPairs
					   << "pairings" << endl;
	this->_ODF11.resize(_numElements);
	this->_ODF12.resize(_numElements);
	this->_ODF21.resize(_numElements);
	this->_ODF22.resize(_numElements);
	this->_localODF11.resize(_numElements);
	this->_localODF12.resize(_numElements);
	this->_localODF21.resize(_numElements);
	this->_localODF22.resize(_numElements);

	if (this->_numPairs < 1) {
		global_log->error() << "No components with dipoles. ODFs not being calculated!" << endl;
	} else if (this->_numPairs > 4) {
		global_log->error()
			<< "Number of pairings for ODF calculation too high. Current maximum number of ODF pairings is 4." << endl;
	}

	this->reset();
}

void ODF::afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) {
	if (simstep > this->_initStatistics && simstep % this->_recordingTimesteps == 0) {
		this->record(particleContainer, global_simulation->getDomain(), domainDecomp, simstep);
	}
}

void ODF::endStep(ParticleContainer* /*particleContainer*/, DomainDecompBase* domainDecomp, Domain* domain,
				  unsigned long simstep) {
	if (simstep > this->_initStatistics && simstep % this->_writeFrequency == 0) {
		this->collect(domainDecomp);
		this->output(domain, simstep);
		this->reset();
	}
}

void ODF::reset() {
	global_log->info() << "[ODF] resetting data sets" << endl;

	for (unsigned long i = 0; i < this->_numElements; i++) {
		this->_ODF11[i] = 0;
		this->_ODF12[i] = 0;
		this->_ODF21[i] = 0;
		this->_ODF22[i] = 0;
		this->_localODF11[i] = 0;
		this->_localODF12[i] = 0;
		this->_localODF21[i] = 0;
		this->_localODF22[i] = 0;
	}
}

void ODF::record(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* /*domainDecomp*/,
				 unsigned long /*simstep*/) {
	for (auto it = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it) {
		if (it->numDipoles() == 1) {
			this->calculateOrientation(particleContainer, domain, *it);
		}
	}
}

void ODF::calculateOrientation(ParticleContainer* particleContainer, Domain* domain, const Molecule& mol1) {
	std::array<double, 3> r1{mol1.r(0), mol1.r(1), mol1.r(2)};

	std::array<double, 3> upVec1{};
	{
		std::array<double, 4> q1{mol1.q().qw(), mol1.q().qx(), mol1.q().qy(), mol1.q().qz()};
		upVec1[0] = 2 * (q1[1] * q1[3] + q1[0] * q1[2]);
		upVec1[1] = 2 * (q1[2] * q1[3] - q1[0] * q1[1]);
		upVec1[2] = 1 - 2 * (q1[1] * q1[1] + q1[2] * q1[2]);
	}
	// records mutual orientation of particle pairs

	// TODO Implement rotation matrices to calculate orientations for dipole direction unit vectors other than [0 0 1];
	double Q2[4], upVec2[3], r12[3], dist1D[3], auxVec1[3], auxVec2[3], projection1[3], projection2[3];

	double cosPhi1, cosPhi2, cosGamma12, Gamma12, norm1, norm2, normr12, shellcutoff;
	double roundingThreshold = 0.0001;
	unsigned long ind_phi1, ind_phi2, ind_gamma, elementID;
	unsigned maximum;
	bool bool1, bool2, bool3;
	auto cid = mol1.getComponentLookUpID();
	shellcutoff = this->_shellCutOff[cid];
	std::array<double, 3> boxMin{r1[0] - shellcutoff, r1[1] - shellcutoff, r1[2] - shellcutoff};
	std::array<double, 3> boxMax{r1[0] + shellcutoff, r1[1] + shellcutoff, r1[2] + shellcutoff};
	for (auto it = particleContainer->regionIterator(boxMin.data(), boxMax.data(), ParticleIterator::ALL_CELLS);
		 it.isValid(); ++it) {
		// inner loop over all particles. this should only loop over the neighbours of the particle from the
		// outer loop to improve efficiency

		if (it->numDipoles() == 1) {
			double distanceSquared = 0.;
			// minimum image convention
			for (int i = 0; i < 3; i++) {
				dist1D[i] = r1[i] - it->r(i);
				if (abs(dist1D[i]) > 0.5 * domain->getGlobalLength(i)) {
					dist1D[i] = domain->getGlobalLength(i) - abs(dist1D[i]);
				}
				distanceSquared += dist1D[i] * dist1D[i];
			}

			if (this->_mixingRule == 1) {
				shellcutoff = 1. / 3 * this->_shellCutOff[cid] + this->_shellCutOff[it->getComponentLookUpID()];
			}

			if (distanceSquared < shellcutoff * shellcutoff && mol1.getID() != it->getID()) {
				normr12 = 0.;

				// reading second molecule's quaternions

				Q2[0] = it->q().qw();
				Q2[1] = it->q().qx();
				Q2[2] = it->q().qy();
				Q2[3] = it->q().qz();

				cosPhi1 = 0.;
				cosPhi2 = 0.;
				cosGamma12 = 0.;
				norm1 = 0.;
				norm2 = 0.;
				ind_phi1 = 0;
				ind_phi2 = 0;
				ind_gamma = 0;
				bool1 = false;
				bool2 = false;
				bool3 = false;

				// calculate distance vector between molecules
				for (unsigned i = 0; i < 3; i++) {
					r12[i] = it->r(i) - r1[i];
					if (r12[i] > 0.5 * domain->getGlobalLength(i)) {
						r12[i] = -(domain->getGlobalLength(i) - r12[i]);
					} else if (r12[i] < -0.5 * domain->getGlobalLength(i)) {
						r12[i] = -(r12[i] - domain->getGlobalLength(i));
					}
					normr12 += r12[i] * r12[i];
				}

				normr12 = sqrt(normr12);

				for (double& i : r12) {
					i /= normr12;
				}

				// calculate vectors pointing in the direction defined by the dipole
				upVec2[0] = 2 * (Q2[1] * Q2[3] + Q2[0] * Q2[2]);
				upVec2[1] = 2 * (Q2[2] * Q2[3] - Q2[0] * Q2[1]);
				upVec2[2] = 1 - 2 * (Q2[1] * Q2[1] + Q2[2] * Q2[2]);

				// calculate projection of the vectors onto plane perpendicular to the distance vector with cross
				// product for calculation of the torque angle gamma
				auxVec1[0] = upVec1[1] * r12[2] - upVec1[2] * r12[1];
				auxVec1[1] = upVec1[2] * r12[0] - upVec1[0] * r12[2];
				auxVec1[2] = upVec1[0] * r12[1] - upVec1[1] * r12[0];

				projection1[0] = r12[1] * auxVec1[2] - r12[2] * auxVec1[1];
				projection1[1] = r12[2] * auxVec1[0] - r12[0] * auxVec1[2];
				projection1[2] = r12[0] * auxVec1[1] - r12[1] * auxVec1[0];

				auxVec2[0] = upVec2[1] * r12[2] - upVec2[2] * r12[1];
				auxVec2[1] = upVec2[2] * r12[0] - upVec2[0] * r12[2];
				auxVec2[2] = upVec2[0] * r12[1] - upVec2[1] * r12[0];

				projection2[0] = r12[1] * auxVec2[2] - r12[2] * auxVec2[1];
				projection2[1] = r12[2] * auxVec2[0] - r12[0] * auxVec2[2];
				projection2[2] = r12[0] * auxVec2[1] - r12[1] * auxVec2[0];

				// calculate cos(phi) and norm of projection vector
				for (unsigned i = 0; i < 3; i++) {
					cosPhi1 += r12[i] * upVec1[i];
					cosPhi2 -= r12[i] * upVec2[i];
					norm1 += projection1[i] * projection1[i];
					norm2 += projection2[i] * projection2[i];
				}

				norm1 = sqrt(norm1);
				norm2 = sqrt(norm2);

				// calculate cos(gamma) as dot product of projections
				for (unsigned i = 0; i < 3; i++) {
					projection1[i] /= norm1;
					projection2[i] /= norm2;
					cosGamma12 += projection1[i] * projection2[i];
				}

				// precaution to prevent numerically intractable values (e.g. VERY close to zero but not zero) and
				// values just barely out of boundaries -1/1, by rounding to 0,-1 or 1 respectively

				if (abs(cosPhi1) < roundingThreshold || abs(abs(cosPhi1) - 1) < roundingThreshold) {
					cosPhi1 = round(cosPhi1);
				}
				if (abs(cosPhi2) < roundingThreshold || abs(abs(cosPhi2) - 1) < roundingThreshold) {
					cosPhi2 = round(cosPhi2);
				}
				if (abs(cosGamma12) < roundingThreshold || abs(abs(cosGamma12) - 1) < roundingThreshold) {
					cosGamma12 = round(cosGamma12);
				}

				Gamma12 = acos(cosGamma12);
				// determine array element
				// NOTE: element 0 of array ODF is unused

				maximum = max(this->_phi1Increments, this->_phi2Increments);
				maximum = max(maximum, this->_gammaIncrements);

				for (unsigned i = 0; i < maximum; i++) {
					if (1. - i * 2. / (double)this->_phi1Increments >= cosPhi1 &&
						cosPhi1 > 1. - (i + 1) * 2. / (double)this->_phi1Increments) {
						ind_phi1 = i;
						bool1 = true;
					}

					if (1. - i * 2. / (double)this->_phi2Increments >= cosPhi2 &&
						cosPhi2 > 1. - (i + 1) * 2. / (double)this->_phi2Increments) {
						ind_phi2 = i;
						bool2 = true;
					}

					if (i * M_PI / (double)this->_gammaIncrements <= Gamma12 &&
						Gamma12 < (i + 1) * M_PI / (double)this->_gammaIncrements) {
						ind_gamma = i + 1;
						bool3 = true;
					}
				}

				if (ind_gamma == this->_gammaIncrements + 1) {
					ind_gamma = this->_gammaIncrements;
				}

				// manually assign bin for cos(...) == M_PI/-1, because loop only includes values < pi
				if (bool1 == 0 && cosPhi1 == -1.) {
					ind_phi1 = this->_phi1Increments - 1;
					bool1 = true;
				}

				if (bool2 == 0 && cosPhi2 == -1.) {
					ind_phi2 = this->_phi2Increments - 1;
					bool2 = true;
				}

				if (bool3 == 0 && Gamma12 == M_PI) {
					ind_gamma = this->_gammaIncrements;
					bool3 = true;
				}

				// notification if anything goes wrong during calculataion
				if (bool1 == 0 || bool2 == 0 || bool3 == 0) {
					global_log->warning() << "Array element in ODF calculation not properly assigned!" << endl;
					global_log->warning() << "Mol-ID 1 = " << mol1.getID() << "  Mol-ID 2 = " << it->getID() << endl;
					global_log->warning()
						<< "upVec1=" << upVec1[0] << " " << upVec1[1] << " " << upVec1[2] << " " << endl;
					global_log->warning()
						<< "upVec2=" << upVec2[0] << " " << upVec2[1] << " " << upVec2[2] << " " << endl;
					global_log->warning() << "r12=" << r12[0] << " " << r12[1] << " " << r12[2] << " " << endl;
					global_log->warning() << "[cosphi1 cosphi2 cosgamma12] = [" << cosPhi1 << " " << cosPhi2 << " "
										  << cosGamma12 << "]" << endl;
					global_log->warning() << "indices are " << ind_phi1 << " " << ind_phi2 << " " << ind_gamma << endl;
				}

				elementID = ind_phi1 * this->_phi2Increments * this->_gammaIncrements +
							(ind_phi2 * this->_gammaIncrements) + ind_gamma;

				// determine component pairing

				if (cid == 0 && it->getComponentLookUpID() == 0) {
					this->_localODF11[elementID]++;
				}

				else if (cid == 0 && it->getComponentLookUpID() == 1) {
					this->_localODF12[elementID]++;
				}

				else if (cid == 1 && it->getComponentLookUpID() == 1) {
					this->_localODF22[elementID]++;
				}

				else {
					this->_localODF21[elementID]++;
				}
			}
		}
	}
}

void ODF::collect(DomainDecompBase* domainDecomp) {
	if (this->_numPairs == 1) {
		domainDecomp->collCommInit(this->_numElements);

		for (unsigned long i = 0; i < this->_numElements; i++) {
			domainDecomp->collCommAppendUnsLong(this->_localODF11[i]);
		}
		domainDecomp->collCommAllreduceSum();

		for (unsigned long i = 0; i < this->_numElements; i++) {
			this->_ODF11[i] = domainDecomp->collCommGetUnsLong();
		}
		domainDecomp->collCommFinalize();
	}

	else {
		domainDecomp->collCommInit(this->_numElements * 4);

		for (unsigned long i = 0; i < this->_numElements; i++) {
			domainDecomp->collCommAppendUnsLong(this->_localODF11[i]);
			domainDecomp->collCommAppendUnsLong(this->_localODF12[i]);
			domainDecomp->collCommAppendUnsLong(this->_localODF22[i]);
			domainDecomp->collCommAppendUnsLong(this->_localODF21[i]);
		}
		domainDecomp->collCommAllreduceSum();

		for (unsigned long i = 0; i < this->_numElements; i++) {
			this->_ODF11[i] = domainDecomp->collCommGetUnsLong();
			this->_ODF12[i] = domainDecomp->collCommGetUnsLong();
			this->_ODF22[i] = domainDecomp->collCommGetUnsLong();
			this->_ODF21[i] = domainDecomp->collCommGetUnsLong();
		}
		domainDecomp->collCommFinalize();
	}
}

void ODF::output(Domain* /*domain*/, long unsigned timestep) {
	global_log->info() << "[ODF] writing output" << std::endl;
	// Setup outfile

	double cosPhi1 = 1.;
	double cosPhi2 = 1. - 2. / this->_phi2Increments;
	double Gamma12 = 0.;
	string prefix;
	ostringstream osstrm;
	osstrm << this->_outputPrefix;
	osstrm.fill('0');
	osstrm.width(7);
	osstrm << right << timestep;
	prefix = osstrm.str();
	osstrm.str("");
	osstrm.clear();

	if (this->_numPairs == 1) {
		string ODF11name = prefix + ".ODF11";
		ofstream outfile(ODF11name.c_str());
		outfile.precision(6);

		outfile << "//Output generated by ODF plugin\n"
				<< "//Angular distribution at time step = " << timestep << " for component pairing pairing 11\n";
		outfile << "cosPhi1\tcosPhi2\tGamma12\tcount\n";
		for (unsigned long i = 0; i < this->_numElements - 1; i++) {
			Gamma12 += M_PI / (double)this->_gammaIncrements;
			if (i % this->_gammaIncrements == 0) {
				cosPhi2 -= 2. / (double)this->_phi2Increments;
				Gamma12 = M_PI / (double)this->_gammaIncrements;
			}
			if (i % (this->_gammaIncrements * this->_phi2Increments) == 0) {
				cosPhi1 -= 2. / (double)this->_phi1Increments;
				cosPhi2 = 1. - 2. / (double)this->_phi2Increments;
			}
			outfile << cosPhi1 << "\t" << cosPhi2 << "\t" << Gamma12 << "\t" << this->_ODF11[i + 1] << "\n";
		}
		outfile.close();
	} else {
		string ODF11name = prefix + ".ODF11";
		string ODF12name = prefix + ".ODF12";
		string ODF22name = prefix + ".ODF22";
		string ODF21name = prefix + ".ODF21";

		ofstream ODF11(ODF11name.c_str());
		ofstream ODF12(ODF12name.c_str());
		ofstream ODF22(ODF22name.c_str());
		ofstream ODF21(ODF21name.c_str());
		ODF11.precision(5);
		ODF12.precision(5);
		ODF22.precision(5);
		ODF21.precision(5);

		ODF11 << "//Output generated by ODF plugin\n"
			  << "//Angular distribution at time step = " << timestep << " for component pairing pairing 11\n";
		ODF11 << "cosPhi1\tcosPhi2\tcosGamma12\tcount\n";
		ODF12 << "//Output generated by ODF plugin\n"
			  << "//Angular distribution at time step = " << timestep << " for component pairing pairing 12\n";
		ODF12 << "cosPhi1\tcosPhi2\tcosGamma12\tcount\n";
		ODF22 << "//Output generated by ODF plugin\n"
			  << "//Angular distribution at time step = " << timestep << " for component pairing pairing 22\n";
		ODF22 << "phi1\tphi2\tgamm12\tcount\n";
		ODF21 << "//Output generated by ODF plugin\n"
			  << "//Angular distribution at time step = " << timestep << " for component pairing pairing 21\n";
		ODF21 << "cosPhi1\tcosPhi2\tcosGamma12\tcount\n";

		for (unsigned long i = 0; i < this->_numElements - 1; i++) {
			Gamma12 += M_PI / (double)this->_gammaIncrements;
			if (i % this->_gammaIncrements == 0) {
				cosPhi2 -= 2. / (double)this->_phi2Increments;
				Gamma12 = M_PI / (double)this->_gammaIncrements;
			}
			if (i % (this->_gammaIncrements * this->_phi2Increments) == 0) {
				cosPhi1 -= 2. / (double)this->_phi1Increments;
				cosPhi2 = 1. - 2. / (double)this->_phi2Increments;
			}

			ODF11 << cosPhi1 << "\t" << cosPhi2 << "\t" << Gamma12 << "\t" << this->_ODF11[i + 1] << "\n";
			ODF12 << cosPhi1 << "\t" << cosPhi2 << "\t" << Gamma12 << "\t" << this->_ODF12[i + 1] << "\n";
			ODF22 << cosPhi1 << "\t" << cosPhi2 << "\t" << Gamma12 << "\t" << this->_ODF22[i + 1] << "\n";
			ODF21 << cosPhi1 << "\t" << cosPhi2 << "\t" << Gamma12 << "\t" << this->_ODF21[i + 1] << "\n";
		}
		ODF11.close();
		ODF12.close();
		ODF22.close();
		ODF21.close();
	}
}
