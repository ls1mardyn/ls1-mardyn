
#include <cmath>

#include "Domain.h"
#include "Simulation.h"
#include "longRange/Homogeneous.h"
#include "molecules/Molecule.h"
#include "particleContainer/ParticleContainer.h"

#include "utils/Logger.h"
#include "utils/mardyn_assert.h"


Homogeneous::Homogeneous(double cutoffRadius, double cutoffRadiusLJ, Domain* domain, ParticleContainer* particleContainer, Simulation* simulation) :
	_components(simulation->getEnsemble()->getComponents()),
	_domain(domain),
	_particleContainer(particleContainer),
	_cutoff(cutoffRadius),
	_cutoffLJ(cutoffRadiusLJ)
{}

void Homogeneous::init() {
	_comp2params = _domain->getComp2Params();
	Log::global_log->info() << "Long range correction for homogeneous systems is used " << std::endl;
	double upotCorrLJ = 0.;
	double virialCorrLJ = 0.;
	double mySelbstTerm = 0.;

	unsigned int numcomp = _components->size();
	unsigned long nummolecules = 0;
	for (unsigned int i = 0; i < numcomp; ++i) {
		Component& ci = (*_components)[i];
		nummolecules += ci.getNumMolecules();
		unsigned int numljcentersi = ci.numLJcenters();
		unsigned int numchargesi = ci.numCharges();
		unsigned int numdipolesi = ci.numDipoles();

		// effective dipoles computed from point charge distributions
		double chargeBalance[3];
		for (unsigned d = 0; d < 3; d++) {
			chargeBalance[d] = 0;
		}
		for (unsigned int si = 0; si < numchargesi; si++) {
			double tq = ci.charge(si).q();
			for (unsigned d = 0; d < 3; d++) {
				chargeBalance[d] += tq * ci.charge(si).r()[d];
			}
		}
		// point dipoles
		for (unsigned int si = 0; si < numdipolesi; ++si) {
			const double tmy = ci.dipole(si).absMy();
			double evect = 0;
			for (unsigned d = 0; d < 3; d++) {
				evect += ci.dipole(si).e()[d] * ci.dipole(si).e()[d];
			}
			double norm = 1.0 / sqrt(evect);
			for (unsigned d = 0; d < 3; d++) {
				chargeBalance[d] += tmy * ci.dipole(si).e()[d] * norm;
			}
		}
		double my2 = 0.0;
		for (unsigned d = 0; d < 3; d++) {
			my2 += chargeBalance[d] * chargeBalance[d];
		}
		mySelbstTerm += my2 * ci.getNumMolecules();

		for (unsigned int j = 0; j < numcomp; ++j) {
			Component& cj = (*_components)[j];
			unsigned int numljcentersj = cj.numLJcenters();
			ParaStrm& params = _comp2params(i, j);
			params.reset_read();
			// LJ centers
			for (unsigned int si = 0; si < numljcentersi; ++si) {
				const double xi = ci.ljcenter(si).rx();
				const double yi = ci.ljcenter(si).ry();
				const double zi = ci.ljcenter(si).rz();
				const double tau1 = sqrt(xi * xi + yi * yi + zi * zi);
				for (unsigned int sj = 0; sj < numljcentersj; ++sj) {
					double xj = cj.ljcenter(sj).rx();
					double yj = cj.ljcenter(sj).ry();
					double zj = cj.ljcenter(sj).rz();
					double tau2 = sqrt(xj * xj + yj * yj + zj * zj);
					if (tau1 + tau2 >= _cutoffLJ) {
						std::ostringstream error_message;
						error_message << "Error calculating cutoff corrections, rc too small" << std::endl;
						MARDYN_EXIT(error_message.str());
					}
					double eps24;
					params >> eps24;
					double sig2;
					params >> sig2;
					double uLJshift6;
					params >> uLJshift6;  // 0 unless TRUNCATED_SHIFTED
					if (uLJshift6 == 0.0) {
						double fac = double(ci.getNumMolecules()) * double(cj.getNumMolecules()) * eps24;
						if (tau1 == 0. && tau2 == 0.) {
							upotCorrLJ += fac * (this->_TICCu(-6, _cutoffLJ, sig2) - this->_TICCu(-3, _cutoffLJ, sig2));
							virialCorrLJ +=
								fac * (this->_TICCv(-6, _cutoffLJ, sig2) - this->_TICCv(-3, _cutoffLJ, sig2));
						} else if (tau1 != 0. && tau2 != 0.) {
							upotCorrLJ += fac * (this->_TISSu(-6, _cutoffLJ, sig2, tau1, tau2) -
												 this->_TISSu(-3, _cutoffLJ, sig2, tau1, tau2));
							virialCorrLJ += fac * (this->_TISSv(-6, _cutoffLJ, sig2, tau1, tau2) -
												   this->_TISSv(-3, _cutoffLJ, sig2, tau1, tau2));
						} else {
							if (tau2 == 0.) {
								tau2 = tau1;
							}
							upotCorrLJ += fac * (this->_TICSu(-6, _cutoffLJ, sig2, tau2) -
												 this->_TICSu(-3, _cutoffLJ, sig2, tau2));
							virialCorrLJ += fac * (this->_TICSv(-6, _cutoffLJ, sig2, tau2) -
												   this->_TICSv(-3, _cutoffLJ, sig2, tau2));
						}
					}
				}
			}
		}
	}
	_upotCorrLJ_no_num_molecules = upotCorrLJ;
	_virialCorrLJ_no_num_molecules = virialCorrLJ;
	_mySelbstTerm_no_num_molecules = mySelbstTerm;
}

void Homogeneous::calculateLongRange() {
	const double globalRho = _domain->getglobalRho();
	const unsigned long globalNumMolecules = _domain->getglobalNumMolecules();
	const double epsilonRF = _domain->getepsilonRF();
	const double fac = M_PI * globalRho / (3. * static_cast<double>(globalNumMolecules));
	const double upotCorrLJ = fac * _upotCorrLJ_no_num_molecules;
	const double virialCorrLJ = -fac * _virialCorrLJ_no_num_molecules;

	const double epsRFInvrc3 = 2. * (epsilonRF - 1.) / ((_cutoff * _cutoff * _cutoff) * (2. * epsilonRF + 1.));
	const double mySelbstTerm = -0.5 * epsRFInvrc3 * _mySelbstTerm_no_num_molecules;

	const double upotCorr = upotCorrLJ + mySelbstTerm;
	const double virialCorr = virialCorrLJ + 3. * mySelbstTerm;

	Log::global_log->debug() << "Far field terms: U_pot_correction = " << upotCorr << " virial_correction = " << virialCorr
					   << std::endl;
	_domain->setUpotCorr(upotCorr);
	_domain->setVirialCorr(virialCorr);
	
	// Pot. energy and virial are also stored in molecule class
	// To maintain continuity (sum of U over N == U_global_domain), the correction values are added to each molecule
	// For the virial, the correction is distributed equally to the diagonal elements
	// NOTE: UPotCorr may need some sophistication as correction could depend on component
	_uPotCorrPerMol = upotCorr/globalNumMolecules;
	const double virialCorrPerMol[9] = {
		virialCorr/(3.*globalNumMolecules),  // rxfx
		virialCorr/(3.*globalNumMolecules),  // ryfy
		virialCorr/(3.*globalNumMolecules),  // rzfz
		0., 								 // rxfy
		0., 								 // rxfz
		0., 								 // ryfz
		0., 								 // ryfx
		0.,  								 // rzfx
		0., 								 // rzfy
	};
	for (auto tempMol = _particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMol.isValid(); ++tempMol) {
		tempMol->Viadd(virialCorrPerMol);
		tempMol->Uadd(_uPotCorrPerMol);
	}
}

double Homogeneous::getUpotCorr(Molecule* /* mol */) {
	// unsigned int cid = mol->componentid();
	// TODO: May need some sophistication as correction could depend on component
	return _uPotCorrPerMol;
}

double Homogeneous::_TICCu(int n, double rc, double sigma2) {
	return -pow(rc, 2 * n + 3) / (pow(sigma2, n) * (2 * n + 3));
}

double Homogeneous::_TICSu(int n, double rc, double sigma2, double tau) {
	return -(pow(rc + tau, 2 * n + 3) - pow(rc - tau, 2 * n + 3)) * rc /
			   (4 * pow(sigma2, n) * tau * (n + 1) * (2 * n + 3)) +
		   (pow(rc + tau, 2 * n + 4) - pow(rc - tau, 2 * n + 4)) /
			   (4 * pow(sigma2, n) * tau * (n + 1) * (2 * n + 3) * (2 * n + 4));
}

double Homogeneous::_TISSu(int n, double rc, double sigma2, double tau1, double tau2) {
	double tauMinus, tauPlus;
	tauPlus = tau1 + tau2;
	tauMinus = tau1 - tau2;
	return -(pow(rc + tauPlus, 2 * n + 4) - pow(rc + tauMinus, 2 * n + 4) - pow(rc - tauMinus, 2 * n + 4) +
			 pow(rc - tauPlus, 2 * n + 4)) *
			   rc / (8 * pow(sigma2, n) * tau1 * tau2 * (n + 1) * (2 * n + 3) * (2 * n + 4)) +
		   (pow(rc + tauPlus, 2 * n + 5) - pow(rc + tauMinus, 2 * n + 5) - pow(rc - tauMinus, 2 * n + 5) +
			pow(rc - tauPlus, 2 * n + 5)) /
			   (8 * pow(sigma2, n) * tau1 * tau2 * (n + 1) * (2 * n + 3) * (2 * n + 4) * (2 * n + 5));
}

double Homogeneous::_TICCv(int n, double rc, double sigma2) { return 2 * n * _TICCu(n, rc, sigma2); }

double Homogeneous::_TICSv(int n, double rc, double sigma2, double tau) {
	return -(pow(rc + tau, 2 * n + 2) - pow(rc - tau, 2 * n + 2)) * rc * rc / (4 * pow(sigma2, n) * tau * (n + 1)) -
		   3 * _TICSu(n, rc, sigma2, tau);
}

double Homogeneous::_TISSv(int n, double rc, double sigma2, double tau1, double tau2) {
	double tauMinus, tauPlus;
	tauPlus = tau1 + tau2;
	tauMinus = tau1 - tau2;
	return -(pow(rc + tauPlus, 2 * n + 3) - pow(rc + tauMinus, 2 * n + 3) - pow(rc - tauMinus, 2 * n + 3) +
			 pow(rc - tauPlus, 2 * n + 3)) *
			   rc * rc / (8 * pow(sigma2, n) * tau1 * tau2 * (n + 1) * (2 * n + 3)) -
		   3 * _TISSu(n, rc, sigma2, tau1, tau2);
}
