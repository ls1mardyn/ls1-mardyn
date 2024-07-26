#include "molecules/Comp2Param.h"

#include <cmath>
#include <functional>

#include "utils/Logger.h"
#include "Simulation.h"


void Comp2Param::initialize(
		const std::vector<Component>& components, std::map<int,std::map<int,MixingRuleBase*>> mixcoeff,
		double epsRF, double rc, double rcLJ)
{
	m_numcomp = components.size();
	m_ssparatbl.redim(m_numcomp, m_numcomp);

	// interaction between LJ centers
	for (unsigned int compi = 0; compi < m_numcomp; ++compi) {
		ParaStrm& pstrmii = m_ssparatbl(compi, compi);
		unsigned int nci = components[compi].numLJcenters();
		double epsi, sigi, epsj, sigj, epsilon24, sigma2, shift6i;
		// single-component interaction
		for (unsigned int centeri = 0; centeri < nci; ++centeri) {
			const LJcenter& ljcenteri = static_cast<const LJcenter&>(components[compi].ljcenter(centeri));
			epsi = ljcenteri.eps();
			sigi = ljcenteri.sigma();
			shift6i = ljcenteri.shift6();
			for (unsigned int centerj = 0; centerj < nci; ++centerj) {
				const LJcenter& ljcenterj = static_cast<const LJcenter&>(components[compi].ljcenter(centerj));
				epsj = ljcenterj.eps();
				sigj = ljcenterj.sigma();
				epsilon24 = 24. * sqrt(epsi * epsj);
				sigma2 = .5 * (sigi + sigj);
				sigma2 *= sigma2;
				pstrmii << epsilon24;
				pstrmii << sigma2;
				pstrmii << shift6i;
			}
		}
		// interaction between different components
		for (unsigned int compj = compi + 1; compj < m_numcomp; ++compj) {
			ParaStrm& pstrmij = m_ssparatbl(compi, compj);
			unsigned int ncj = components[compj].numLJcenters();
			const auto mixingrule = mixcoeff[compi][compj];
			// Generic mixing functions
			std::function<double(double, double)> mixingSigma;
			std::function<double(double, double)> mixingEpsilon;
			// Get parameters
			if (mixingrule->getType() == "LB") {
				const double eta = mixingrule->getParameters().at(0);
				const double xi  = mixingrule->getParameters().at(1);
				mixingSigma = [=](double epsi, double epsj) { return xi * sqrt(epsi * epsj); };
				mixingEpsilon = [=](double epsi, double epsj) { return eta * (sigi + sigj); };
#ifndef NDEBUG
				Log::global_log->debug() << "Mixing : cid+1(compi)=" << compi+1 << " <--> cid+1(compj)=" << compj+1 << ": xi=" << xi << ", eta=" << eta << std::endl;
#endif
			} else {
				Log::global_log->error() << "Mixing: Only LB rule supported" << std::endl;
				Simulation::exit(1);
			}
			double shift6combined, sigperrc2, sigperrc6;
			for (unsigned int centeri = 0; centeri < nci; ++centeri) {
				const LJcenter& ljcenteri = static_cast<const LJcenter&>(components[compi].ljcenter(centeri));
				epsi = ljcenteri.eps();
				sigi = ljcenteri.sigma();
				for (unsigned int centerj = 0; centerj < ncj; ++centerj) {
					const LJcenter& ljcenterj = static_cast<const LJcenter&>(components[compj].ljcenter(centerj));
					epsj = ljcenterj.eps();
					sigj = ljcenterj.sigma();
					epsilon24 = 24. * mixingEpsilon(epsi, epsj);
					sigma2 = 0.5 * mixingSigma(sigi, sigj);
					sigma2 *= sigma2;
					sigperrc2 = sigma2 / (rcLJ * rcLJ);
					sigperrc6 = sigperrc2 * sigperrc2 * sigperrc2;
					shift6combined = epsilon24 * (sigperrc6 - sigperrc6 * sigperrc6);
					pstrmij << epsilon24;
					pstrmij << sigma2;
					pstrmij << shift6combined;
#ifndef NDEBUG
					Log::global_log->debug() << "Component " << compi << ": eps24=" << epsilon24 << " sig2=" << sigma2 << " shift6=" << shift6combined << std::endl;
#endif
				}
			}
			ParaStrm& pstrmji = m_ssparatbl(compj, compi);
			for (unsigned int centerj = 0; centerj < ncj; ++centerj) {
				const LJcenter& ljcenterj = static_cast<const LJcenter&>(components[compj].ljcenter(centerj));
				epsj = ljcenterj.eps();
				sigj = ljcenterj.sigma();
				for (unsigned int centeri = 0; centeri < nci; ++centeri) {
					const LJcenter& ljcenteri = static_cast<const LJcenter&>(components[compi].ljcenter(centeri));
					epsi = ljcenteri.eps();
					sigi = ljcenteri.sigma();
					epsilon24 = 24. * mixingEpsilon(epsi, epsj);
					sigma2 = 0.5 * mixingSigma(sigi, sigj);
					sigma2 *= sigma2;
					sigperrc2 = sigma2 / (rcLJ * rcLJ);
					sigperrc6 = sigperrc2 * sigperrc2 * sigperrc2;
					shift6combined = epsilon24 * (sigperrc6 - sigperrc6 * sigperrc6);
					pstrmji << epsilon24;
					pstrmji << sigma2;
					pstrmji << shift6combined;
				}
			}
		}
	}

	for (unsigned int compi = 0; compi < m_numcomp; ++compi) {
		unsigned nei = components[compi].numCharges();
		unsigned int nqi = components[compi].numQuadrupoles();
		unsigned int ndi = components[compi].numDipoles();
		for (unsigned int compj = 0; compj < m_numcomp; ++compj) {
			ParaStrm& pstrmij = m_ssparatbl(compi, compj);
			unsigned nej = components[compj].numCharges();
			unsigned int nqj = components[compj].numQuadrupoles();
			unsigned int ndj = components[compj].numDipoles();
			for (unsigned ei = 0; ei < nei; ei++) {
				const Charge& chargei = static_cast<const Charge&> (components[compi].charge(ei));
				double cvali = chargei.q();
				// Charge-Charge
				for (unsigned ej = 0; ej < nej; ej++) {
					const Charge& chargej = static_cast<const Charge&> (components[compj].charge(ej));
					double cvalj = chargej.q();
					double q1q2per4pie0 = cvali * cvalj; // 4pi*epsilon0 = Einheit!!!
					pstrmij << q1q2per4pie0;
				}
				// Charge-Quadrupole
				for (unsigned qj = 0; qj < nqj; qj++) {
					const Quadrupole& quadrupolej = static_cast<const Quadrupole&> (components[compj].quadrupole(qj));
					double absqj = quadrupolej.absQ();
					double qQ05per4pie0 = 0.5 * cvali * absqj;
					pstrmij << qQ05per4pie0;
				}
				// Charge-Dipole
				for (unsigned int dj = 0; dj < ndj; dj++) {
					const Dipole& dipolej = static_cast<const Dipole&> (components[compj].dipole(dj));
					double absmyj = dipolej.absMy();
					double minusqmyper4pie0 = -cvali * absmyj;
					pstrmij << minusqmyper4pie0;
				}
			}
			for (unsigned int qi = 0; qi < nqi; ++qi) {
				const Quadrupole& quadrupolei = static_cast<const Quadrupole&> (components[compi].quadrupole(qi));
				double absqi = quadrupolei.absQ();
				// Quadrupole-Charge
				for (unsigned ej = 0; ej < nej; ej++) {
					const Charge& chargej = static_cast<const Charge&> (components[compj].charge(ej));
					double cvalj = chargej.q();
					double qQ05per4pie0 = 0.5 * cvalj * absqi;
					pstrmij << qQ05per4pie0;
				}
				// Quadrupole-Quadrupole
				for (unsigned int qj = 0; qj < nqj; ++qj) {
					const Quadrupole& quadrupolej = static_cast<const Quadrupole&> (components[compj].quadrupole(qj));
					double absqj = quadrupolej.absQ();
					double q2075 = .75 * absqi * absqj;
					pstrmij << q2075;
				}
				// Quadrupole-Dipole
				for (unsigned int dj = 0; dj < ndj; ++dj) {
					const Dipole& dipolej = static_cast<const Dipole&> (components[compj].dipole(dj));
					double absmyj = dipolej.absMy();
					double qmy15 = 1.5 * absqi * absmyj;
					pstrmij << qmy15;
				}
			}
			for (unsigned int di = 0; di < ndi; ++di) {
				const Dipole& dipolei = static_cast<const Dipole&> (components[compi].dipole(di));
				double absmyi = dipolei.absMy();
				double epsRFInvrc3 = 2. * (epsRF - 1.) / ((rc * rc * rc) * (2. * epsRF + 1.));
				// Dipole-Charge
				for (unsigned ej = 0; ej < nej; ej++) {
					const Charge& chargej = static_cast<const Charge&> (components[compj].charge(ej));
					double cvalj = chargej.q();
					double minusqmyper4pie0 = -cvalj * absmyi;
					pstrmij << minusqmyper4pie0;
				}
				// Dipole-Quadrupole
				for (unsigned int qj = 0; qj < nqj; ++qj) {
					const Quadrupole& quadrupolej = static_cast<const Quadrupole&> (components[compj].quadrupole(qj));
					double absqj = quadrupolej.absQ();
					double myq15 = 1.5 * absmyi * absqj;
					pstrmij << myq15;
				}
				// Dipole-Dipole
				for (unsigned int dj = 0; dj < ndj; ++dj) {
					const Dipole& dipolej = static_cast<const Dipole&> (components[compj].dipole(dj));
					double absmyj = dipolej.absMy();
					double my2 = absmyi * absmyj;
					pstrmij << my2;
					double rffac = my2 * epsRFInvrc3;
					pstrmij << rffac;
				}
			}
		}
	}
}
