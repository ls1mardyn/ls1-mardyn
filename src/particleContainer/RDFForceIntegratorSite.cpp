/*
 * RDFForceIntegratorSite.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: tijana
 */

#include "RDFForceIntegratorSite.h"

RDFForceIntegratorSite::RDFForceIntegratorSite(
		ParticleContainer* moleculeContainer, double rc, double d,
		std::vector<std::vector<double> >* globalADist,
		std::vector<std::vector<std::vector<double> > >* globalSiteADist) :
	RDFForceIntegrator(moleculeContainer, rc, d, globalADist, globalSiteADist) {
	// TODO Auto-generated constructor stub
	_rho = 0;
}

RDFForceIntegratorSite::~RDFForceIntegratorSite() {
	// TODO Auto-generated destructor stub
}

double RDFForceIntegratorSite::integrateRDFSite(Molecule* mol,
		double* normal_dim, int* boundary, int plane, unsigned int site,
		double* force, bool add_influence) {

	std::vector<double> globalAcc = (*_globalADist)[0];
	std::vector<std::vector<double> > globalSiteAcc = (*_globalSiteADist)[0];

	double V = (_rmax[0] - _rmin[0]) * (_rmax[1] - _rmin[1]) * (_rmax[2]
			- _rmin[2]);

	double potential = 0;

	// number density of the domain
	//int numMolecules = _moleculeContainer->getNumberOfParticles();

	//double rho = _numMolecules / (V);

	int bin; // bin for the radius that rdf is read for
	double g;
	for (int i = 0; i < 3; i++)
		force[i] = 0;
	for (double z = normal_dim[0] + _d / 2; z <= normal_dim[1]; z += _d) {
		for (double x = _d / 2; x <= std::sqrt( normal_dim[1] *  normal_dim[1] - z * z); x += _d) {
			double r = std::sqrt(x * x + z * z);
			if (mol->numSites() == 1) {

				bin = (int) (r * globalAcc.size() / _rc - 0.5);
				g = globalAcc[bin];
				double sig2 = mol->getSigma() * mol->getSigma();
				double r2 = r * r;
				double lj6 = sig2 * sig2 * sig2 / (r2 * r2 * r2);
				double lj12_6 = 6 * 4 * mol->getEps() * (lj6 * lj6 - lj6); //times 6 like potforce returns
				double ljf = 24 * mol->getEps() * (lj6 - 2 * lj6 * lj6) / r;

				double f[3] = { 0, 0, 0 };
				f[plane] = boundary[plane] * 2 * PI * _rho * g * ljf * z * x
						* _d * _d / r;
				force[plane] += f[plane];
				potential += boundary[plane] * 2 * PI * _rho * g * lj12_6 * x
						* _d * _d;
				if (add_influence)
					mol->Fljcenteradd(site, f);

				if (boundary[0] == -1 && plane == 0) {
					mol->addLeftxRdfInfluence(site, f);
				}
			} else {

				bin = (int) (r * globalSiteAcc[0].size() / (_rc + 2
						* mol->ljcenter_disp(site)) - 0.5);

				// if multiple LJ centers, use site-site rdf
				// iterate through sites, treat cell center as a site

				for (unsigned int j = 0; j < mol->numLJcenters(); j++) {

					g = globalSiteAcc[site * mol->numLJcenters() + j][bin];

					if (site != j)
						g /= 2;

					double sig2 = mol->getSigma() * mol->getSigma();
					double r2 = r * r;
					double lj6 = sig2 * sig2 * sig2 / (r2 * r2 * r2);

					double ljf = 24 * mol->getEps() * (lj6 - 2 * lj6 * lj6) / r;
					double lj12_6 = 6 * 4 * mol->getEps() * (lj6 * lj6 - lj6); //times 6 like potforce returns

					double f[3] = { 0, 0, 0 };
					f[plane] = boundary[plane] * 2 * PI * _rho * g * ljf * z * x
							* _d * _d / r;
					force[plane] += f[plane];
					potential += boundary[plane] * 2 * PI * _rho * g * lj12_6
							* x * _d * _d;
					if (add_influence)
						mol->Fljcenteradd(site, f);


					if (boundary[0] == -1 && plane == 0) {
						mol->addLeftxRdfInfluence(site, f);

					}
				}

			}
		}
	}
	return potential;
}

double RDFForceIntegratorSite::processMolecule(Molecule* currentMolecule,
		double* force, bool add_influence, bool unit_test) {
	double pot = 0;
	//std::cout<<"here 1"<<std::endl;
	for (unsigned int site = 0; site < currentMolecule->numSites(); site++) {
		double rm[3] = { currentMolecule->r(0), currentMolecule->r(1),
				currentMolecule->r(2) };

		double r[3] = { currentMolecule->r(0)
				+ currentMolecule->site_d(site)[0], currentMolecule->r(1)
				+ currentMolecule->site_d(site)[1], currentMolecule->r(2)
				+ currentMolecule->site_d(site)[2] };

		// if this is a halo molecule, skip it
		if (rm[0] < _rmin[0] || rm[1] < _rmin[1] || rm[2] < _rmin[2] || rm[0]
				> _rmax[0] || rm[1] > _rmax[1] || rm[2] > _rmax[2])
			return 0;

		// if the molecule is not within rc of the boundary, skip it
		if (rm[0] > _low_limit[0] && rm[0] < _high_limit[0] && rm[1]
				> _low_limit[1] && rm[1] < _high_limit[1] && rm[2]
				> _low_limit[2] && rm[2] < _high_limit[2])
			return 0;

		if (currentMolecule->numSites() != currentMolecule->numLJcenters()) {
			std::cout
					<< "Molecule consists of something other than LJ centers. In this case RDF cannot be used.";
			return 0;
		}

		int boundary[3] = { 0, 0, 0 };
		if (rm[0] < _low_limit[0])
			boundary[0] = -1;
		else if (rm[0] >= _high_limit[0])
			boundary[0] = 1;
		if (rm[1] < _low_limit[1])
			boundary[1] = -1;
		else if (rm[1] >= _high_limit[1])
			boundary[1] = 1;
		if (rm[2] < _low_limit[2])
			boundary[2] = -1;
		else if (rm[2] >= _high_limit[2])
			boundary[2] = 1;

		//if (boundary[0] == 0 && boundary[1] == 0 && boundary[2] == 0)

		if (boundary[0] != -1 && boundary[0] != 1)
			return 0;

		// if molecule close to the boundary, add RDF force
		// integration limits for axes

		double normal_lim[2] = { 0, _rc };

		// box size for the numerical quadrature
		//_dn = _dr = currentMolecule->getSigma() / 20;

		if (boundary[0] == -1) {
			normal_lim[0] = std::abs(r[0] - _rmin[0]);
			if (r[0] < _rmin[0])
				normal_lim[1] -= _rmin[0] - r[0];
			//			if (r[0] < _rmin[0])
			//				normal_lim[0] = std::abs(2 * r[0] - _rmin[0]);

		} else if (boundary[0] == 1) {
			normal_lim[0] = std::abs(_rmax[0] - r[0]);
			if (r[0] > _rmax[0])
				normal_lim[1] -= r[0] - _rmax[0];
			//			if (r[0] > _rmax[0])
			//				normal_lim[0] = std::abs(_rmax[0] - 2 * r[0]);
		}
		//std::cout<<"here2"<<std::endl;

		pot += integrateRDFSite(currentMolecule, normal_lim, boundary, 0, site,
				force, add_influence);
		//if (std::abs(force[0]) < 0.00001) std::cout<<"rm: "<<rm[0]<<" r: "<<r[0]<<std::endl;
//
//		if (r[0] < 3.23 && r[0] > 3.21)
//			std::cout<<"r: "<<r[0]<<" f: "<<force[0]<<std::endl;
//		if (r[0] < 3.21 && r[0] > 3.2)
//			std::cout<<"r: "<<r[0]<<" f: "<<force[0]<<std::endl;
		// xy plane is the boundary
		// integrate xy
		/*
		 for (unsigned int i = 0; boundary[2] != 0 && i
		 < componentIds.size(); i++) {
		 if (boundary[2] == -1) {
		 zlim[0] = r2 - rc;
		 zlim[1] = (r2 > rmin[2] ? rmin[2] : 2 * r2 - rmin[2]);
		 normal_lim[0] = abs(r2 - rmin[2]);
		 } else if (boundary[2] == 1) {
		 zlim[0] = (r2 < rmax[2] ? rmax[2] : 2 * r2 - rmax[2]);
		 zlim[1] = r2 + rc;
		 normal_lim[0] = abs(rmax[2] - r2);
		 }
		 double h = zlim[1] - zlim[0];
		 double a = std::sqrt(h * (2 * rc - h));
		 if (boundary[0] == 0)
		 integrateRDFSite(normal_lim, currentMolecule, rc, dn, dr,
		 globalADist[currentMolecule->componentid()
		 * totalComponents + i],
		 globalSiteADist[currentMolecule->componentid()
		 * totalComponents + i], 2, currSum, site,
		 boundary);
		 else {
		 if (boundary[0] == -1) {
		 xlim[0] = rmin[0];
		 xlim[1] = r0 + a;
		 } else if (boundary[0] == 1) {
		 xlim[0] = r0 - a;
		 xlim[1] = rmax[0];
		 }

		 ylim[0] = r1 - a;
		 ylim[1] = r1 + a;
		 integrateRDFCartesian(xlim, ylim, zlim, currentMolecule,
		 rc, dx, dy, dz,
		 globalADist[currentMolecule->componentid()
		 * totalComponents + i],
		 globalSiteADist[currentMolecule->componentid()
		 * totalComponents + i], 2, currSum, site,
		 boundary);
		 }

		 }

		 // xz plane is the boundary
		 // integrate xz plane
		 for (unsigned int i = 0; boundary[1] != 0 && i
		 < componentIds.size(); i++) {
		 if (boundary[1] == -1) {
		 ylim[1] = (r1 > rmin[1] ? rmin[1] : 2 * r1 - rmin[1]);
		 ylim[0] = r1 - rc;
		 normal_lim[0] = abs(r1 - rmin[1]);
		 } else if (boundary[1] == 1) {
		 ylim[0] = (r1 < rmax[1] ? rmax[1] : 2 * r1 - rmax[1]);
		 ylim[1] = r1 + rc;
		 normal_lim[0] = abs(rmax[1] - r1);
		 }
		 double h = ylim[1] - ylim[0];
		 double a = std::sqrt(h * (2 * rc - h));
		 if (boundary[0] == -1) {
		 xlim[0] = rmin[0];
		 xlim[1] = r0 + a;
		 } else if (boundary[0] == 1) {
		 xlim[0] = r0 - a;
		 xlim[1] = rmax[0];
		 } else {
		 xlim[0] = r0 - a;
		 xlim[1] = r0 + a;
		 }
		 if (boundary[2] == -1) {
		 zlim[0] = rmin[2];
		 zlim[1] = r2 + rc;
		 } else if (boundary[2] == 1) {
		 zlim[0] = r2 - rc;
		 zlim[1] = rmax[2];
		 } else {
		 zlim[0] = r2 - a;
		 zlim[1] = r2 + a;
		 }
		 if (boundary[0] == 0 && boundary[2] == 0)
		 integrateRDFSite(normal_lim, currentMolecule, rc, dr, dn,
		 globalADist[currentMolecule->componentid()
		 * totalComponents + i],
		 globalSiteADist[currentMolecule->componentid()
		 * totalComponents + i], 1, currSum, site,
		 boundary);
		 else
		 integrateRDFCartesian(xlim, ylim, zlim, currentMolecule,
		 rc, dx, dy, dz,
		 globalADist[currentMolecule->componentid()
		 * totalComponents + i],
		 globalSiteADist[currentMolecule->componentid()
		 * totalComponents + i], 1, currSum, site,
		 boundary);
		 }

		 delete currSum;
		 */
	}
	return pot;
}

double RDFForceIntegratorSite::traverseMolecules() {
	Molecule* currentMolecule;
	double force[3] = { 0, 0, 0 };
	double total_pot = 0;
	int num_bins = 400;
	double length = (_moleculeContainer->getBoundingBoxMax(0)
			- _moleculeContainer->getBoundingBoxMin(0)) / num_bins;

	//double V = (_rmax[0] - _rmin[0]) * (_rmax[1] - _rmin[1]) * length;
	double V = (_rmax[0] - _rmin[0]) * (_rmax[1] - _rmin[1]) * (_rmax[2] - _rmin[2]);
	double* rhos = new double[num_bins];
	for (int i = 1; i <= num_bins; i++) {
		double bottom[3] = { (i - 1) * length,
				_moleculeContainer->getBoundingBoxMin(1),
				_moleculeContainer->getBoundingBoxMin(2) };
		double top[3] = { i * length, _moleculeContainer->getBoundingBoxMax(1),
				_moleculeContainer->getBoundingBoxMax(2) };
		int local_num = _moleculeContainer->countParticles(
				_moleculeContainer->begin()->componentid(), bottom, top);
		rhos[i - 1] = local_num / V;
	}

	_numMolecules = _moleculeContainer->countParticles(
			_moleculeContainer->begin()->componentid(), _rmin, _rmax);
	// iterate through molecules and add rdf influence
	for (currentMolecule = _moleculeContainer->begin(); currentMolecule
			!= _moleculeContainer->end(); currentMolecule
			= _moleculeContainer->next()) {
		int rho_idx = (int) (currentMolecule->r(0) / length + 0.5);
		_rho = _numMolecules/V;//rhos[rho_idx];
		total_pot += processMolecule(currentMolecule, force, true);
	}
	return total_pot;
}

void RDFForceIntegratorSite::integrateRDFSiteCartesian(double xlim[2],
		double ylim[2], double zlim[2], Molecule* mol, int plane,
		unsigned int site, int boundary[3]) {

	std::vector<double> globalAcc = (*_globalADist)[0];
	std::vector<std::vector<double> > globalSiteAcc = (*_globalSiteADist)[0];

	// volume of the domain
	double V = (_rmax[0] - _rmin[0]) * (_rmax[1] - _rmin[1]) * (_rmax[2]
			- _rmin[2]);

	// number density of the domain
	//int numMolecules = _moleculeContainer->getNumberOfParticles();
	_numMolecules = 9826;
	double rho = _numMolecules / (V);

	// molecule position
	double molr[3] = { mol->r(0) + mol->site_d(site)[0], mol->r(1)
			+ mol->site_d(site)[1], mol->r(2) + +mol->site_d(site)[2] };

	int bin; // bin for the radius that rdf is read for
	double currf[3], absr, normal, radial, g = 0;

	// dividing part of the sphere outside the bounding box into cells of size
	// dx, dy, dz
	for (double x = xlim[0]; x + _d / 2 <= xlim[1]; x += _d) {
		for (double y = ylim[0]; y + _d / 2 <= ylim[1]; y += _d) {
			for (double z = zlim[0]; z + _d / 2 <= zlim[1]; z += _d) {
				// distance of cell center to molecule
				double r[3] = { molr[0] - x - _d / 2, molr[1] - y - _d / 2,
						molr[2] - z - _d / 2 };

				absr = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
				// check if the middle of the cell is within the cutoff radius,
				// if not, this cell will not contribute
				if (absr > _rc)
					continue;

				if (mol->numSites() == 1) {
					bin = (int) (absr * globalAcc.size() / _rc - 0.5);

					g = globalAcc[bin];
					double currPot = 0;
					PotForceLJ(r, absr * absr, 24 * mol->getEps(),
							mol->getSigma() * mol->getSigma(), currf, currPot);

					double f[3] = { 0, 0, 0 };

					for (int d = 0; d < 3; d++) {
						f[d] = rho * g * currf[d] * _d * _d * _d;
					}

					//mol->Fljcenteradd(site, f);
					if (boundary[0] == -1 && plane == 0) {
						mol->addLeftxRdfInfluence(site, f);
					}

				} else {
					// if multiple LJ centers, use site-site rdf
					// iterate through sites, treat cell center as a site

					for (unsigned int j = 0; j < mol->numLJcenters(); j++) {

						// rdf (probability) value
						bin = (int) (absr * globalSiteAcc[site
								* mol->numLJcenters() + site].size() / _rc
								- 0.5);
						if (site > j)
							g
									= globalSiteAcc[j * mol->numLJcenters()
											+ site][bin];
						else
							g
									= globalSiteAcc[site * mol->numLJcenters()
											+ j][bin];
						if (site != j)
							g /= 2;
						double currPot = 0;
						PotForceLJ(r, absr * absr, 24 * mol->getEps(),
								mol->getSigma() * mol->getSigma(), currf,
								currPot);

						//currPot *= 2 * PI * rho * g * radial * dx * dy * dz / 6;
						double f[3] = { 0, 0, 0 };

						for (int d = 0; d < 3; d++) {
							f[d] = rho * g * currf[d] * _d * _d * _d;
						}
						//if (f[0] > 0 && mol->id() == 18) cout<<"site: "<<site<<" r: "<<absr<<" value: "<<f[0]<<endl;

						//mol->Fljcenteradd(site, f);
						if (boundary[0] == -1 && plane == 0) {
							mol->addLeftxRdfInfluence(site, f);
						}

					}
				}
			}
		}

	}

}
