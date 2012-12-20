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

	double potential = 0;

	int bin; // bin for the radius that rdf is read for
	double g;
	for (int i = 0; i < 3; i++)
		force[i] = 0;
	for (double z = normal_dim[0] + _d / 2; z <= normal_dim[1]; z += _d) {
		for (double x = _d / 2; x <= std::sqrt(
				normal_dim[1] * normal_dim[1] - z * z); x += _d) {
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
				potential += 2 * PI * _rho * g * lj12_6 * x * _d * _d;
				if (add_influence) {
					mol->Fljcenteradd(site, f);
				}

				// track the force coming from the rdf boundary
				if (boundary[0] == -1 && plane == 0 && add_influence) {
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
					f[plane] = boundary[plane] * 2 * PI * _rho * g * ljf * z
							* x * _d * _d / r;
					force[plane] += f[plane];
					potential += 2 * PI * _rho * g * lj12_6 * x * _d * _d;
					if (add_influence)
						mol->Fljcenteradd(site, f);

					// track the force coming from the rdf boundary
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

		/** For efficiency, now only consider open boundaries in
		 * x direction
		 */
		if (boundary[0] != -1 && boundary[0] != 1)
			return 0;

		// if molecule close to the boundary, add RDF force
		// integration limits for axes

		double normal_lim[2] = { 0, _rc };
		double xlim[2] = { 0, 0 };
		double ylim[2] = { 0, 0 };
		double zlim[2] = { 0, 0 };
		double a, h;
		if (boundary[0] == -1) {
			normal_lim[0] = std::abs(r[0] - _rmin[0]);
			xlim[0] = r[0] - _rc;
			xlim[1] = 0;

			// bound the integration to up to rc
			if (r[0] < _rmin[0]) {
				normal_lim[1] -= _rmin[0] - r[0];
				xlim[0] += _rmin[0] - r[0];
			}

		} else if (boundary[0] == 1) {
			xlim[0] = _rmax[0];
			xlim[1] = r[0] + _rc;
			normal_lim[0] = std::abs(_rmax[0] - r[0]);

			// bound the integration to up to rc
			if (r[0] > _rmax[0]) {
				normal_lim[1] -= r[0] - _rmax[0];
				xlim[1] -= r[0] - _rmax[0];
			}

		}
		if (boundary[0] != 0) {
			h = std::abs(xlim[1] - xlim[0]);
			a = std::sqrt(h * (2 * _rc - h));
			ylim[0] = r[1] - a;
			ylim[1] = r[1] + a;

			zlim[0] = r[2] - a;
			zlim[1] = r[2] + a;
		}

		pot += integrateRDFSite(currentMolecule, normal_lim, boundary, 0, site,
				force, false);

		// this is to add the potential if the site of the molecule is outside the domain
		// the forces cancel each other out so this was used in the previous call to
		// integrate RDFSite
		// see thesis for more details
//
//		if (boundary[0] == -1 && r[0] < _rmin[0]) {
//			normal_lim[0] = 0;
//			normal_lim[1] = std::abs(r[0] - _rmin[0]);
//			pot += 2 * integrateRDFSite(currentMolecule, normal_lim, boundary,
//					0, site, force, false);
//		}
//
//		if (boundary[0] == 1 && r[0] > _rmax[0]) {
//			normal_lim[0] = 0;
//			normal_lim[1] = std::abs(r[0] - _rmax[0]);
//			pot += 2 * integrateRDFSite(currentMolecule, normal_lim, boundary,
//					0, site, force, false);
//
//		}
		//		pot += integrateRDFSiteCartesian(xlim, ylim, zlim, currentMolecule, 0, site,
		//				boundary, force, add_influence);

		// xy plane is the boundary
		// integrate xy


		/*
		 * THe code from here on is if the molecule is in vicinity of more than one boundary
		 * It was not properly tested as I always worked with boundaries in x direction
		 * I checked it with breakpoints some months ago
		 */
		/*
		if (boundary[2] != 0) {
			if (boundary[2] == -1) {
				zlim[0] = r[2] - _rc;
				zlim[1] = _rmin[2];
				normal_lim[0] = std::abs(r[2] - _rmin[2]);
				normal_lim[1] = _rc;
				if (r[2] < _rmin[2]) {
					normal_lim[1] -= _rmin[2] - r[2];
					zlim[2] += _rmin[2] - r[2];
				}
			} else if (boundary[2] == 1) {
				zlim[0] = _rmax[2];
				zlim[1] = r[2] + _rc;
				normal_lim[0] = std::abs(_rmax[0] - r[0]);
				normal_lim[1] = _rc;
				if (r[2] > _rmax[2]) {
					normal_lim[2] -= r[2] - _rmax[2];
					zlim[1] -= r[2] - _rmax[2];
				}

			}
			double h = zlim[1] - zlim[0];
			double a = std::sqrt(h * (2 * _rc - h));
			if (boundary[0] == 0)
				integrateRDFSite(currentMolecule, normal_lim, boundary, 2,
						site, force, add_influence);
			else {
				if (boundary[0] == -1) {
					xlim[0] = _rmin[0];
					xlim[1] = r[0] + a;
				} else if (boundary[0] == 1) {
					xlim[0] = r[0] - a;
					xlim[1] = _rmax[0];
				}

				ylim[0] = r[1] - a;
				ylim[1] = r[1] + a;
				integrateRDFSiteCartesian(xlim, ylim, zlim, currentMolecule, 2,
						site, boundary, force, false);
			}

		}

		// xz plane is the boundary
		// integrate xz plane
		if (boundary[1] != 0) {
			if (boundary[1] == -1) {
				ylim[0] = r[1] - _rc;
				ylim[1] = _rmin[1];
				normal_lim[0] = std::abs(r[1] - _rmin[1]);
				normal_lim[1] = _rc;
				if (r[1] < _rmin[1]) {
					normal_lim[1] -= _rmin[1] - r[1];
					ylim[2] += _rmin[1] - r[1];
				}
			} else if (boundary[1] == 1) {
				ylim[0] = _rmax[1];
				ylim[1] = r[1] + _rc;
				normal_lim[0] = std::abs(_rmax[1] - r[1]);
				normal_lim[1] = _rc;
				if (r[1] > _rmax[1]) {
					normal_lim[1] -= r[1] - _rmax[1];
					ylim[1] -= r[1] - _rmax[1];
				}
			}
			double h = ylim[1] - ylim[0];
			double a = std::sqrt(h * (2 * _rc - h));
			if (boundary[0] == -1) {
				xlim[0] = _rmin[0];
				xlim[1] = r[0] + a;
			} else if (boundary[0] == 1) {
				xlim[0] = r[0] - a;
				xlim[1] = _rmax[0];
			} else {
				xlim[0] = r[0] - a;
				xlim[1] = r[0] + a;
			}
			if (boundary[2] == -1) {
				zlim[0] = _rmin[2];
				zlim[1] = r[2] + a;
			} else if (boundary[2] == 1) {
				zlim[0] = r[2] - a;
				zlim[1] = _rmax[2];
			} else {
				zlim[0] = r[2] - a;
				zlim[1] = r[2] + a;
			}
			if (boundary[0] == 0 && boundary[2] == 0)
				integrateRDFSite(currentMolecule, normal_lim, boundary, 1,
						site, force, add_influence);
			else
				integrateRDFSiteCartesian(xlim, ylim, zlim, currentMolecule, 1,
						site, boundary, force, false);
		}
		*/
	}
	return pot;
}

double RDFForceIntegratorSite::traverseMolecules() {
	Molecule* currentMolecule;
	double force[3] = { 0, 0, 0 };
	double total_pot = 0;

	// get number density
	double V = (_rmax[0] - _rmin[0]) * (_rmax[1] - _rmin[1]) * (_rmax[2]
			- _rmin[2]);
	_numMolecules = _moleculeContainer->countParticles(
			_moleculeContainer->begin()->componentid(), _rmin, _rmax);
	_rho = _numMolecules / V;

	// iterate through molecules and add rdf influence
	for (currentMolecule = _moleculeContainer->begin(); currentMolecule
			!= _moleculeContainer->end(); currentMolecule
			= _moleculeContainer->next()) {

		total_pot += processMolecule(currentMolecule, force, true);
	}
	return total_pot;
}

double RDFForceIntegratorSite::integrateRDFSiteCartesian(double xlim[2],
		double ylim[2], double zlim[2], Molecule* mol, int plane,
		unsigned int site, int boundary[3], double* force, bool add_influence) {

	std::vector<double> globalAcc = (*_globalADist)[0];
	std::vector<std::vector<double> > globalSiteAcc = (*_globalSiteADist)[0];

	// molecule position
	double molr[3] = { mol->r(0) + mol->site_d(site)[0], mol->r(1)
			+ mol->site_d(site)[1], mol->r(2) + +mol->site_d(site)[2] };

	int bin; // bin for the radius that rdf is read for
	double currf[3], absr, normal, radial, g = 0, potential = 0;

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
						f[d] = _rho * g * currf[d] * _d * _d * _d;
						force[d] += f[d];
					}
					potential += _rho * g * currPot * _d * _d * _d;
					if (add_influence)
						mol->Fljcenteradd(site, f);
					if (boundary[0] == -1 && plane == 0) {
						mol->addLeftxRdfInfluence(site, f);
					}

				} else {
					// if multiple LJ centers, use site-site rdf
					// iterate through sites, treat cell center as a site

					for (unsigned int j = 0; j < mol->numLJcenters(); j++) {

						// rdf (probability) value
						bin = (int) (absr * globalSiteAcc[site
								* mol->numLJcenters() + site].size() / (_rc + 2
								* mol->ljcenter_disp(site)) - 0.5);
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

						potential += _rho * g * currPot * _d * _d * _d;
						double f[3] = { 0, 0, 0 };

						for (int d = 0; d < 3; d++) {
							f[d] = _rho * g * currf[d] * _d * _d * _d;
							force[d] += f[d];
						}

						if (add_influence)
							mol->Fljcenteradd(site, f);
						if (boundary[0] == -1 && plane == 0) {
							mol->addLeftxRdfInfluence(site, f);
						}

					}
				}
			}
		}

	}
	return potential;
}
