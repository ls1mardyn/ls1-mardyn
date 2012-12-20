/*
 * RDFForceIntegratorExact.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: tijana
 */

#include "RDFForceIntegratorExact.h"

double RDFForceIntegratorExact::getGaussianRandomNumber() {

	if (first_unif) {

		// in this case: generate new numbers
		double s = 2.0;
		double v[2] = { 0, 0 };

		while (s >= 1.0) {
			unif_rand[0] = ((double) rand()) / ((double) RAND_MAX);
			unif_rand[1] = ((double) rand()) / ((double) RAND_MAX);
			v[0] = (2.0 * unif_rand[0] - 1.0);
			v[1] = (2.0 * unif_rand[1] - 1.0);
			s = v[0] * v[0] + v[1] * v[1];
		}
		unif_rand[0] = v[0] * sqrt((-2.0 * log(s)) / s);
		unif_rand[1] = v[1] * sqrt((-2.0 * log(s)) / s);

		// change to the other variable in the next call
		first_unif = false;
		return unif_rand[0];

		// otherwise: change to the other random number
	} else {
		first_unif = true;
		return unif_rand[1];
	}
}

double* RDFForceIntegratorExact::precomputeScalingFactorsX(bool unit_test) {
	if (called_x && !unit_test)
		return _scaling_factors_x;
	if (called_x == false && !unit_test)
		called_x = true;

	std::cout << "started precomputing scaling factors" << std::endl;
	std::vector<double> globalAcc = (*_globalADist)[0];
	//std::vector<std::vector<double> > globalSiteAcc = (*_globalSiteADist)[0];

	_d_level = _d;
	_n_levels = (int) ((_rc + _extension) / _d_level + 0.5) + 1;
	_n_n = (int) (2 * _extension / _d + 0.5) + 1;

	_n_r = (int) ((_rc + _extension) / _d + 0.5) + 1;
	_n_alpha = (int) (360 / _d_alpha + 0.5);
	_scaling_factors_x = new double[_n_levels * _n_n * _n_r];

	for (int i = 0; i < _n_levels * _n_n * _n_r; i++)
		_scaling_factors_x[i] = -10; //initialize to some nonsense value

	double above, below, curr_g, normal, radial, curr_normal, curr_radial,
			curr_other, alpha, scale, curr_r, level, phi;
	int curr_bin;

	// level is where the molecule is
	for (int idx_level = 0; idx_level < _n_levels; idx_level++) {
		level = idx_level * _d_level;

		// iterate throug normal and radial direction
		for (int idx_n = 0; idx_n < _n_n; idx_n++) {
			for (int idx_r = 0; idx_r < _n_r; idx_r++) {
				normal = -_extension + idx_n * _d - _d / 2;
				radial = idx_r * _d;
				above = 0;
				below = 0;

				for (int idx_alpha = 0; idx_alpha < _n_alpha; idx_alpha++) {
					alpha = idx_alpha * _d_alpha;

					// normal coordinate of the center of mass
					curr_normal = normal + _extension * sin(alpha * PI / 180);

					for (int idx_phi = 0; idx_phi < _n_alpha; idx_phi++) {
						phi = idx_phi * _d_alpha;

						// radial and otehr coordinate of the center of mass
						curr_radial = radial + _extension * cos(
								alpha * PI / 180) * cos(phi * PI / 180);
						curr_other = _extension * cos(alpha * PI / 180) * sin(
								phi * PI / 180);

						// distance of the center of mass
						curr_r = std::sqrt(
								(curr_normal - level) * (curr_normal - level)
										+ curr_radial * curr_radial
										+ curr_other * curr_other);

						curr_g = 0;
						// if center of mass within rc from molecule, it contributes
						if (curr_r < _rc) {

							curr_bin = (int) (curr_r * globalAcc.size() / (_rc
									+ 2 * _extension) - 0.5);
							// get rdf for this center-center interaction
							curr_g = globalAcc[curr_bin];

						}
						// for unit test
						if (unit_test)
							curr_g = curr_r;

						// add contribution for below or above the boundary
						if (curr_normal > 0) {
							below += curr_g;
						} else {
							above += curr_g;
						}
					}

				}
				// compute scaling factor
				if (above + below != 0)
					scale = above / (above + below);
				else
					scale = 0;

				_scaling_factors_x[idx_level * _n_n * _n_r + idx_n * _n_r
						+ idx_r] = scale;

				// simple scaling factor for unit test, don't consider this
				if (unit_test)
					_scaling_factors_x[idx_level * _n_n * _n_r + idx_n * _n_r
							+ idx_r] = above + below;
			}
		}
	}

	return _scaling_factors_x;
}

RDFForceIntegratorExact::RDFForceIntegratorExact(
		ParticleContainer* moleculeContainer, double rc, double d,
		std::vector<std::vector<double> >* globalADist,
		std::vector<std::vector<std::vector<double> > >* globalSiteADist) :
	RDFForceIntegrator(moleculeContainer, rc, d, globalADist, globalSiteADist) {
	// TODO Auto-generated constructor stub
	called_x = false;
	_d_level = _n_levels = _n_n = _n_r = _n_alpha = _rho = 0;

	timestep = 0;

	_scaling_factors_x = NULL;
	first_unif = true;
	_extension = moleculeContainer->begin()->ljcenter_disp(0);
	_d_alpha = 5;

	std::vector<std::vector<double> > globalSiteAcc = (*_globalSiteADist)[0];

	// volume of the domain
	double V = (_rmax[0] - _rmin[0]) * (_rmax[1] - _rmin[1]) * (_rmax[2]
			- _rmin[2]);

	// number density of the domain
	_numMolecules = _moleculeContainer->countParticles(
			_moleculeContainer->begin()->componentid(), _rmin, _rmax);
	_rho = _numMolecules / (V);
	precomputeScalingFactorsX();
}

RDFForceIntegratorExact::~RDFForceIntegratorExact() {
	// TODO Auto-generated destructor stub

}
double RDFForceIntegratorExact::processMolecule(Molecule* currentMolecule,
		double* force, bool add_influence, bool unit_test) {

	// total potential on this molecule
	double pot = 0;

	// if the scaling factors were not precomputed, do that now
	if (!called_x)
		precomputeScalingFactorsX();

	// center of molecule
	double rm[3] = { currentMolecule->r(0), currentMolecule->r(1),
			currentMolecule->r(2) };

	// boundary set to -1 if left, 1 if right
	int boundary[3] = { 0, 0, 0 };
	if (rm[0] < _low_limit[0])
		boundary[0] = -1;
	else if (rm[0] > _high_limit[0])
		boundary[0] = 1;
	if (rm[1] < _low_limit[1])
		boundary[1] = -1;
	else if (rm[1] > _high_limit[1])
		boundary[1] = 1;
	if (rm[2] < _low_limit[2])
		boundary[2] = -1;
	else if (rm[2] > _high_limit[2])
		boundary[2] = 1;

	for (unsigned int site = 0; site < currentMolecule->numSites(); site++) {
		for (int d = 0; d < 3; d++)
			force[d] = 0;

		// if this is a halo molecule, skip it
		if (rm[0] < _rmin[0] || rm[1] < _rmin[1] || rm[2] < _rmin[2] || rm[0]
				> _rmax[0] || rm[1] > _rmax[1] || rm[2] > _rmax[2])
			return 0;

		// if the molecule is not within rc of the boundary, skip it
		if (rm[0] > _low_limit[0] && rm[0] < _high_limit[0] && rm[1]
				> _low_limit[1] && rm[1] < _high_limit[1] && rm[2]
				> _low_limit[2] && rm[2] < _high_limit[2])
			return 0;

		// FIXME: what do i actually do here?
		if (currentMolecule->numSites() != currentMolecule->numLJcenters()) {
			std::cout
					<< "Molecule consists of something other than LJ centers. In this case RDF cannot be used.";
			return 0;
		}

		// ATTENTION: change this if the boundary should not be in x direction
		if (boundary[0] != -1 && boundary[0] != 1)
			continue;

		// if molecule close to the boundary, add RDF force
		// integration limits for axes

		double xlim[2] = { 0, 0 };
		double ylim[2] = { 0, 0 };
		double zlim[2] = { 0, 0 };

		// box size for the numerical quadrature
		double a, h;

		// set integration limits
		if (boundary[0] == -1) {
			xlim[0] = rm[0] - _rc - _extension;
			xlim[1] = _rmin[0] + _extension;

		} else if (boundary[0] == 1) {
			xlim[0] = rm[0] + _rc + _extension;
			xlim[1] = _rmax[0] - _extension;
		}

		// hight of the spherical cap and the other dimension of it
		//h = std::min(_rc + _extension, std::abs(xlim[1] - xlim[0]));
		h = std::abs(xlim[1] - xlim[0]) - _extension;
		a = std::sqrt(h * (2 * (_rc + _extension) - h));
		ylim[0] = rm[1] - a;
		ylim[1] = rm[1] + a;
		zlim[0] = rm[2] - a;
		zlim[1] = rm[2] + a;

		// do the integration
		pot += integrateRDFSiteCartesian(xlim, ylim, zlim, currentMolecule, 0,
				site, boundary, false, force, unit_test);

		//double rnd = -1 + 2 * ((double) rand() / (double) RAND_MAX);
		//		double rnd = -1 + 2  * ((double) rand() / (double) RAND_MAX);//getGaussianRandomNumber();
		//		double randf[3] = {0, 0, 0};
		//		double other_randf[3] = {0, 0, 0};
		//		double rand_axis[3] = {(double) rand() / (double) RAND_MAX, (double) rand() / (double) RAND_MAX,(double) rand() / (double) RAND_MAX};
		//		double normalizer = std::sqrt(rand_axis[0] * rand_axis[0] + rand_axis[1] * rand_axis[1] + rand_axis[2] * rand_axis[2]);
		//		for (int d = 0; d < 3; d++)
		//			rand_axis[d]/= normalizer;
		//		double angle = 10 * (double) rand() / (double) RAND_MAX;
		//		angle = angle * 3.14/180;
		//		if (rnd > 0.9) {
		//			Quaternion q = currentMolecule->q();
		//			q.multiply_left(Quaternion(cos(angle / 2), rand_axis[0] * sin(angle / 2),
		//					rand_axis[1] * sin(angle / 2), rand_axis[2] * sin(angle / 2)));
		//			currentMolecule->setq(q);
		//		}
		//if (timestep % 10 == 0)
		//		randf[0] = rnd * force[0];
		//		other_randf[0] =  (1 - rnd) * force[0];
		//force[0] += rnd * force[0];
		//currentMolecule->Fadd(force);

		//force[0] += randf[0];


	}

	return pot;
}
double RDFForceIntegratorExact::traverseMolecules() {

	Molecule* currentMolecule = _moleculeContainer->begin();
	_extension = currentMolecule->ljcenter_disp(0);
	_d_alpha = 5;

	// total potential
	double total_pot = 0;

	// site-site rdfs
	std::vector<std::vector<double> > globalSiteAcc = (*_globalSiteADist)[0];

	// volume of the domain
	double V = (_rmax[0] - _rmin[0]) * (_rmax[1] - _rmin[1]) * (_rmax[2]
			- _rmin[2]);

	// number density of the domain

	_numMolecules = _moleculeContainer->countParticles(
			_moleculeContainer->begin()->componentid(), _rmin, _rmax);
	_rho = _numMolecules / (V);

	// iterate through molecules and add rdf influence
	for (currentMolecule = _moleculeContainer->begin(); currentMolecule
			!= _moleculeContainer->end(); currentMolecule
			= _moleculeContainer->next()) {
		double force[3] = { 0, 0, 0 };

		total_pot += this->processMolecule(currentMolecule, force, true);

	}
	timestep++;
	return total_pot;
}

double RDFForceIntegratorExact::integrateRDFSiteCartesian(double xlim[2],
		double ylim[2], double zlim[2], Molecule* mol, int plane,
		unsigned int site, int boundary[3], bool add_influence,
		double* return_force, bool unit_test) {

	// molecule - moleucle rdf
	std::vector<double> globalAcc = (*_globalADist)[0];

	// site-site rdf
	std::vector<std::vector<double> > globalSiteAcc = (*_globalSiteADist)[0];

	// molecule position
	double pot = 0;
	double molr[3] = { mol->r(0), mol->r(1), mol->r(2) };
	double siter[3] = { mol->r(0) + mol->site_d(site)[0], mol->r(1)
			+ mol->site_d(site)[1], mol->r(2) + mol->site_d(site)[2] };

	int bin; // bin for the radius that rdf is read for
	double currf[3], absmold, allowed_dist, small_y, small_z, abssited;
	double g[2] = { 0, 0 };
	double scale[2] = { 1, 1 };
	// dividing part of the sphere outside the bounding box into cells of size
	// dx, dy, dz
	if (plane == 0) {
		double h_rc = std::abs(xlim[1] - xlim[0]) - 2 * _extension;
		allowed_dist = h_rc * (2 * _rc - h_rc); // for rounded boundaries
	}

	double x;

	double rnd = (double) rand() / (double) RAND_MAX;

	// make sure it is integrated in the same direction for both upper and lower boundary
	for (x = xlim[0] - boundary[0] * _d / 2; boundary[0] * x >= boundary[0]
			* xlim[1]; x -= boundary[0] * _d) {

		for (double y = ylim[0] + _d / 2; y <= ylim[1]; y += _d) {

			for (double z = zlim[0] + _d / 2; z <= zlim[1]; z += _d) {

				// distance of cell center to molecule
				double mold[3] = { molr[0] - x, molr[1] - y, molr[2] - z };

				absmold = std::sqrt(
						mold[0] * mold[0] + mold[1] * mold[1] + mold[2]
								* mold[2]);

				// check if the middle of the cell is within the cutoff radius,
				// if not, this cell will not contribute
				if (absmold > _rc + _extension)
					continue;

				double sited[3] = { siter[0] - x, siter[1] - y, siter[2] - z };

				abssited = std::sqrt(
						sited[0] * sited[0] + sited[1] * sited[1] + sited[2]
								* sited[2]);

				// check rounded boundaries
				if (plane == 0) {
					small_y = std::abs(y - molr[1]) - allowed_dist;
					small_z = std::abs(z - molr[2]) - allowed_dist;
				}

				if (small_y > 0 && small_z > 0 && std::sqrt(
						small_y * small_y + small_z * small_z) > _extension) {
					continue;
				}

				// if in the region for scaling factors and lower x boundary, get scaling factor
				if (plane == 0 && boundary[0] == -1 && std::abs(x) < _extension
						&& !unit_test) {

					int idx_level = (int) (molr[0] / _d_level + 0.5);
					int idx_normal = (int) ((x + _extension) / _d + 0.5);
					int idx_r = (int) (std::sqrt(
							(y - molr[1]) * (y - molr[1]) + (z - molr[2]) * (z
									- molr[2])) / _d + 0.5);
					scale[0] = scale[1] = _scaling_factors_x[idx_level * _n_n
							* _n_r + idx_normal * _n_r + idx_r];

				}

				// if in the region for scaling factors and higher x boundary, get scaling factor
				if (plane == 0 && boundary[0] == 1 && std::abs(_rmax[0] - x)
						<= _extension && !unit_test) {

					int idx_level = (int) ((_rmax[0] - molr[0]) / _d_level
							+ 0.5);

					int idx_normal = (int) ((_rmax[0] - x + _extension) / _d
							+ 0.5);
					int idx_r = (int) (std::sqrt(
							(y - molr[1]) * (y - molr[1]) + (z - molr[2]) * (z
									- molr[2])) / _d + 0.5);

					scale[0] = scale[1] = _scaling_factors_x[idx_level * _n_n
							* _n_r + idx_normal * _n_r + idx_r];
				}

				// if multiple LJ centers, use site-site rdf
				// get both site-site rdfs

				bin = (int) (abssited * globalSiteAcc[site
						* mol->numLJcenters() + site].size() / (_rc + 2
						* _extension) - 0.5);

				// consider the infinitesimal element to be site 0
				g[0] = globalSiteAcc[site * mol->numLJcenters() + 0][bin];

				// consider the infinitesimal element to be site 1
				g[1] = globalSiteAcc[site * mol->numLJcenters() + 1][bin];

				if (site != 0)
					g[0] /= 2;

				if (site != 1)
					g[1] /= 2;

				g[0] *= scale[0];
				g[1] *= scale[1];

				// for the unit test, do not consider this
				if (unit_test)
					g[0] = g[1] = abssited;

				double currPot = 0;
				// get current force
				PotForceLJ(sited, abssited * abssited, 24 * mol->getEps(),
						mol->getSigma() * mol->getSigma(), currf, currPot);

				// force on the site coming from both other sites accross the boundary
				double f0[3] = { 0, 0, 0 };
				double f1[3] = { 0, 0, 0 };

				for (int d = 0; d < 1; d++) {
					f0[d] = _rho * g[0] * currf[d] * _d * _d * _d;
					f1[d] = _rho * g[1] * currf[d] * _d * _d * _d;

				}
				pot += _rho * g[0] * currPot * _d * _d * _d;
				pot += _rho * g[1] * currPot * _d * _d * _d;

				if (add_influence) {
					// add force to site
					mol->Fljcenteradd(site, f0);
					mol->Fljcenteradd(site, f1);

				}

				// add to asses the quality of the scheme
				// note that if usher is called this will be added even when usher is searching
				// for the position
				if (boundary[0] == -1 && plane == 0) {
					mol->addLeftxRdfInfluence(site, f0);
					mol->addLeftxRdfInfluence(site, f1);
				}
				for (int d = 0; d < 3; d++) {

					return_force[d] += f0[d];
					return_force[d] += f1[d];
				}

			}

		}

	}

	return pot;
}
