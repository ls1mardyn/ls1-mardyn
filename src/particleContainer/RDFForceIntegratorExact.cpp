/*
 * RDFForceIntegratorExact.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: tijana
 */

#include "RDFForceIntegratorExact.h"


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
			curr_other, alpha, scale, curr_r, level, phi, site_normal,
			site_radial, site_other, site_r;
	int curr_bin;

	// level is where the molecule is
	for (int idx_level = 0; idx_level < _n_levels; idx_level++) {
		level = idx_level * _d_level;

		for (int idx_n = 0; idx_n < _n_n; idx_n++) {
			for (int idx_r = 0; idx_r < _n_r; idx_r++) {
				normal = -_extension + idx_n * _d - _d / 2;
				radial = idx_r * _d;
				above = 0;
				below = 0;

				for (int idx_alpha = 0; idx_alpha < _n_alpha; idx_alpha++) {
					alpha = idx_alpha * _d_alpha;

					curr_normal = normal + _extension * sin(alpha * PI / 180);

					for (int idx_phi = 0; idx_phi < _n_alpha; idx_phi++) {
						phi = idx_phi * _d_alpha;
						curr_radial = radial + _extension * cos(
								alpha * PI / 180) * cos(phi * PI / 180);
						curr_other = _extension * cos(alpha * PI / 180) * sin(
								phi * PI / 180);

						curr_r = std::sqrt(
								(curr_normal - level) * (curr_normal - level)
										+ curr_radial * curr_radial
										+ curr_other * curr_other);

						curr_g = 0;
						if (curr_r < _rc) {

							curr_bin = (int) (curr_r * globalAcc.size() / (_rc
									+ 2 * _extension) - 0.5);
							curr_g = globalAcc[curr_bin];

						}
						if (unit_test)
							curr_g = curr_r;
						if (curr_normal > 0) {
							below += curr_g;
						} else {
							above += curr_g;
						}

						//std::cout<<"curr_z "<<curr_normal<<" level "<<level<<" g "<<curr_g<<" curr_r "<<curr_r<<std::endl;
					}

				}
				if (above + below != 0)
					scale = above / (above + below);
				else
					scale = 0;

				_scaling_factors_x[idx_level * _n_n * _n_r + idx_n * _n_r
						+ idx_r] = scale;
				if (unit_test)
					_scaling_factors_x[idx_level * _n_n * _n_r + idx_n * _n_r
							+ idx_r] = above + below;
			}
		}
	}

	std::cout << "finished precomputing scaling factors" << std::endl;
	return _scaling_factors_x;
}

double RDFForceIntegratorExact::checkScalingFactor(int idx_level, int idx_n,
		int idx_r) {

}

void RDFForceIntegratorExact::getScalingFactor(double* mol_r, double* site_r,
		double x, double y, double z, int site_i, double* scale) {

}

RDFForceIntegratorExact::RDFForceIntegratorExact(
		ParticleContainer* moleculeContainer, double rc, double d,
		std::vector<std::vector<double> >* globalADist,
		std::vector<std::vector<std::vector<double> > >* globalSiteADist) :
	RDFForceIntegrator(moleculeContainer, rc, d, globalADist, globalSiteADist) {
	// TODO Auto-generated constructor stub
	called_x = false;
	_extension = 0;
	_d_alpha = 0;
	_d_level = 0;
	_n_levels = 0;
	_n_n = 0;
	_n_r = 0;
	_n_alpha = 0;
	_rho = 0;
	_g_start = 0;
	timestep = 0;
	rhos = NULL;

	_scaling_factors_x = NULL;
}

RDFForceIntegratorExact::~RDFForceIntegratorExact() {
	// TODO Auto-generated destructor stub

}
double RDFForceIntegratorExact::processMolecule(Molecule* currentMolecule,
		double* force, bool add_influence, bool unit_test) {

	// total potential on this molecule
	double pot = 0;

	// this is the total force on this molecule, initialize to zero


	// if the scaling factors were not precomputed, do that now
	if (!called_x)
		precomputeScalingFactorsX();

	// center of molecule
	double rm[3] = { currentMolecule->r(0), currentMolecule->r(1),
			currentMolecule->r(2) };

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
		// current site
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

		// FIXME: what do i actually do here?
		if (currentMolecule->numSites() != currentMolecule->numLJcenters()) {
			std::cout
					<< "Molecule consists of something other than LJ centers. In this case RDF cannot be used.";
			return 0;
		}

		// set where the boundary is


		//if (boundary[0] == 0 && boundary[1] == 0 && boundary[2] == 0)
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
			//std::cout<<"xlim1 "<<xlim[1]<<" xlim0: "<<xlim[0]<<std::endl;

		}

		h = std::min(_rc + _extension, std::abs(xlim[1] - xlim[0]));
		a = std::sqrt(h * (2 * (_rc + _extension) - h));
		ylim[0] = rm[1] - a;
		ylim[1] = rm[1] + a;
		zlim[0] = rm[2] - a;
		zlim[1] = rm[2] + a;

		// do the integration
		pot += integrateRDFSiteCartesian(xlim, ylim, zlim, currentMolecule, 0,
				site, boundary, add_influence, force, unit_test);

//		double rnd = -1 + 2 * ((double) rand() / (double) RAND_MAX);
//
//		double randf[3] = {0, 0, 0};
//		randf[0] = rnd * 0.1 * force[0];
//		currentMolecule->Fljcenteradd(site, randf);
//		force[0] += randf[0];


		//		if (force[0] == 0 && currentMolecule->id() == 4789) {
		//			std::cout<<currentMolecule->id()<<"rm: "<<rm[0]<<" r: "<<r[0]<<" xlim: "<<xlim[0]<<" "<<xlim[1]<<" ylim: "<<ylim[0]<<" "<<ylim[1]<<" zlim: "<<zlim[0]<<" "<<zlim[1]<<std::endl;
		//		}

		if (rm[0] < 0.01 || rm[0] > _rmax[0] - 0.01) {
			double diff, diffr;
			if (boundary[0] == -1) {
				diff = rm[0];
				diffr = r[0];
			}
			if (boundary[0] == 1) {
				diff = _rmax[0] - rm[0];
				diffr = _rmax[0] - r[0];
			}
			std::cout << "rm: " << diff << " boundary " << boundary[0]
					<< " sited " << diffr << " f: " << force[0] << " xlim "
					<< xlim[0] << " " << xlim[1] << " diffx " << h << std::endl;
		}
	}

	return pot;
}
double RDFForceIntegratorExact::traverseMolecules() {

	Molecule* currentMolecule = _moleculeContainer->begin();
	_extension = currentMolecule->ljcenter_disp(0);
	_d_alpha = 5;
	//_dx = _dy = _dz = _dr = _dn = currentMolecule->getSigma() / 5;
	double total_pot = 0;
	std::vector<std::vector<double> > globalSiteAcc = (*_globalSiteADist)[0];

	unsigned int num = 0;
	for (unsigned int i = 0; i < globalSiteAcc[0].size(); i++) {
		if (globalSiteAcc[0][i] != 0)
			break;
		num++;
	}
	_g_start = num * ((_rc + 2 * _extension) / globalSiteAcc[0].size()) - 2
			* _extension;
	std::cout << _g_start << std::endl;

	// volume of the domain
	double V = (_rmax[0] - _rmin[0]) * (_rmax[1] - _rmin[1]) * (_rmax[2]
			- _rmin[2]);

	// number density of the domain
	//int numMolecules = _moleculeContainer->getNumberOfParticles();
	_numMolecules = _moleculeContainer->countParticles(
			_moleculeContainer->begin()->componentid(), _rmin, _rmax);
	_rho = _numMolecules / (V);
	std::cout << "num " << _numMolecules << std::endl;

	int num_bins = 800;
	int sample = 50;
	if (timestep % sample == 0) {
		rhos = new double[num_bins];
		for (int i = 0; i < num_bins; i++)
			rhos[i] = 0;
	}

	double length = (_moleculeContainer->getBoundingBoxMax(0)
			- _moleculeContainer->getBoundingBoxMin(0)) / num_bins;
	double binV = length * (_rmax[1] - _rmin[1]) * (_rmax[2] - _rmin[2]);
	srand ( (unsigned)time(NULL));
	for (int i = 1; i <= num_bins; i++) {
		double bottom[3] = { (i - 1) * length,
				_moleculeContainer->getBoundingBoxMin(1),
				_moleculeContainer->getBoundingBoxMin(2) };
		double top[3] = { i * length, _moleculeContainer->getBoundingBoxMax(1),
				_moleculeContainer->getBoundingBoxMax(2) };
		int local_num = _moleculeContainer->countParticles(
				_moleculeContainer->begin()->componentid(), bottom, top);
		rhos[i - 1] += local_num;

	}

	// iterate through molecules and add rdf influence
	for (currentMolecule = _moleculeContainer->begin(); currentMolecule
			!= _moleculeContainer->end(); currentMolecule
			= _moleculeContainer->next()) {
		double force[3] = { 0, 0, 0 };
		int idx = (int) (currentMolecule->r(0) / length + 0.5);
		//		if (timestep % (sample - 1) == 0 && timestep != 0)
		//			_rho = rhos[idx] / (sample * binV);
		//		if (idx == 799) {
		//			std::cout << _rho * binV << std::endl;
		//		}
		total_pot += this->processMolecule(currentMolecule, force, true);

	}
	timestep++;
	return total_pot;
}

double RDFForceIntegratorExact::integrateRDFSiteCartesian(double xlim[2],
		double ylim[2], double zlim[2], Molecule* mol, int plane,
		unsigned int site, int boundary[3], bool add_influence,
		double* return_force, bool unit_test) {

	std::vector<double> globalAcc = (*_globalADist)[0];
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
		allowed_dist = h_rc * (2 * _rc - h_rc);//std::abs(ylim[0] - molr[1]) - _extension;
	}

	int cnt = 0;
	double x;

	double rnd = (double) rand() / (double) RAND_MAX;

	// make sure it is integrated in the same direction for both upper and lower boundary
	for (x = xlim[0] - boundary[0] * _d / 2; boundary[0] * x >= boundary[0]
			* xlim[1]; x -= boundary[0] * _d) {

		cnt++;
		int cnt2 = 0;
		for (double y = ylim[0] + _d / 2; y <= ylim[1]; y += _d) {

			for (double z = zlim[0] + _d / 2; z <= zlim[1]; z += _d) {
				cnt2++;
				// distance of cell center to molecule
				double mold[3] = { molr[0] - x, molr[1] - y, molr[2] - z };

				absmold = std::sqrt(
						mold[0] * mold[0] + mold[1] * mold[1] + mold[2]
								* mold[2]);

				if (absmold > _rc + _extension)
					continue;

				double sited[3] = { siter[0] - x, siter[1] - y, siter[2] - z };

				abssited = std::sqrt(
						sited[0] * sited[0] + sited[1] * sited[1] + sited[2]
								* sited[2]);

				if (plane == 0) {
					small_y = std::abs(y - molr[1]) - allowed_dist;
					small_z = std::abs(z - molr[2]) - allowed_dist;
				}

				if (small_y > 0 && small_z > 0 && std::sqrt(
						small_y * small_y + small_z * small_z) > _extension) {
					continue;
				}

				// check if the middle of the cell is within the cutoff radius,
				// if not, this cell will not contribute


				if (plane == 0 && boundary[0] == -1 && std::abs(x) < _extension
						&& !unit_test) {

					int idx_level = (int) (molr[0] / _d_level + 0.5);
					int idx_normal = (int) ((x + _extension) / _d + 0.5);
					int idx_r = (int) (std::sqrt(
							(y - molr[1]) * (y - molr[1]) + (z - molr[2]) * (z
									- molr[2])) / _d + 0.5);
					scale[0] = scale[1] = _scaling_factors_x[idx_level * _n_n
							* _n_r + idx_normal * _n_r + idx_r];

					//checkScalingFactor(idx_level, idx_normal, idx_r);
					//getScalingFactor(molr, siter, x, y, z, site, scale);

				}

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
					//					std::cout << "scale pos: "
					//							<< _scaling_factors_xpos[idx_level * _n_n * _n_r
					//									+ idx_normal * _n_r + idx_r] << " "
					//							<< "scale neg: " << _scaling_factors_xneg[idx_level
					//							* _n_n * _n_r + idx_normal * _n_r + idx_r]
					//							<< std::endl;

					//checkScalingFactor(idx_level, idx_normal, idx_r);
				}

				// if multiple LJ centers, use site-site rdf
				// iterate through sites, treat cell center as a site


				// rdf (probability) value
				//				if (mol->id() == 4789 && scale[0] != 0) {
				//					std::cout<<scale[0]<<std::endl;
				//				}
				bin = (int) (abssited * globalSiteAcc[site
						* mol->numLJcenters() + site].size() / (_rc + 2
						* _extension) - 0.5);

				g[0] = globalSiteAcc[site * mol->numLJcenters() + 0][bin];
				g[1] = globalSiteAcc[site * mol->numLJcenters() + 1][bin];
				if (site != 0)
					g[0] /= 2;

				if (site != 1)
					g[1] /= 2;

				g[0] *= scale[0];
				g[1] *= scale[1];

				if (unit_test)
					g[0] = g[1] = abssited;
				double currPot = 0;
				PotForceLJ(sited, abssited * abssited, 24 * mol->getEps(),
						mol->getSigma() * mol->getSigma(), currf, currPot);

				//currPot *= 2 * PI * rho * g * radial * dx * dy * dz / 6;
				double f0[3] = { 0, 0, 0 };
				double f1[3] = { 0, 0, 0 };

				for (int d = 0; d < 1; d++) {
					f0[d] = _rho * g[0] * currf[d] * _d * _d * _d;
					f1[d] = _rho * g[1] * currf[d] * _d * _d * _d;

				}
				pot += _rho * g[0] * currPot * _d * _d * _d;
				pot += _rho * g[1] * currPot * _d * _d * _d;
				//if (f[0] > 0 && mol->id() == 18) cout<<"site: "<<site<<" r: "<<absr<<" value: "<<f[0]<<endl;
				if (add_influence) {


//					double f0_rest[3] = {f0[0], f0[1], f0[2]};
//					double f1_rest[3] = {f1[0], f1[1], f1[2]};
//					f0[0] *= rnd;
//					f0_rest[0] *= 1 - rnd;
//					f1[0] *= rnd;
//					f1_rest[0] *= 1 - rnd;
					mol->Fljcenteradd(site, f0);
					mol->Fljcenteradd(site, f1);
//					mol->Fljcenteradd(0, f0);
//					mol->Fljcenteradd(0, f1);
//					mol->Fljcenteradd(1, f0_rest);
//					mol->Fljcenteradd(1, f1_rest);

					//					mol->Fadd(f0);
					//					mol->Fadd(f1);
				}

				if (boundary[0] == -1 && plane == 0) {
					mol->addLeftxRdfInfluence(site, f0);
					mol->addLeftxRdfInfluence(site, f1);
				}
				for (int d = 0; d < 3; d++) {

					return_force[d] += f0[d];
					return_force[d] += f1[d];
				}

				//				if (cnt == 18 && cnt2 == 120 && scale[0] != 1 && scale[0] != 0
				//						&& (molr[0] < 0.01 || molr[0] > _rmax[0] - 0.01)) {
				//					double diff, diffr;
				//					if (boundary[0] == -1) {
				//						diff = molr[0];
				//						diffr = siter[0];
				//					}
				//					if (boundary[0] == 1) {
				//						diff = _rmax[0] - molr[0];
				//						diffr = _rmax[0] - siter[0];
				//					}
				//					std::cout << "rm: " << diff << " boundary " << boundary[0]
				//							<< " mold " << absmold << " siter " << diffr
				//							<< " realsiter " << siter[0] <<" sited[0] " << sited[0] << " abssited "
				//							<< abssited << " x " << x << " g "
				//							<< g[site] / scale[0] <<" bin " << bin << " f0 " <<f0[0] << std::endl;
				//				}

			}

		}

	}

	return pot;
}
