/*
 * RDFForceIntegratorExact.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: tijana
 */

#include "RDFForceIntegratorExact.h"
bool RDFForceIntegratorExact::called = false;
double RDFForceIntegratorExact::_dr = 0;
double RDFForceIntegratorExact::_dn = 0;
double RDFForceIntegratorExact::_extension = 0;
double RDFForceIntegratorExact::_d_alpha = 0;
double RDFForceIntegratorExact::_d_level = 0;
int RDFForceIntegratorExact::_n_levels = 0;
int RDFForceIntegratorExact::_n_n = 0;
int RDFForceIntegratorExact::_n_r = 0;
int RDFForceIntegratorExact::_n_alpha = 0;
double RDFForceIntegratorExact::_dx = 0;
double RDFForceIntegratorExact::_dy = 0;
double RDFForceIntegratorExact::_dz = 0;
double RDFForceIntegratorExact::_rho = 0;
double RDFForceIntegratorExact::_g_start = 0;

double* RDFForceIntegratorExact::_scaling_factors = NULL;

std::vector<std::vector<double> >
		RDFForceIntegratorExact::globalNondecliningDist = std::vector<
				std::vector<double> >();
std::vector<std::vector<double> >
		RDFForceIntegratorExact::globalNondecliningADist = std::vector<
				std::vector<double> >();
std::vector<std::vector<std::vector<double> > >
		RDFForceIntegratorExact::globalNondecliningSiteDist = std::vector<
				std::vector<std::vector<double> > >();
std::vector<std::vector<std::vector<double> > >
		RDFForceIntegratorExact::globalNondecliningSiteADist = std::vector<
				std::vector<std::vector<double> > >();
std::vector<double> RDFForceIntegratorExact::rmids = std::vector<double>();

void RDFForceIntegratorExact::precomputeScalingFactors() {
	if (called)
		return;
	if (called == false)
		called = true;

	std::string
			rdf_file =
					"/home/tijana/Desktop/thesis/tijana/Ethan_10k_epsilon/prolonged/rdf_nondeclining/rdf_Ethan_10k_eps_double_prolonged_rc4_0-0.000090000.rdf";
	std::vector<std::string> file_names;
	file_names.push_back(rdf_file);

	// componentids that appear in the container, number of sites per such component
	std::vector<unsigned int> componentIds;
	std::vector<int> numSites;
	Molecule* currentMolecule;

	// find which component ids appear
	for (currentMolecule = _moleculeContainer->begin(); currentMolecule
			!= _moleculeContainer->end(); currentMolecule
			= _moleculeContainer->next()) {
		unsigned int i;
		for (i = 0; i < componentIds.size(); i++) {
			if (currentMolecule->componentid() == componentIds[i])
				break;
		}
		if (i == componentIds.size()) {
			componentIds.push_back(currentMolecule->componentid());
			numSites.push_back(currentMolecule->numSites());
		}
	}

	unsigned int totalComponents = componentIds.size();

	if (file_names.size() != totalComponents * totalComponents)
		std::cout << "Wrong number of rdf input files. There are "
				<< totalComponents * totalComponents
				<< " component combinations requiring the same number of rdf input files."
				<< std::endl;

	// read rdf files for different component pairs
	for (unsigned int i = 0; i < totalComponents; i++) {
		for (unsigned int j = 0; j < totalComponents; j++) {
			globalNondecliningDist.push_back(std::vector<double>());
			globalNondecliningADist.push_back(std::vector<double>());
			globalNondecliningSiteDist.push_back(std::vector<
					std::vector<double> >());
			globalNondecliningSiteADist.push_back(std::vector<std::vector<
					double> >());
			RDF::readRDFInputFile(file_names[i * totalComponents + j].c_str(),
					i, j, numSites[i], numSites[j], &rmids,
					&globalNondecliningDist[i * totalComponents + j],
					&globalNondecliningADist[i * totalComponents + j],
					&globalNondecliningSiteDist[i * totalComponents + j],
					&globalNondecliningSiteADist[i * totalComponents + j]);
		}
	}

	std::cout << "started precomputing scaling factors" << std::endl;
	//std::vector<double> globalAcc = (*_globalADist)[0];
	std::vector<std::vector<double> > globalSiteAcc = (*_globalSiteADist)[0];

	_d_level = _dn;
	_n_levels = (int) ((_rc + _extension) / _d_level + 0.5) + 1;
	_n_n = (int) (2 * _extension / _dn + 0.5) + 1;
	_d_alpha = 30;
	_n_r = (int) ((_rc + _extension) / _dr + 0.5) + 1;
	_n_alpha = (int) (360 / _d_alpha + 0.5);
	_scaling_factors = new double[_n_levels * _n_n * _n_r];

	for (int i = 0; i < _n_levels * _n_n * _n_r; i++)
		_scaling_factors[i] = -10; //initialize to some nonsense value

	double above, below, curr_g, normal, radial, curr_normal, curr_radial,
			curr_other, alpha, scale, curr_r, level, phi, site_normal,
			site_radial, site_other, site_r;
	int curr_bin;

	// level is where the molecule is
	for (int idx_level = 0; idx_level < _n_levels; idx_level++) {
		level = idx_level * _d_level;

		for (int idx_n = 0; idx_n < _n_n; idx_n++) {
			for (int idx_r = 0; idx_r < _n_r; idx_r++) {

				normal = -_extension + idx_n * _dn - _dn/2;
				radial = idx_r * _dr;
				above = 0;
				below = 0;

				for (int idx_alpha = 0; idx_alpha < _n_alpha; idx_alpha++) {
					alpha = idx_alpha * _d_alpha;

					curr_normal = normal + _extension * sin(alpha * PI / 180);

					for (int idx_phi = 0; idx_phi < _n_alpha; idx_phi++) {
						phi = idx_phi * _d_alpha;
						curr_radial = radial + _extension * cos(alpha * PI
								/ 180) * cos(phi * PI / 180);
						curr_other = _extension * cos(alpha * PI / 180) * sin(
								phi * PI / 180);

						curr_r = std::sqrt((curr_normal - level) * (curr_normal
								- level) + curr_radial * curr_radial
								+ curr_other * curr_other);

						curr_g = 0;
						if (curr_r < _rc) {

							curr_bin = (int) (curr_r
									* globalNondecliningADist[0].size() / (_rc
									+ 2 * _extension) - 0.5);
							curr_g = globalNondecliningADist[0][curr_bin];
							//std::cout<<"curr_bin "<<curr_bin<<" g "<<curr_g<<std::endl;
						}
						if (curr_normal > 0) {
							below += curr_g;
						} else {
							above += curr_g;
						}
						//std::cout<<"curr_z "<<curr_z<<" level "<<level<<" g "<<curr_g<<" curr_r "<<curr_r<<std::endl;
					}

				}
				if (above + below != 0)
					scale = above / (above + below);
				else
					scale = 0;

				_scaling_factors[idx_level * _n_n * _n_r + idx_n * _n_r + idx_r]
						= scale;
			}
		}
	}

	std::cout << "finished precomputing scaling factors" << std::endl;
}

double RDFForceIntegratorExact::checkScalingFactor(int idx_level, int idx_n,
		int idx_r) {
	double level = idx_level * _d_level;

	double normal = -_extension + idx_n * _dn - _dn/2;
	double radial = idx_r * _dr;
	double above = 0;
	double below = 0;
	//std::cout<<"in the method "<<level<<" "<<normal<<" "<<radial<<std::endl;
	int curr_bin;
	double alpha, phi, curr_normal, curr_radial, curr_other, curr_g, scale, curr_r;
	for (int idx_alpha = 0; idx_alpha < _n_alpha; idx_alpha++) {
		alpha = idx_alpha * _d_alpha;

		curr_normal = normal + _extension * sin(alpha * PI / 180);

		for (int idx_phi = 0; idx_phi < _n_alpha; idx_phi++) {
			phi = idx_phi * _d_alpha;
			curr_radial = radial + _extension * cos(alpha * PI / 180) * cos(phi
					* PI / 180);
			curr_other = _extension * cos(alpha * PI / 180) * sin(phi * PI
					/ 180);

			curr_r = std::sqrt((curr_normal - level) * (curr_normal - level)
					+ curr_radial * curr_radial + curr_other * curr_other);

			curr_g = 0;
			if (curr_r < _rc) {

				curr_bin = (int) (curr_r * globalNondecliningADist[0].size()
						/ (_rc + 2 * _extension) - 0.5);
				curr_g = globalNondecliningADist[0][curr_bin];
				//std::cout<<"curr_bin "<<curr_bin<<" g "<<curr_g<<std::endl;
			}
			if (curr_normal > 0) {
				below += curr_g;
			} else {
				above += curr_g;
			}
			//std::cout<<"curr_z "<<curr_z<<" level "<<level<<" g "<<curr_g<<" curr_r "<<curr_r<<std::endl;
		}

	}
	if (above + below != 0)
		scale = above / (above + below);
	else
		scale = 0;

	_scaling_factors[idx_level * _n_n * _n_r + idx_n * _n_r + idx_r] = scale;

}

void RDFForceIntegratorExact::getScalingFactor(double* mol_r, double* site_r,
		double x, double y, double z, int site_i, double* scale) {

	double above[2], below[2], curr_g[2], curr_x, curr_y, curr_z, curr_r,
			alpha, phi, other_site_x, other_site_y, other_site_z, other_site_r,
			ext_sin, ext_cos, sin_phi, cos_phi;
	int curr_bin;

	above[0] = above[1] = below[0] = below[1] = 0;

	for (int idx_alpha = 0; idx_alpha < _n_alpha; idx_alpha++) {
		alpha = idx_alpha * _d_alpha;
		ext_sin = _extension * sin(alpha * PI / 180);
		ext_cos = _extension * cos(alpha * PI / 180);
		curr_x = x + ext_sin;
		other_site_x = x + 2 * ext_sin;

		for (int idx_phi = 0; idx_phi < _n_alpha; idx_phi++) {
			phi = idx_phi * _d_alpha;
			cos_phi = cos(phi * PI / 180);
			sin_phi = sin(phi * PI / 180);
			curr_y = y + ext_cos * cos_phi;

			curr_z = z + ext_cos * sin_phi;

			other_site_y = y + 2 * ext_cos * cos_phi;

			other_site_z = z + 2 * ext_cos * sin_phi;

			curr_r = std::sqrt((curr_x - mol_r[0]) * (curr_x - mol_r[0])
					+ (curr_y - mol_r[1]) * (curr_y - mol_r[1]) + (curr_z
					- mol_r[2]) * (curr_z - mol_r[2]));

			other_site_r = std::sqrt((other_site_x - site_r[0]) * (other_site_x
					- site_r[0]) + (other_site_y - site_r[1]) * (other_site_y
					- site_r[1]) + (other_site_z - site_r[2]) * (other_site_z
					- site_r[2]));

			curr_g[0] = curr_g[1] = 0;
			if (curr_r < _rc) {

				curr_bin = (int) (curr_r * globalNondecliningADist[0].size()
						/ (_rc + 2 * _extension) - 0.5);
				curr_g[0] = globalNondecliningADist[0][curr_bin];//globalNondecliningSiteADist[0][site_i * 2 + 0][curr_bin];
				curr_g[1] = globalNondecliningADist[0][curr_bin];//globalNondecliningSiteADist[0][site_i * 2 + 1][curr_bin];

				//std::cout<<"curr_bin "<<curr_bin<<" g "<<curr_g<<std::endl;
			}
			if (curr_x > 0) {
				below[0] += curr_g[0];
				below[1] += curr_g[1];
			} else {
				above[0] += curr_g[0];
				above[1] += curr_g[1];
			}
			//std::cout<<"curr_z "<<curr_z<<" level "<<level<<" g "<<curr_g<<" curr_r "<<curr_r<<std::endl;
		}

	}
	if (above[0] + below[0] != 0)
		scale[0] = above[0] / (above[0] + below[0]);
	else
		scale[0] = 0;

	if (above[1] + below[1] != 0)
		scale[1] = above[1] / (above[1] + below[1]);
	else
		scale[1] = 0;

}

RDFForceIntegratorExact::RDFForceIntegratorExact(
		ParticleContainer* moleculeContainer, double rc, std::vector<
				std::vector<double> >* globalADist, std::vector<std::vector<
				std::vector<double> > >* globalSiteADist) :
	RDFForceIntegrator(moleculeContainer, rc, globalADist, globalSiteADist) {
	// TODO Auto-generated constructor stub

}

RDFForceIntegratorExact::~RDFForceIntegratorExact() {
	// TODO Auto-generated destructor stub
}

void RDFForceIntegratorExact::traverseMolecules() {
	Molecule* currentMolecule = _moleculeContainer->begin();
	_extension = currentMolecule->ljcenter_disp(0);
	_dx = _dy = _dz = _dr = _dn = currentMolecule->getSigma() / 20;
	precomputeScalingFactors();
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
	_numMolecules = 9826;
	_rho = _numMolecules / (V);

	// iterate through molecules and add rdf influence
	for (currentMolecule = _moleculeContainer->begin(); currentMolecule
			!= _moleculeContainer->end(); currentMolecule
			= _moleculeContainer->next()) {
		for (unsigned int site = 0; site < currentMolecule->numSites(); site++) {
			double rm[3] = { currentMolecule->r(0), currentMolecule->r(1),
					currentMolecule->r(2) };

			double r[3] = { currentMolecule->r(0) + currentMolecule->site_d(
					site)[0], currentMolecule->r(1) + currentMolecule->site_d(
					site)[1], currentMolecule->r(2) + currentMolecule->site_d(
					site)[2] };
			// if this is a halo molecule, skip it
			if (rm[0] < _rmin[0] || rm[1] < _rmin[1] || rm[2] < _rmin[2]
					|| rm[0] > _rmax[0] || rm[1] > _rmax[1] || rm[2] > _rmax[2])
				continue;

			// if the molecule is not within rc of the boundary, skip it
			if (rm[0] > _low_limit[0] && rm[0] < _high_limit[0] && rm[1]
					> _low_limit[1] && rm[1] < _high_limit[1] && rm[2]
					> _low_limit[2] && rm[2] < _high_limit[2])
				continue;

			if (currentMolecule->numSites() != currentMolecule->numLJcenters()) {
				std::cout
						<< "Molecule consists of something other than LJ centers. In this case RDF cannot be used.";
				return;
			}

			int boundary[3];
			if (r[0] < _low_limit[0])
				boundary[0] = -1;
			else if (r[0] >= _high_limit[0])
				boundary[0] = 1;
			if (r[1] < _low_limit[1])
				boundary[1] = -1;
			else if (r[1] >= _high_limit[1])
				boundary[1] = 1;
			if (r[2] < _low_limit[2])
				boundary[2] = -1;
			else if (r[2] >= _high_limit[2])
				boundary[2] = 1;

			//if (boundary[0] == 0 && boundary[1] == 0 && boundary[2] == 0)
			if (boundary[0] != -1)
				continue;

			// if molecule close to the boundary, add RDF force
			// integration limits for axes

			double xlim[2] = { 0, 0 };
			double ylim[2] = { 0, 0 };
			double zlim[2] = { 0, 0 };

			// box size for the numerical quadrature
			double a, h;

			if (boundary[0] == -1) {
				xlim[0] = rm[0] - _rc - _extension;
				xlim[1] = _rmin[0] + _extension;
				h = xlim[1] - xlim[0] - _extension;
				a = std::sqrt(h * (2 * (_rc + _extension) - h));
				ylim[0] = rm[1] - a;
				ylim[1] = rm[1] + a;
				zlim[0] = rm[2] - a;
				zlim[1] = rm[2] + a;

			} else if (boundary[0] == 1) {
				//normal_lim[0] = abs(_rmax[0] - r[0]);
			}

			integrateRDFSiteCartesian(xlim, ylim, zlim, currentMolecule, 0,
					site, boundary);
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
	}
}

void RDFForceIntegratorExact::integrateRDFSiteCartesian(double xlim[2],
		double ylim[2], double zlim[2], Molecule* mol, int plane,
		unsigned int site, int boundary[3]) {

	std::vector<double> globalAcc = (*_globalADist)[0];
	std::vector<std::vector<double> > globalSiteAcc = (*_globalSiteADist)[0];

	// molecule position

	double molr[3] = { mol->r(0), mol->r(1), mol->r(2) };
	double siter[3] = { mol->r(0) + mol->site_d(site)[0], mol->r(1)
			+ mol->site_d(site)[1], mol->r(2) + mol->site_d(site)[2] };

	int bin; // bin for the radius that rdf is read for
	double currf[3], absmold, allowed_dist, small_y, small_z, abssited;
	double g[2] = { 0, 0 };
	double scale[2] = { 1, 1 };
	// dividing part of the sphere outside the bounding box into cells of size
	// dx, dy, dz
	allowed_dist = std::abs(ylim[0] - molr[1]) - _extension;
	//int idx_level, idx_normal, idx_r;
	for (double x = xlim[0] + _dx / 2; x <= xlim[1]; x += _dx) {
		for (double y = ylim[0] + _dy / 2; y <= ylim[1]; y += _dy) {
			for (double z = zlim[0] + _dz / 2; z <= zlim[1]; z += _dz) {
				// distance of cell center to molecule
				double mold[3] = { molr[0] - x, molr[1] - y, molr[2] - z };

				absmold = std::sqrt(mold[0] * mold[0] + mold[1] * mold[1]
						+ mold[2] * mold[2]);

				if (absmold > _rc + _extension)
					continue;

				double sited[3] = { siter[0] - x, siter[1] - y, siter[2] - z };

				abssited = std::sqrt(sited[0] * sited[0] + sited[1] * sited[1]
						+ sited[2] * sited[2]);

				small_y = std::abs(y - molr[1]) - allowed_dist;
				small_z = std::abs(z - molr[2]) - allowed_dist;

				if (small_y > 0 && small_z > 0 && std::sqrt(small_y * small_y
						+ small_z * small_z) > _extension) {
					continue;
				}

				// check if the middle of the cell is within the cutoff radius,
				// if not, this cell will not contribute


				if (plane == 0 && std::abs(x) < _extension) {

					int idx_level = (int) (molr[0] / _dn + 0.5);
					int idx_normal = (int) ((x + _extension) / _dn + 0.5);
					int idx_r = (int) (std::sqrt((y - molr[1]) * (y - molr[1])
							+ (z - molr[2]) * (z - molr[2])) / _dn + 0.5);
					//std::cout<<"at calling "<<molr[0]<<" "<<x<<" "<<std::sqrt((y - molr[1]) * (y - molr[1])
						//	+ (z - molr[2]) * (z - molr[2]))<<std::endl;
					scale[0] = scale[1] = _scaling_factors[idx_level * _n_n
							* _n_r + idx_normal * _n_r + idx_r];
					//checkScalingFactor(idx_level, idx_normal, idx_r);
					//getScalingFactor(molr, siter, x, y, z, site, scale);

				}

				// if multiple LJ centers, use site-site rdf
				// iterate through sites, treat cell center as a site


				// rdf (probability) value
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

				double currPot = 0;
				PotForceLJ(sited, abssited * abssited, 24 * mol->getEps(),
						mol->getSigma() * mol->getSigma(), currf, currPot);

				//currPot *= 2 * PI * rho * g * radial * dx * dy * dz / 6;
				double f0[3] = { 0, 0, 0 };
				double f1[3] = { 0, 0, 0 };

				for (int d = 0; d < 3; d++) {
					f0[d] = _rho * g[0] * currf[d] * _dx * _dy * _dz;
					f1[d] = _rho * g[1] * currf[d] * _dx * _dy * _dz;
				}
				//if (f[0] > 0 && mol->id() == 18) cout<<"site: "<<site<<" r: "<<absr<<" value: "<<f[0]<<endl;

				//mol->Fljcenteradd(site, f);
				//if (boundary[0] == -1 && plane == 0) {
				mol->addLeftxRdfInfluence(site, f0);
				mol->addLeftxRdfInfluence(site, f1);
				//}

			}

		}

	}

}
/*
 * if (mol->numSites() == 1) {
 bin = (int) (abssited * globalAcc.size() / (_rc + 2
 * _extension) - 0.5);

 g = globalAcc[bin];
 double currPot = 0;
 PotForceLJ(sited, abssited * abssited, 24 * mol->getEps(),
 mol->getSigma() * mol->getSigma(), currf, currPot);

 double f[3] = { 0, 0, 0 };

 for (int d = 0; d < 3; d++) {
 f[d] = rho * g * currf[d] * _dx * _dy * _dz;
 }

 //mol->Fljcenteradd(site, f);
 if (boundary[0] == -1 && plane == 0) {
 mol->addLeftxRdfInfluence(site, f);
 }

 } */
