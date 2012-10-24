/*
 * RDFForceIntegratorExtendedSite.cpp
 *
 *  Created on: Aug 30, 2012
 *      Author: tijana
 */

#include "RDFForceIntegratorExtendedSite.h"

RDFForceIntegratorExtendedSite::RDFForceIntegratorExtendedSite(
		ParticleContainer* moleculeContainer, double rc, double d,
		std::vector<std::vector<double> >* globalADist,
		std::vector<std::vector<std::vector<double> > >* globalSiteADist, std::string rdf_file) :
	RDFForceIntegrator(moleculeContainer, rc, d, globalADist, globalSiteADist) {
	// TODO Auto-generated constructor stub

	rdf_file_nondeclining = rdf_file;

	// initialize everything to 0 or equivalent
	_n_n = _n_r = _n_levels = _n_alpha = 0;
	_extension = _d_alpha = _d_level = 0;
	_scaling_factors = NULL;
	called = false;
	globalNondecliningDist = std::vector<std::vector<double> >();
	globalNondecliningADist = std::vector<std::vector<double> >();
	globalNondecliningSiteDist
			= std::vector<std::vector<std::vector<double> > >();
	globalNondecliningSiteADist
			= std::vector<std::vector<std::vector<double> > >();
	rmids = std::vector<double>();
	Molecule* currentMolecule = _moleculeContainer->begin();
		_extension = currentMolecule->ljcenter_disp(0);
	precomputeScalingFactors();
}

RDFForceIntegratorExtendedSite::~RDFForceIntegratorExtendedSite() {
	// TODO Auto-generated destructor stub
}

void RDFForceIntegratorExtendedSite::precomputeScalingFactors() {
	// making sure this method is only called once since factors don't change
	// througout the simulation
	if (called)
		return;
	if (called == false)
		called = true;

	std::vector<std::string> file_names;
	std::cout<<"using the nondeclining rdf file: "<<rdf_file_nondeclining<<std::endl;
	file_names.push_back(rdf_file_nondeclining);

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
			globalNondecliningSiteDist.push_back(
					std::vector<std::vector<double> >());
			globalNondecliningSiteADist.push_back(
					std::vector<std::vector<double> >());
			RDF::readRDFInputFile(file_names[i * totalComponents + j].c_str(),
					i, j, numSites[i], numSites[j], &rmids,
					&globalNondecliningDist[i * totalComponents + j],
					&globalNondecliningADist[i * totalComponents + j],
					&globalNondecliningSiteDist[i * totalComponents + j],
					&globalNondecliningSiteADist[i * totalComponents + j]);
		}
	}

	std::cout << "started precomputing scaling factors" << std::endl;


	_d_level = _d; // levels spacing (level is the boundary level)
	_n_levels = (int) ((_rc + 2 * _extension) / _d_level + 0.5) + 1; // num levels
	_n_n = (int) (2 * _extension / _d + 0.5) + 1; // number in normal direction
	_n_r = (int) ((_rc + 2 * _extension) / _d + 0.5) + 1; // num radial direction
	_d_alpha = 10; // angle spacing (degrees)
	_n_alpha = (int) (360 / _d_alpha + 0.5); // num angles
	_scaling_factors = new double[_n_levels * _n_n * _n_r]; // array for the factors

	std::cout << "nn, nr, nalpha, nlevels" << _n_n << " " << _n_r << " "
			<< _n_alpha << " " << _n_levels << std::endl;
	for (int i = 0; i < _n_levels * _n_n * _n_r; i++)
		_scaling_factors[i] = -10; //initialize to some nonsense value

	double above, below, curr_g, z, x, curr_z, curr_x, curr_y, alpha, scale,
			curr_r, mol_x, mol_y, mol_r, level, mol_z, phi;
	int curr_bin;

	// site we are integration for is located at 0
	// iterate through levels (level is where the boundary is)
	for (int idx_level = 0; idx_level < _n_levels; idx_level++) {
		level = -_extension + idx_level * _d_level; // where the boundary is

		// iterate: normal distance of the imaginary site
		for (int idx_n = 0; idx_n < _n_n; idx_n++) {
			z = level - _extension + idx_n * _d + _d / 2; // site normal coord

			// iterate: radial distance of the imaginary site
			for (int idx_r = 0; idx_r < _n_r; idx_r++) {
				x = idx_r * _d + _d / 2; // site radial coord

				above = 0;
				below = 0;
				// probing the sphere, first angle
				for (int idx_alpha = 0; idx_alpha < _n_alpha; idx_alpha++) {
					alpha = idx_alpha * _d_alpha;

					curr_z = z + 2 * _extension * sin(alpha * PI / 180); // normal coord of other imaginary site
					mol_z = z + _extension * sin(alpha * PI / 180); // normal coord of the center


					// second angle
					for (int idx_phi = 0; idx_phi < _n_alpha; idx_phi++) {
						phi = idx_phi * _d_alpha; // angle in plane

						// the other two coordinates of the other site
						curr_x = x + 2 * _extension * cos(phi * PI / 180)
								* cos(alpha * PI / 180);
						curr_y = 2 * _extension * sin(phi * PI / 180) * cos(
								alpha * PI / 180);

						// the other two coordinates of the center
						mol_x = x + _extension * cos(phi * PI / 180) * cos(
								alpha * PI / 180);
						mol_y = _extension * sin(phi * PI / 180) * cos(
								alpha * PI / 180);

						// distance of the site and center
						curr_r = std::sqrt(
								curr_x * curr_x + curr_z * curr_z + curr_y
										* curr_y);
						mol_r = std::sqrt(
								mol_x * mol_x + mol_y * mol_y + mol_z * mol_z);

						curr_g = 1;

						// if site distance smaller than maximal, read rdf
						if (curr_r < _rc + 2 * _extension) {

							curr_bin = (int) (curr_r
									* globalNondecliningSiteADist[0][0].size()
									/ (_rc + 2 * _extension) - 0.5);
							curr_g
									= globalNondecliningSiteADist[0][0][curr_bin];
							//std::cout<<"curr_bin "<<curr_bin<<" g "<<curr_g<<std::endl;
						}
						// if not in the circle of radius _rc + _extension above level
						if (mol_z < level) {
							below += curr_g;
						} else {
							above += curr_g;
						}

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

double RDFForceIntegratorExtendedSite::processMolecule(
		Molecule* currentMolecule, double* force, bool add_influence,
		bool unit_test) {
	double pot = 0;
	for (unsigned int site = 0; site < currentMolecule->numSites(); site++) {
		// center (molecule) location
		double rm[3] = { currentMolecule->r(0), currentMolecule->r(1),
				currentMolecule->r(2) };

		// site location
		double r[3] = { currentMolecule->r(0)
				+ currentMolecule->site_d(site)[0], currentMolecule->r(1)
				+ currentMolecule->site_d(site)[1], currentMolecule->r(2)
				+ currentMolecule->site_d(site)[2] };

		// if this is a halo molecule, skip it
		if (rm[0] < _rmin[0] || rm[1] < _rmin[1] || rm[2] < _rmin[2] || rm[0]
				> _rmax[0] || rm[1] > _rmax[1] || rm[2] > _rmax[2])
			continue;

		// if the molecule is not within rc of the boundary, skip it
		if (rm[0] > _low_limit[0] && rm[0] < _high_limit[0] && rm[1]
				> _low_limit[1] && rm[1] < _high_limit[1] && rm[2]
				> _low_limit[2] && rm[2] < _high_limit[2])
			continue;

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
			continue;

		// if molecule close to the boundary, add RDF force
		// integration limits for axes

		double normal_lim[2] = { 0, 0 };

		if (boundary[0] == -1) {
			normal_lim[0] = r[0] - _rmin[0] - _extension;
			normal_lim[1] = _rc + _extension + r[0] - rm[0];
			if (r[0] < _rmin[0])
				normal_lim[1] -= _rmin[0] - r[0];
			//normal_lim[0] = std::abs(2 * r[0] - _rmin[0] - _extension);
		} else if (boundary[0] == 1) {
			normal_lim[0] = _rmax[0] - _extension - r[0];
			normal_lim[1] = _rc + _extension + rm[0] - r[0];
			if (r[0] > _rmax[0])
				normal_lim[1] -= r[0] - _rmax[0];
		}
		//if (currentMolecule->r(0) > 10) std::cout<<"boundary: "<<boundary[0]<<" "<<boundary[1]<<" "<<boundary[2]<<std::endl;
		double f[3] = { 0, 0, 0 };
		pot += integrateRDFSite(currentMolecule, normal_lim, boundary, 0, site, f,
				add_influence);
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
					<< " sited " << diffr << " f: " << f[0] << std::endl;
		}
		// xy plane is the boundary
	}
	return pot;
}

double RDFForceIntegratorExtendedSite::integrateRDFSite(Molecule* mol,
		double* normal_dim, int* boundary, int plane, unsigned int site,
		double* force, bool add_influence) {
	std::vector<double> globalAcc = (*_globalADist)[0];
	std::vector<std::vector<double> > globalSiteAcc = (*_globalSiteADist)[0];
	double potential = 0;
	double V = (_rmax[0] - _rmin[0]) * (_rmax[1] - _rmin[1]) * (_rmax[2]
			- _rmin[2]);

	// number density of the domain
	//int numMolecules = _moleculeContainer->getNumberOfParticles();

	double rho = _numMolecules / (V);

	double r, level, scale = 1;

	int bin; // bin for the radius that rdf is read for
	double g;

	// iterate thorugh normal direction
	for (double z = normal_dim[0] + _d / 2; z <= normal_dim[1]; z += _d) {

		// iterate through radial direction
		for (double x = _d / 2; x <= std::sqrt(
				normal_dim[1] * normal_dim[1] - z * z); x += _d) {
			// distance between the imaginary and real site
			r = std::sqrt(x * x + z * z);
			level = normal_dim[0] + _extension; // boundary level (site is always located at 0)
			scale = 1;

			// check if in the region that needs scaling
			if (z - normal_dim[0] < 2 * _extension) {

				int idx_r = (int) (x / _d + 0.5); // radial distance of the site index
				int idx_level = (int) ((_extension + level) / _d_level + 0.5); // boundary level index
				int idx_n = (int) ((z - normal_dim[0]) / _d); // normal distance index
				scale = _scaling_factors[idx_level * _n_n * _n_r + idx_n * _n_r
						+ idx_r];

			}

			// if single-centered moleucle
			if (mol->numSites() == 1) {
				// read moleucle-molecule rdf
				bin = (int) (r * globalAcc.size() / (_rc + 2 * _extension)
						- 0.5);
				g = globalAcc[bin];
				g *= scale;

				double sig2 = mol->getSigma() * mol->getSigma();
				double r2 = r * r;
				double lj6 = sig2 * sig2 * sig2 / (r2 * r2 * r2);
				double ljf = 24 * mol->getEps() * (lj6 - 2 * lj6 * lj6) / r;
				double lj12_6 = 6 * 4 * mol->getEps() * (lj6 * lj6 - lj6);

				double f[3] = { 0, 0, 0 };
				// compute force
				f[plane] = boundary[0] * 2 * PI * rho * g * ljf * z * x * _d
						* _d / r;
				force[plane] += f[plane];

				// add to molecule
				if (add_influence)
					mol->Fljcenteradd(site, f);
				potential += 2 * PI * rho * g * lj12_6 * x
															* _d * _d;
				// add for assesing scheme quality
				if (boundary[0] == 1 && plane == 0) {
					mol->addLeftxRdfInfluence(site, f);
				}

			// if multicentered molecule
			} else {

				bin = (int) (r * globalSiteAcc[0].size() / (_rc + 2
						* _extension) - 0.5);

				// if multiple LJ centers, use site-site rdf
				// iterate through sites, treat cell center as a site

				for (unsigned int j = 0; j < mol->numLJcenters(); j++) {

					g = globalSiteAcc[site * mol->numLJcenters() + j][bin];
					//}
					if (site != j)
						g /= 2;

					g *= scale;

					double sig2 = mol->getSigma() * mol->getSigma();
					double r2 = r * r;
					double lj6 = sig2 * sig2 * sig2 / (r2 * r2 * r2);

					double ljf = 24 * mol->getEps() * (lj6 - 2 * lj6 * lj6) / r;

					double f[3] = { 0, 0, 0 };
					double lj12_6 = 6 * 4 * mol->getEps() * (lj6 * lj6 - lj6);
					potential += 2 * PI * rho * g * lj12_6 * x
											* _d * _d;
					f[plane] = boundary[0] * 2 * PI * rho * g * ljf * z * x
							* _d * _d / r;
					force[plane] += f[plane];
					if (add_influence) {
						mol->Fljcenteradd(site, f);

					}

					if (boundary[0] == -1 && plane == 0) {
						mol->addLeftxRdfInfluence(site, f);
					}
				}

			}
		}
	}
	return potential;
}

double RDFForceIntegratorExtendedSite::traverseMolecules() {
	Molecule* currentMolecule = _moleculeContainer->begin();
	_extension = currentMolecule->ljcenter_disp(0);
	// box size for the numerical quadrature
	//_dn = _dr = currentMolecule->getSigma() / 20;
	// iterate through molecules and add rdf influence
	double total_pot = 0;
	precomputeScalingFactors();
	_numMolecules = _moleculeContainer->countParticles(
			_moleculeContainer->begin()->componentid(), _rmin, _rmax);
	std::cout << "nn, nr, nalpha, nlevels" << _n_n << " " << _n_r << " "
			<< _n_alpha << " " << _n_levels << std::endl;
	for (currentMolecule = _moleculeContainer->begin(); currentMolecule
			!= _moleculeContainer->end(); currentMolecule
			= _moleculeContainer->next()) {
		double force[3] = {0,0,0};
		total_pot += this->processMolecule(currentMolecule, force, true, false);
	}
	return total_pot;

}

