/*
 * RDFForceIntegratorExtendedSite.h
 *
 *  Created on: Aug 30, 2012
 *      Author: tijana
 */

#ifndef RDFFORCEINTEGRATOREXTENDEDSITE_H_
#define RDFFORCEINTEGRATOREXTENDEDSITE_H_

#include "RDFForceIntegrator.h"
#include "RDF.h"

class RDFForceIntegratorExtendedSite: public RDFForceIntegrator {
public:
	/*
	 * constructor
	 * @param moleculeContainer container for molecules
	 * @param rc cutoff radius
	 * @param d quadrature spacing
	 * @param globalADist molecule-molecule rdf (accumulated one)
	 * @param globalSiteADist site-site rdf (accumulated one)
	 */
	RDFForceIntegratorExtendedSite(ParticleContainer* moleculeContainer,
			double rc, double d, std::vector<std::vector<double> >* globalADist,
			std::vector<std::vector<std::vector<double> > >* globalSiteADist, std::string rdf_file_nondeclining);
	virtual ~RDFForceIntegratorExtendedSite();

	/*
	 * Traverses molecules and calls processMolecule for each one
	 */
	double traverseMolecules();

	/*
	 * Determines if the molecule is close to the boundary, if yes, calls integrateRDFSite
	 * @param currentMolecule molecule to which is proessed
	 * @param force here the total force on this molecule will be stored
	 * @param add_influence if the force should be added to the molecule or not (not added if just assesing the quality of the scheme)
	 * @param unit_test if it is a unit test, use linear rdf instead of the one from file
	 */
	double processMolecule(Molecule* currentMolecule, double* force, bool add_influence = true, bool unit_test = false);

private:
	double _normal_lim[2]; // limits for the integration in normal and radial direction
	double _extension; // displacement between the site and the center of mass
	double* _scaling_factors; // scaling factors according to level, normal, radial are saved here
	double	_d_alpha, _d_level; // spacing of levels and angles (used for scaling factors)
	int _n_r, _n_n, _n_levels, _n_alpha; // number of radial and normal points (used for scaling factors)
	bool called; // if precomputeScalingFactors() was called algready

	/*
	 * Integrates the force for one site
	 * @param currentMolecule: molecule that contains the site
	 * @param normal_dim: integration limits
	 * @param boundary: vector with three values, +1 means upper, -1 lower boundary, 0 it is not the boundary
	 * @param plane which direction is currently integrated for (x, y or z)
	 * @param force this will contain the total force on the molecule
	 * @param add_influence if the influence of the rdf force should be added to the molecule or not
	 */
	double integrateRDFSite(Molecule* currentMolecule, double* normal_dim,
			int* boundary, int plane, unsigned int site, double* force, bool add_influence);

	/*
	 * Called at the first traversal of the molecules. Precomputes scaling factors according
	 * to the distance of the site to the boundary, distance to the simulation domain and radial distance
	 * of the site
	 */
	void precomputeScalingFactors();

	// nondeclining rdf arrays, used for the scaling factors
    std::vector<std::vector<double> > globalNondecliningDist;
	std::vector<std::vector<double> > globalNondecliningADist;
	std::vector<std::vector<std::vector<double> > > globalNondecliningSiteDist;
	std::vector<std::vector<std::vector<double> > > globalNondecliningSiteADist;
	std::vector<double> rmids;
	std::string rdf_file_nondeclining;


};

#endif /* RDFFORCEINTEGRATOREXTENDEDSITE_H_ */
