/*
 * RDFForceIntegratorExact.h
 *
 *  Created on: Aug 30, 2012
 *      Author: tijana
 */

#ifndef RDFFORCEINTEGRATOREXACT_H_
#define RDFFORCEINTEGRATOREXACT_H_

#include "RDFForceIntegrator.h"
#include "molecules/potforce.h"
#include "RDF.h"

class RDFForceIntegratorExact: public RDFForceIntegrator {
public:
	RDFForceIntegratorExact(ParticleContainer* moleculeContainer, double rc,
			std::vector<std::vector<double> >* globalADist, std::vector<
					std::vector<std::vector<double> > >* globalSiteADist);
	virtual ~RDFForceIntegratorExact();

	void traverseMolecules();

	void getScalingFactor(double* mol_r, double* site_r,
			double x, double y, double z, int site_i, double* scale);

private:
	static double _dx, _dy, _dz, _extension, _dn, _dr, *_scaling_factors, _d_alpha, _d_level;
	static int _n_r, _n_n, _n_levels, _n_alpha;
	static bool called;
	static void precomputeScalingFactors();
	static std::vector<std::vector<double> > globalNondecliningDist;
	static std::vector<std::vector<double> > globalNondecliningADist;
	static std::vector<std::vector<std::vector<double> > >
			globalNondecliningSiteDist;
	static std::vector<std::vector<std::vector<double> > >
			globalNondecliningSiteADist;
	static std::vector<double> rmids;

	void integrateRDFSiteCartesian(double xlim[2], double ylim[2],
			double zlim[2], Molecule* mol, int plane, unsigned int site,
			int boundary[3]);

};

#endif /* RDFFORCEINTEGRATOREXACT_H_ */
