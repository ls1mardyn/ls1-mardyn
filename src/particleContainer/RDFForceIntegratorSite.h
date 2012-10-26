/*
 * RDFForceIntegratorSite.h
 *
 *  Created on: Aug 30, 2012
 *      Author: tijana
 */

#ifndef RDFFORCEINTEGRATORSITE_H_
#define RDFFORCEINTEGRATORSITE_H_

#include "RDFForceIntegrator.h"
#include "molecules/potforce.h"

class RDFForceIntegratorSite: public RDFForceIntegrator {
public:
	/*
	 * Constructor
	 *
	 * @param moleculeContainer linkedcells container for the molecules
	 * @param rc cutoff radius
	 * @param d quadrature size
	 * @param globalADist vector of molecule-molecule rdfs, one vector per component (sof ar only worked with 1 component)
	 * @param globalSiteADist vectors of site-site rdfs,per component (so far considered 1 component), 0-0, 0-1, 1-0, 1-1 order
	 */
	RDFForceIntegratorSite(ParticleContainer* moleculeContainer, double rc, double d, std::vector<std::vector<double> >* globalADist,
			std::vector<std::vector<std::vector<double> > >* globalSiteADist);
	virtual ~RDFForceIntegratorSite();

	/*
	 * Traverses all the molecules close to the boundaries, computes force and potetnial
	 * @return total potential on teh molecules
	 */
	double traverseMolecules();

	/*
	 * processes one molecule
	 * @param currentMoleucle molecule to integrate for
	 * @param force output force on the molecule
	 * @param add_influence if true the force will be added to the molecule
	 * @return potential of the molecule coming from the boundaries
	 */
	double processMolecule(Molecule* currentMolecule, double* force, bool add_influence = true, bool unit_test = false);


private:
	double _rho;
	/*
	 * RDF integration in polar coordinates, should be used if only one rdf boundary
	 * @param currentMolecule molecule it is integrated for
	 * @param normal_dim integration limits
	 * @param boundary (+1 if high, -1 if low, for each dimension)
	 * @param site site it is integrated forr
	 * @param force output force on this molecule
	 * @param add_influence if true the force will be added to the molecule
	 * @return potential of the molecule coming from the boundary
	 */
	double integrateRDFSite(Molecule* currentMolecule, double* normal_dim, int* boundary, int plane, unsigned int site, double* force, bool add_influence);

	/*
	 * RDF integration in cartesian coordinates
	 * Should be used if molecule close to more than one boundary
	 * Not tested thoroughly
	 * @param xlim limits in x direction
	 * @param ylim limits in y direction
	 * @param zlim limits in z direction
	 * @param mol molecule it is integrated for
	 * @param plane which direction the boundary is in (currently integrated for)
	 * @param site which site it is integrated for
	 * @param boundary (+1 if high, -1 if low boudnary, in each direction)
	 * @param force output force on this molecule
	 * @param add_influence if true the force will be added to the molecule
	 * @return potential on this molecule coming from this boudnary
	 */
	double integrateRDFSiteCartesian(double xlim[2], double ylim[2],
			double zlim[2], Molecule* mol, int plane, unsigned int site,
			int boundary[3], double* force, bool add_influence);

};

#endif /* RDFFORCEINTEGRATORSITE_H_ */
