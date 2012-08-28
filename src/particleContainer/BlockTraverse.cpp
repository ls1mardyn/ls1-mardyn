/***************************************************************************
 *   Copyright (C) 2010 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include <vector>
#include <cmath>

#include "BlockTraverse.h"
#include "molecules/Molecule.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "ParticleCell.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "RDF.h"
#include "molecules/potforce.h"

using namespace std;
using Log::global_log;

//################################################
//############ PUBLIC METHODS ####################
//################################################


void BlockTraverse::processCell(ParticleCell& cell, double& cutoffRadiusSquare,
		double& LJCutoffRadiusSquare, double& tersoffCutoffRadiusSquare,
		ParticlePairsHandler* particlePairsHandler) {
	vector<Molecule*>::iterator molIter1;
	vector<Molecule*>::iterator molIter2;
	double distanceVector[3];

	for (molIter1 = cell.getParticlePointers().begin(); molIter1
			!= cell.getParticlePointers().end(); molIter1++) {
		Molecule& molecule1 = **molIter1;
		unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching
		molIter2 = molIter1;
		molIter2++; // no self interaction

		for (; molIter2 != cell.getParticlePointers().end(); molIter2++) {
			Molecule& molecule2 = **molIter2;
			assert(&molecule1 != &molecule2);
			double dd = molecule2.dist2(molecule1, distanceVector);
			double force[3] = { 0, 0, 0 };
			if (dd < cutoffRadiusSquare) {
				particlePairsHandler->processPair(molecule1, molecule2,
						distanceVector, MOLECULE_MOLECULE, dd, (dd
								< LJCutoffRadiusSquare), force);
				if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd
						< tersoffCutoffRadiusSquare)) {
					particlePairsHandler->preprocessTersoffPair(molecule1,
							molecule2, false);
				}
			}
		}
	}
}

void BlockTraverse::processCellPair(ParticleCell &cell1, ParticleCell& cell2,
		double& cutoffRadiusSquare, double& LJCutoffRadiusSquare,
		double& tersoffCutoffRadiusSquare,
		ParticlePairsHandler* particlePairsHandler) {
	vector<Molecule*>::iterator molIter1;
	vector<Molecule*>::iterator molIter2;
	double distanceVector[3];

	for (molIter1 = cell1.getParticlePointers().begin(); molIter1
			!= cell1.getParticlePointers().end(); molIter1++) {
		Molecule& molecule1 = **molIter1;
		unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

		for (molIter2 = cell2.getParticlePointers().begin(); molIter2
				!= cell2.getParticlePointers().end(); molIter2++) {
			Molecule& molecule2 = **molIter2;
			assert(&molecule1 != &molecule2);
			double dd = molecule2.dist2(molecule1, distanceVector);
			double force[3] = { 0, 0, 0 };
			if (dd < cutoffRadiusSquare) {
				particlePairsHandler->processPair(molecule1, molecule2,
						distanceVector, MOLECULE_MOLECULE, dd, (dd
								< LJCutoffRadiusSquare), force);
				if ((num_tersoff > 0) && (molecule2.numTersoff() > 0) && (dd
						< tersoffCutoffRadiusSquare)) {
					particlePairsHandler->preprocessTersoffPair(molecule1,
							molecule2, false);
				}
			}
		}
	}
}

BlockTraverse::BlockTraverse(ParticleContainer* moleculeContainer, vector<
		ParticleCell>& cells, vector<unsigned long>& innerCellIndices, vector<
		unsigned long>& boundaryCellIndices,
		vector<unsigned long>& haloCellIndices,
		vector<vector<unsigned long> >& forwardNeighbourOffsets, vector<vector<
				unsigned long> >& backwardNeighbourOffsets) :
	_moleculeContainer(moleculeContainer), _cells(cells), _innerCellIndices(
			innerCellIndices), _boundaryCellIndices(boundaryCellIndices),
			_haloCellIndices(haloCellIndices), _forwardNeighbourOffsets(
					&forwardNeighbourOffsets), _backwardNeighbourOffsets(
					&backwardNeighbourOffsets), _allocatedOffsets(false) {
}

BlockTraverse::BlockTraverse(ParticleContainer* moleculeContainer, vector<
		ParticleCell>& cells, vector<unsigned long>& innerCellIndices, vector<
		unsigned long>& boundaryCellIndices,
		vector<unsigned long>& haloCellIndices) :
	_moleculeContainer(moleculeContainer), _cells(cells), _innerCellIndices(
			innerCellIndices), _boundaryCellIndices(boundaryCellIndices),
			_haloCellIndices(haloCellIndices), _forwardNeighbourOffsets(0),
			_backwardNeighbourOffsets(0), _allocatedOffsets(true) {
	_forwardNeighbourOffsets = new vector<vector<unsigned long> > ;
	_backwardNeighbourOffsets = new vector<vector<unsigned long> > ;
}

BlockTraverse::~BlockTraverse() {
	if (_allocatedOffsets) {
		delete _forwardNeighbourOffsets;
		delete _backwardNeighbourOffsets;
	}
}

void BlockTraverse::assignOffsets(
		vector<unsigned long>& forwardNeighbourOffsets,
		vector<unsigned long>& backwardNeighbourOffsets) {
	_forwardNeighbourOffsets->assign(_cells.size(), forwardNeighbourOffsets);
	_backwardNeighbourOffsets->assign(_cells.size(), backwardNeighbourOffsets);
}

double BlockTraverse::integrateRDFSite(double normal_dim[2], Molecule* mol,
		double rc, double dz, double dx, vector<double> globalAcc, vector<
				vector<double> > globalSiteAcc, int plane, unsigned int site,
		int boundary[3]) {

	// volume of the domain
	double V = (_moleculeContainer->getBoundingBoxMax(0)
			- _moleculeContainer->getBoundingBoxMin(0))
			* (_moleculeContainer->getBoundingBoxMax(1)
					- _moleculeContainer->getBoundingBoxMin(1))
			* (_moleculeContainer->getBoundingBoxMax(2)
					- _moleculeContainer->getBoundingBoxMin(2));

	// number density of the domain
	//int numMolecules = _moleculeContainer->getNumberOfParticles();
	int numMolecules = 18522;
	double rho = numMolecules / (V);

	// potential
	double pot = 0;

	int bin; // bin for the radius that rdf is read for
	double g;
	for (double z = normal_dim[0] + dz / 2; z <= rc; z += dz) {
		for (double x = dx / 2; x <= std::sqrt(rc * rc - z * z); x += dx) {
			double r = std::sqrt(x * x + z * z);
			double currPot = 0;
			if (mol->numSites() == 1) {
				bin = (int) (r * globalAcc.size() / rc - 0.5);
				g = globalAcc[bin];
				double sig2 = mol->getSigma() * mol->getSigma();
				double r2 = r * r;
				double lj6 = sig2 * sig2 * sig2 / (r2 * r2 * r2);
				double lj = 4 * mol->getEps() * (lj6 * lj6 - lj6);
				double ljf = 24 * mol->getEps() * (lj6 - 2 * lj6 * lj6) / r;
				currPot = 2 * PI * rho * g * lj * x * dx * dz;
				double f[3] = { 0, 0, 0 };
				f[plane] = -2 * PI * rho * g * ljf * z * x * dx * dz / r;
				//mol->Fljcenteradd(site, f);

				if (boundary[0] == -1 && plane == 0) {
					mol->addLeftxRdfInfluence(site, f);
				}
			} else {
				bin = (int) (r * globalSiteAcc[0].size() / rc - 0.5);

				// if multiple LJ centers, use site-site rdf
				// iterate through sites, treat cell center as a site

				for (unsigned int j = 0; j < mol->numLJcenters(); j++) {

					g = globalSiteAcc[site * mol->numLJcenters() + j][bin];

					if (site != j)
						g /= 2;

					double sig2 = mol->getSigma() * mol->getSigma();
					double r2 = r * r;
					double lj6 = sig2 * sig2 * sig2 / (r2 * r2 * r2);
					double lj = 4 * mol->getEps() * (lj6 * lj6 - lj6);
					double ljf = 24 * mol->getEps() * (lj6 - 2 * lj6 * lj6) / r;
					currPot = 2 * PI * rho * g * lj * x * dx * dz;
					double f[3] = { 0, 0, 0 };
					f[plane] = -2 * PI * rho * g * ljf * z * x * dx * dz / r;
					//mol->Fljcenteradd(site, f);

					if (boundary[0] == -1 && plane == 0) {
						mol->addLeftxRdfInfluence(site, f);
					}
					pot += currPot;
				}

			}
		}
	}
	return pot;
}

void BlockTraverse::traverseRDFBoundaryCartesian(
		vector<vector<double> >* globalADist,
		vector<vector<vector<double> > >* globalSiteADist,
		ParticlePairsHandler* particlePairsHandler) {
	double rc = _moleculeContainer->getCutoff();
	// placeholder
	vector<double> rmids;

	int countBoundary = 0;
	int countTotal = 0;

	double rmin[3]; // lower corner of the process-specific domain //ENABLE_MPI
	double rmax[3];
	double low_limit[3];
	double high_limit[3];
	for (int i = 0; i < 3; i++) {
		rmin[i] = _moleculeContainer->getBoundingBoxMin(i);
		low_limit[i] = rmin[i] + rc;
		rmax[i] = _moleculeContainer->getBoundingBoxMax(i);
		high_limit[i] = rmax[i] - rc;
	}
	Molecule* currentMolecule;

	// iterate through molecules and add rdf influence
	for (currentMolecule = _moleculeContainer->begin(); currentMolecule
			!= _moleculeContainer->end(); currentMolecule
			= _moleculeContainer->next()) {
		for (unsigned int site = 0; site < currentMolecule->numSites(); site++) {
			double rm0 = currentMolecule->r(0);
			double rm1 = currentMolecule->r(1);
			double rm2 = currentMolecule->r(2);
			double r0 = currentMolecule->r(0)
					+ currentMolecule->site_d(site)[0];
			double r1 = currentMolecule->r(1)
					+ currentMolecule->site_d(site)[1];
			double r2 = currentMolecule->r(2)
					+ currentMolecule->site_d(site)[2];
			// if this is a halo molecule, skip it
			if (rm0 < rmin[0] || rm1 < rmin[1] || rm2 < rmin[2] || rm0
					> rmax[0] || rm1 > rmax[1] || rm2 > rmax[2])
				continue;

			countTotal++;

			if (rm0 > low_limit[0] && rm0 < high_limit[0] && rm1 > low_limit[1]
					&& rm1 < high_limit[1] && rm2 > low_limit[2] && rm2
					< high_limit[2])
				continue;

			if (currentMolecule->numSites() != currentMolecule->numLJcenters()) {
				std::cout
						<< "Molecule consists of something other than LJ centers. In this case RDF cannot be used.";
				return;
			}
			int boundary[3] = { 0, 0, 0 };
			if (r0 < low_limit[0])
				boundary[0] = -1;
			else if (r0 >= high_limit[0])
				boundary[0] = 1;
			if (r1 < low_limit[1])
				boundary[1] = -1;
			else if (r1 >= high_limit[1])
				boundary[1] = 1;
			if (r2 < low_limit[2])
				boundary[2] = -1;
			else if (r2 >= high_limit[2])
				boundary[2] = 1;

			//if (boundary[0] == 0 && boundary[1] == 0 && boundary[2] == 0)
			if (boundary[0] != -1)
				continue;

			// if molecule close to the boundary, add RDF force

			countBoundary++;

			// integration limits for axes
			double xlim[2] = { 0, 0 };
			double ylim[2] = { 0, 0 };
			double zlim[2] = { 0, 0 };
			double normal_lim[2] = { 5.6925, rc };
			// box size for the numerical quadrature
			double dx, dy, dz, dn, dr;
			dx = dy = dz = dn = dr = currentMolecule->getSigma() / 20;//1.0;//rc / 10;

			//double* currSum = new double[3];

			// yz plane is the boundary
			// integrate yz
			//for (unsigned int i = 0; boundary[0] != 0 && i
			//	< componentIds.size(); i++) {

			if (boundary[0] == -1) {
				xlim[1] = (r0 > rmin[0] ? rmin[0] : 2 * r0 - rmin[0]);
				xlim[0] = r0 - rc;
				normal_lim[0] = abs(r0 - rmin[0]);
				if (r0 < 0)
					normal_lim[0] = abs(2 * r0 - rmin[0]);
			} else if (boundary[0] == 1) {
				xlim[0] = (r0 < rmax[0] ? rmax[0] : 2 * r0 - rmax[0]);
				xlim[1] = r0 + rc;
				normal_lim[0] = abs(rmax[0] - r0);
			}
			//if (normal_lim[0] < 0.999) normal_lim[0] = 0.999;
			double h = xlim[1] - xlim[0];
			double a = std::sqrt(h * (2 * rc - h));

			ylim[0] = r1 - a;
			ylim[1] = r1 + a;
			zlim[0] = r2 - a;
			zlim[1] = r2 + a;

			integrateRDFSite(normal_lim, currentMolecule, rc, dr, dn,
					(*globalADist)[0], (*globalSiteADist)[0], 0, site, boundary);
			//integrateRDFCartesian(xlim, ylim, zlim, currentMolecule, rc, dx,
			//		dy, dz, (*globalADist)[0], (*globalSiteADist)[0], 0,
			//		site, boundary);

			//}

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

double BlockTraverse::integrateRDFCartesian(double xlim[2], double ylim[2],
		double zlim[2], Molecule* mol, double rc, double dx, double dy,
		double dz, vector<double> globalAcc,
		vector<vector<double> > globalSiteAcc, int plane, unsigned int site,
		int boundary[3]) {

	// volume of the domain
	double V = (_moleculeContainer->getBoundingBoxMax(0)
			- _moleculeContainer->getBoundingBoxMin(0))
			* (_moleculeContainer->getBoundingBoxMax(1)
					- _moleculeContainer->getBoundingBoxMin(1))
			* (_moleculeContainer->getBoundingBoxMax(2)
					- _moleculeContainer->getBoundingBoxMin(2));

	// number density of the domain
	//int numMolecules = _moleculeContainer->getNumberOfParticles();
	int numMolecules = 18522;
	double rho = numMolecules / (V);

	// molecule position
	double molr[3] = { mol->r(0) + mol->site_d(site)[0], mol->r(1)
			+ mol->site_d(site)[1], mol->r(2) + +mol->site_d(site)[2] };

	unsigned int other_site = (site == 0 ? 1 : 0);

	// potential
	double pot = 0;

	int bin, other_bin; // bin for the radius that rdf is read for
	double currf[3], absr, normal, radial, g = 0, currPot = 0;

	// dividing part of the sphere outside the bounding box into cells of size
	// dx, dy, dz
	for (double x = xlim[0]; x + dx / 2 <= xlim[1]; x += dx) {
		for (double y = ylim[0]; y + dy / 2 <= ylim[1]; y += dy) {
			for (double z = zlim[0]; z + dz / 2 <= zlim[1]; z += dz) {
				// distance of cell center to molecule
				double r[3] = { molr[0] - x - dx / 2, molr[1] - y - dy / 2,
						molr[2] - z - dz / 2 };

				double other_r[3] = { mol->r(0) + mol->site_d(other_site)[0]
						- x - dx / 2, mol->r(1) + mol->site_d(other_site)[1]
						- y - dy / 2, mol->r(2) + mol->site_d(other_site)[2]
						- z - dz / 2 };
				absr = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
				// check if the middle of the cell is within the cutoff radius,
				// if not, this cell will not contribute
				if (absr > rc)
					continue;

				if (mol->numSites() == 1) {
					bin = (int) (absr * globalAcc.size() / rc - 0.5);

					g = globalAcc[bin];

					PotForceLJ(r, absr * absr, 24 * mol->getEps(),
							mol->getSigma() * mol->getSigma(), currf, currPot);

					currPot *= 2 * PI * rho * g * radial * dx * dy * dz / 6;
					double f[3] = { 0, 0, 0 };

					for (int d = 0; d < 3; d++) {
						f[d] = rho * g * currf[d] * dx * dy * dz;
					}

					//mol->Fljcenteradd(site, f);
					if (boundary[0] == -1 && plane == 0) {
						mol->addLeftxRdfInfluence(site, f);
					}

					pot += currPot;
				} else {
					// if multiple LJ centers, use site-site rdf
					// iterate through sites, treat cell center as a site

					for (unsigned int j = 0; j < mol->numLJcenters(); j++) {

						// rdf (probability) value
						bin
								= (int) (absr * globalSiteAcc[site
										* mol->numLJcenters() + site].size()
										/ rc - 0.5);
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

						PotForceLJ(r, absr * absr, 24 * mol->getEps(),
								mol->getSigma() * mol->getSigma(), currf,
								currPot);

						//currPot *= 2 * PI * rho * g * radial * dx * dy * dz / 6;
						double f[3] = { 0, 0, 0 };

						for (int d = 0; d < 3; d++) {
							f[d] = rho * g * currf[d] * dx * dy * dz;
						}
						//if (f[0] > 0 && mol->id() == 18) cout<<"site: "<<site<<" r: "<<absr<<" value: "<<f[0]<<endl;

						//mol->Fljcenteradd(site, f);
						if (boundary[0] == -1 && plane == 0) {
							mol->addLeftxRdfInfluence(site, f);
						}
						//pot += currPot;
					}
				}
			}
		}

	}

	return pot;
}

void BlockTraverse::traversePairs(ParticlePairsHandler* particlePairsHandler,
		std::vector<std::string> rdf_file_names, int simstep, vector<vector<
				double> >* globalADist,
		vector<vector<vector<double> > >* globalSiteADist) {

	double _cutoffRadius = _moleculeContainer->getCutoff();
	double _LJCutoffRadius = _moleculeContainer->getLJCutoff();
	double _tersoffCutoffRadius = _moleculeContainer->getTersoffCutoff();
	vector<vector<unsigned long> >& forwardNeighbourOffsets =
			*_forwardNeighbourOffsets;
	vector<vector<unsigned long> >& backwardNeighbourOffsets =
			*_backwardNeighbourOffsets;

	particlePairsHandler->init();

	// XXX comment
	double distanceVector[3];
	// loop over all cells
	vector<Molecule*>::iterator molIter1;
	vector<Molecule*>::iterator molIter2;

#ifndef NDEBUG
	// reset forces and momenta to zero
	global_log->debug()
			<< "Resetting forces and momenta, disconnecting Tersoff pairs."
			<< endl;
#endif
	{
		double zeroVec[3] = { 0.0, 0.0, 0.0 };

		// TODO: check if the reset is done twice as leaving this part has no difference on the result.
		Molecule *moleculePtr;
		for (moleculePtr = _moleculeContainer->begin(); moleculePtr
				!= _moleculeContainer->end(); moleculePtr
				= _moleculeContainer->next()) {
			Molecule& molecule1 = *moleculePtr;
			molecule1.setF(zeroVec);
			molecule1.setM(zeroVec);
			molecule1.clearTersoffNeighbourList();
		}
	}

	vector<unsigned long>::iterator cellIndexIter;
	vector<unsigned long>::iterator neighbourOffsetsIter;

	// sqare of the cutoff radius
	double cutoffRadiusSquare = _cutoffRadius * _cutoffRadius;
	double LJCutoffRadiusSquare = _LJCutoffRadius * _LJCutoffRadius;
	double tersoffCutoffRadiusSquare = _tersoffCutoffRadius
			* _tersoffCutoffRadius;

#ifndef NDEBUG
	global_log->debug() << "Processing pairs and preprocessing Tersoff pairs."
			<< endl;
#endif

	// loop over all inner cells and calculate forces to forward neighbours
	for (cellIndexIter = _innerCellIndices.begin(); cellIndexIter
			!= _innerCellIndices.end(); cellIndexIter++) {

		unsigned long cellIndex = *cellIndexIter;

		ParticleCell& currentCell = _cells[cellIndex];
		if (currentCell.getMoleculeCount() < 1)
			continue;

		// forces between molecules in the cell
		processCell(currentCell, cutoffRadiusSquare, LJCutoffRadiusSquare,
				tersoffCutoffRadiusSquare, particlePairsHandler);

		// loop over all neighbours
		for (neighbourOffsetsIter = forwardNeighbourOffsets[cellIndex].begin(); neighbourOffsetsIter
				!= forwardNeighbourOffsets[cellIndex].end(); neighbourOffsetsIter++) {
			ParticleCell& neighbourCell = _cells[cellIndex
					+ *neighbourOffsetsIter];
			if (neighbourCell.getMoleculeCount() > 0) {
				processCellPair(currentCell, neighbourCell, cutoffRadiusSquare,
						LJCutoffRadiusSquare, tersoffCutoffRadiusSquare,
						particlePairsHandler);
			}
		}
	}

	// loop over halo cells and detect Tersoff neighbours within the halo
	// this is relevant for the angle summation
	for (cellIndexIter = _haloCellIndices.begin(); cellIndexIter
			!= _haloCellIndices.end(); cellIndexIter++) {
		unsigned long cellIndex = *cellIndexIter;
		ParticleCell& currentCell = _cells[cellIndex];
		if (currentCell.getMoleculeCount() < 1)
			continue;

		for (molIter1 = currentCell.getParticlePointers().begin(); molIter1
				!= currentCell.getParticlePointers().end(); molIter1++) {
			Molecule& molecule1 = **molIter1;
			if (molecule1.numTersoff() < 1)
				continue;
			molIter2 = molIter1;
			molIter2++;
			for (; molIter2 != currentCell.getParticlePointers().end(); molIter2++) {
				Molecule& molecule2 = **molIter2;
				assert(&molecule1 != &molecule2);
				if (molecule2.numTersoff() > 0) {
					double dd = molecule2.dist2(molecule1, distanceVector);
					if (dd < tersoffCutoffRadiusSquare)
						particlePairsHandler->preprocessTersoffPair(molecule1,
								molecule2, true);
				}
			}

			for (neighbourOffsetsIter
					= forwardNeighbourOffsets[cellIndex].begin(); neighbourOffsetsIter
					!= forwardNeighbourOffsets[cellIndex].end(); neighbourOffsetsIter++) {
				int j = cellIndex + *neighbourOffsetsIter;
				if ((j < 0) || (j >= (int) (_cells.size())))
					continue;
				ParticleCell& neighbourCell = _cells[j];
				if (!neighbourCell.isHaloCell())
					continue;
				for (molIter2 = neighbourCell.getParticlePointers().begin(); molIter2
						!= neighbourCell.getParticlePointers().end(); molIter2++) {
					Molecule& molecule2 = **molIter2;
					if (molecule2.numTersoff() < 1)
						continue;
					double dd = molecule2.dist2(molecule1, distanceVector);
					if (dd < tersoffCutoffRadiusSquare)
						particlePairsHandler->preprocessTersoffPair(molecule1,
								molecule2, true);
				}
			}
		}
	}

	// loop over all boundary cells and calculate forces to forward and backward neighbours
	for (cellIndexIter = _boundaryCellIndices.begin(); cellIndexIter
			!= _boundaryCellIndices.end(); cellIndexIter++) {
		unsigned long cellIndex = *cellIndexIter;
		ParticleCell& currentCell = _cells[cellIndex];

		if (currentCell.getMoleculeCount() < 1)
			continue;

		// forces between molecules in the cell
		processCell(currentCell, cutoffRadiusSquare, LJCutoffRadiusSquare,
				tersoffCutoffRadiusSquare, particlePairsHandler);

		// loop over all forward neighbours
		for (neighbourOffsetsIter = forwardNeighbourOffsets[cellIndex].begin(); neighbourOffsetsIter
				!= forwardNeighbourOffsets[cellIndex].end(); neighbourOffsetsIter++) {
			//cout<<"here!!!!!!!!!!!!!!!!!!"<<endl;
			//cout<<"offset"<< *neighbourOffsetsIter<<"idx "<<cellIndex<<endl;
			ParticleCell& neighbourCell = _cells[cellIndex
					+ *neighbourOffsetsIter];

			if (neighbourCell.getMoleculeCount() < 1)
				continue;

			// loop over all particles in the cell
			for (molIter1 = currentCell.getParticlePointers().begin(); molIter1
					!= currentCell.getParticlePointers().end(); molIter1++) {
				Molecule& molecule1 = **molIter1;
				unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

				for (molIter2 = neighbourCell.getParticlePointers().begin(); molIter2
						!= neighbourCell.getParticlePointers().end(); molIter2++) {
					Molecule& molecule2 = **molIter2;

					double dd = molecule2.dist2(molecule1, distanceVector);
					if (dd < cutoffRadiusSquare) {
						PairType pairType = MOLECULE_MOLECULE;
						if (neighbourCell.isHaloCell()
								&& !molecule1.isLessThan(molecule2)) {
							/* Do not sum up values twice. */
							pairType = MOLECULE_HALOMOLECULE;
						}
						double force[3] = { 0, 0, 0 };
						particlePairsHandler->processPair(molecule1, molecule2,
								distanceVector, pairType, dd, (dd
										< LJCutoffRadiusSquare), force, simstep);
						if ((num_tersoff > 0) && (molecule2.numTersoff() > 0)
								&& (dd < tersoffCutoffRadiusSquare)) {
							particlePairsHandler->preprocessTersoffPair(
									molecule1, molecule2, (pairType
											== MOLECULE_HALOMOLECULE));
						}
					}
				}
			}
		}

		// loop over all backward neighbours. calculate only forces
		// to neighbour cells in the halo region, all others already have been calculated
		for (neighbourOffsetsIter = backwardNeighbourOffsets[cellIndex].begin(); neighbourOffsetsIter
				!= backwardNeighbourOffsets[cellIndex].end(); neighbourOffsetsIter++) {

			ParticleCell& neighbourCell = _cells[cellIndex
					+ *neighbourOffsetsIter];

			if (neighbourCell.isHaloCell() && neighbourCell.getMoleculeCount()
					> 0) {
				// loop over all particles in the cell
				for (molIter1 = currentCell.getParticlePointers().begin(); molIter1
						!= currentCell.getParticlePointers().end(); molIter1++) {
					Molecule& molecule1 = **molIter1;
					unsigned int num_tersoff = molecule1.numTersoff(); // important for loop unswitching

					for (molIter2 = neighbourCell.getParticlePointers().begin(); molIter2
							!= neighbourCell.getParticlePointers().end(); molIter2++) {
						Molecule& molecule2 = **molIter2;
						double force[3] = { 0, 0, 0 };
						double dd = molecule2.dist2(molecule1, distanceVector);
						if (dd < cutoffRadiusSquare) {
							PairType
									pairType =
											molecule1.isLessThan(molecule2) ? MOLECULE_MOLECULE
													: MOLECULE_HALOMOLECULE;
							particlePairsHandler->processPair(molecule1,
									molecule2, distanceVector, pairType, dd,
									(dd < LJCutoffRadiusSquare), force, simstep);
							if ((num_tersoff > 0) && (molecule2.numTersoff()
									> 0) && (dd < tersoffCutoffRadiusSquare)) {
								particlePairsHandler->preprocessTersoffPair(
										molecule1, molecule2, (pairType
												== MOLECULE_HALOMOLECULE));
							}
						}
					}
				}
			}
		}

	}

#ifndef NDEBUG
	global_log->debug() << "processing Tersoff potential." << endl;
#endif
	double params[15];
	double delta_r = 0.;
	bool knowparams = false;

	for (cellIndexIter = _innerCellIndices.begin(); cellIndexIter
			!= _boundaryCellIndices.end(); cellIndexIter++) {

		if (cellIndexIter == _innerCellIndices.end())
			cellIndexIter = _boundaryCellIndices.begin();

		ParticleCell& currentCell = _cells[*cellIndexIter];
		for (molIter1 = currentCell.getParticlePointers().begin(); molIter1
				!= currentCell.getParticlePointers().end(); molIter1++) {
			Molecule& molecule1 = **molIter1;

			if (molecule1.numTersoff() < 1)
				continue;

			if (!knowparams) {
				delta_r = molecule1.tersoffParameters(params);
				knowparams = true;
			}
			particlePairsHandler->processTersoffAtom(molecule1, params, delta_r);
		}
	}

	if (rdf_file_names.size() > 0) {
		this->traverseRDFBoundaryCartesian(globalADist, globalSiteADist,
				particlePairsHandler);

	}
	particlePairsHandler->finish();
}
