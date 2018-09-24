/**
 * @file AutoPasFullMolecule.h
 * @author seckler
 * @date 20.09.18
 */

#pragma once

#include <autopas/utils/SoAType.h>

#include "FullMolecule.h"


/**
 * class that implements additional functions to make the molecule compatible with autopas
 */
class AutoPasSimpleMolecule : public FullMolecule {
public:
	explicit AutoPasSimpleMolecule(unsigned long id = 0, Component *component = nullptr,
	                    double rx = 0., double ry = 0., double rz = 0.,
	                    double vx = 0., double vy = 0., double vz = 0.,
	                    double q0 = 1., double q1 = 1., double q2 = 0., double q3 = 0.,
	                    double Dx = 0., double Dy = 0., double Dz = 0.
	) : FullMolecule(id, component, rx, ry, rz, vx, vy, vz, q0, q1, q2, q3, Dx, Dy, Dz) {}

	AutoPasSimpleMolecule(const AutoPasSimpleMolecule &m) = default;

	/**
	 * get the position of the particle as std::array
	 * @return
	 */
	const std::array<double, 3> &getR() const;

	/**
	 * Get the force acting on the particle as std::array.
	 * @return
	 */
	const std::array<double, 3> &getF() const;

	/**
	 * Set the force acting on the particle as std::array.
	 * @param F force acting on the particle
	 */
	void setF(const std::array<double, 3>& F);

	/**
	 * Adds F to the force acting on the particle.
	 * @param F force to be added
	 */
	void addF(const std::array<double, 3>& F);

	/**
	 * Substracts F to the force acting on the particle.
	 * @param F force to be substracted
	 */
	void subF(const std::array<double, 3>& F);

	// somehow the definition in MoleculeInterface is not found, so we add it here, again...
	void setF(double F[3]) override;

	/**
	 * test whether the particle is in the box specified by rmin and rmax
	 * @param rmin lower corner of the box (inclusive)
	 * @param rmax upper corner of the box (exclusive)
	 * @return true if it is in the box, false otherwise
	 */
	bool inBox(const std::array<double, 3> &rmin, const std::array<double, 3>& rmax) const;

	// somehow the definition in MoleculeInterface is not found, so we add it here, again...
	bool inBox(const double rmin[3], const double rmax[3]) const override {
		return MoleculeInterface::inBox(rmin, rmax);
	}

	/**
	 * print info of molecule
	 * @return
	 */
	std::string toString();

	/**
	 * Enums used as ids for accessing and creating a dynamically sized SoA.
	 */
	enum AttributeNames : int { id, posX, posY, posZ, forceX, forceY, forceZ };

	/**
	 * typedef needed for soa's of autopas
	 */
	typedef autopas::utils::SoAType<size_t, double, double, double, double, double, double>::Type SoAArraysType;
};
