/*************************************************************************
 * Copyright (C) 2010 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
 *                                                                       *
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or (at *
 * your option) any later version.                                       *
 *                                                                       *
 * This program is distributed in the hope that it will be useful, but   *
 * WITHOUT ANY WARRANTY; without even the implied warranty of            * 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      *
 * General Public License for more details.                              *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program; if not, write to the Free Software           *
 * Foundation, 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.   *
 *************************************************************************/

#ifndef PRESSUREGRADIENT_H_
#define PRESSUREGRADIENT_H_

#include <deque>
#include <map>

class DomainDecompBase;
class ParticleContainer;

class PressureGradient {
public:
	PressureGradient(int rank);

	//! @brief assigns a coset ID to a component (ID)
	void assignCoset(unsigned int cid, unsigned int cosetid) { _universalComponentSetID[cid] = cosetid; }
	//! @brief sets the information on the acceleration model for one coset
	void specifyComponentSet(unsigned int cosetid, double v[3], double tau, double ainit[3], double timestep);
	//! @brief sets the number of timesteps between two updates of the uniform acceleration
	void setUCAT(unsigned int uCAT) { this->_universalConstantAccelerationTimesteps = uCAT; }
	//! @brief returns the number of timesteps between two updates of the uniform acceleration
	unsigned int getUCAT() { return this->_universalConstantAccelerationTimesteps; }
	/// sets the zeta value for the flow controller
	void setZetaFlow(double zeta) {
		this->_universalConstantTau = false;
		this->_universalZetaFlow = zeta;
	}
	void specifyTauPrime(double tauPrime, double dt);
	void adjustTau(double dt);
	//! @brief are there any cosets?
	bool isAcceleratingUniformly() {
		return (this->_universalTau.size() > 0)
				&& (this->_universalConstantAccelerationTimesteps > 0);
	}
	//! @brief updates the intensity and direction of the uniform acceleration
	void determineAdditionalAcceleration(
		DomainDecompBase* domainDecomp,
		ParticleContainer* molCont, double dtConstantAcc
	);
	//! @brief returns the acceleration map (necessary for passing data to the integrator)
	std::map<unsigned int, double>* getUAA() { return this->_universalAdditionalAcceleration; }
	//! @brief returns the cosetid of a component (0 for unaccelerated components)
	unsigned int getComponentSet(unsigned int cid) {
		if(_universalComponentSetID.find(cid) == _universalComponentSetID.end())
			return 0;
		else
			return this->_universalComponentSetID[cid];
	}
	//! @brief returns the directed velocity for a component set
	//! @param cosetid ID of the component set
	//! @param d x direction (0), y direction (1), or z direction (2)
	double getDirectedVelocity(unsigned int cosetid, unsigned short int d);
	//! @brief returns the absolute external acceleration for a component set
	//! @param cosetid ID of the component set
	double getUniformAcceleration(unsigned int cosetid);
	//! @brief returns the external acceleration for a component set
	//! @param cosetid ID of the component set
	//! @param d x direction (0), y direction (1), or z direction (2)
	double getUniformAcceleration(unsigned int cosetid, unsigned short int d);
	//! @brief returns the difference between the desired velocity and the global average velocity
	double getMissingVelocity(unsigned int cosetid, unsigned short int d);
	//! @brief total number of particles that belong to the specified component set
	double getCosetN(unsigned int cosetid) { return this->_globalN[cosetid]; }
	unsigned int maxCoset() { return this->_universalTau.size(); }

	//! @brief returns the component -> set ID map
	std::map<unsigned int, unsigned int> getComponentSets() { return this->_universalComponentSetID; }

	std::map<unsigned int, double> getTau() { return this->_universalTau; }
	double getTau(unsigned int set) { return this->_universalTau[set]; }

	double* getTargetVelocity(unsigned int set);
	double* getAdditionalAcceleration(unsigned int set);

private:
	unsigned int _localRank;

	/// calculate new value of the uniform acceleration each # timesteps
	unsigned int _universalConstantAccelerationTimesteps;
	/// assigns a component set ID to some of the components
	std::map<unsigned int, unsigned int> _universalComponentSetID;
	/// local number of molecules that belong to a given component set ID
	std::map<unsigned int, unsigned long> _localN;
	/// global number of molecules that belong to a given component set ID
	std::map<unsigned int, unsigned long> _globalN;
	/// local sum of the velocity vectors corresponding to a given component set ID
	std::map<unsigned int, long double> _localVelocitySum[3];
	/// global sum of the velocity vectors corresponding to a given component set ID
	std::map<unsigned int, long double> _globalVelocitySum[3];
	/// uniform acceleration
	std::map<unsigned int, double> _universalAdditionalAcceleration[3];
	/// target average velocity for the molecules of a coset
	std::map<unsigned int, double> _globalTargetVelocity[3];
	/// delay variable tau of a coset
	std::map<unsigned int, double> _universalTau;
	/// is the tau parameter constant
	bool _universalConstantTau;
	/// zeta parameter of the flow regulation
	double _universalZetaFlow;
	/// tau prime (t') parameter of the flow regulation
	double _universalTauPrime;
	/// queue of previously recorded velocity sums
	std::map<unsigned int, std::deque<long double> > _globalPriorVelocitySums[3];
	/// number of items in the velocity queue
	std::map<unsigned int, unsigned int> _globalVelocityQueuelength;
};
#endif

