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
#include <string>

#include <iostream>
using namespace std;


class DomainDecompBase;
class ParticleContainer;

class PressureGradient {
public:
	PressureGradient(int rank);

	//! @brief assigns a coset ID to a component (ID)
	void assignCoset(unsigned int cid, unsigned int cosetid) { _universalComponentSetID[cid] = cosetid; }
	//! @brief sets the information on the acceleration model for one coset
	void specifyComponentSet(unsigned int cosetid, double v[3], double tau, double ainit[3], double timestep);
	//! @brief sets type of movement for each component (ID)
	void specifyMovementStyle(unsigned int cid, std::string moveStyle) { this->_movementCompID[cid] = moveStyle; }
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
	bool isAcceleratingInstantaneously(unsigned numberOfComp);
	
	//! @brief updates the intensity and direction of the uniform acceleration
	void determineAdditionalAcceleration(
		DomainDecompBase* domainDecomp,
		ParticleContainer* molCont, double dtConstantAcc
	);
	//! @brief returns the type of movement of each component (ID)
	std::string getMovementStyle(unsigned int cid){ return this->_movementCompID[cid]; }
	//! @brief returns the acceleration map (necessary for passing data to the integrator)
	std::map<unsigned int, double>* getUAA() { return this->_universalAdditionalAcceleration; }
	//! @brief returns the cid of a component due to its type of movement
	unsigned getCidMovement(std::string moveStyle, unsigned numberOfComp){
	      for(unsigned i = 1; i < (numberOfComp+1); i++){
		if(getMovementStyle(i) == moveStyle){
		    return i;
		}
	      }
	      return 0;
	} 
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
	//! @brief calculates globalN and globalVelocitySum
	void prepare_getMissingVelocity(DomainDecompBase* domainDecomp, ParticleContainer* molCont, unsigned int cosetid, unsigned numberOfComp);
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
	
	//TEST: Spring Force Influence
	void specifySpringInfluence(unsigned long minSpringID, unsigned long maxSpringID, double averageYPos, double springConst);
	bool isSpringDamped() {
		return (this->_springConst > 0)
				&& (this->_maxSpringID > this->_minSpringID);
	}
	unsigned long getMinSpringID() { return this->_minSpringID; }
	unsigned long getMaxSpringID() { return this->_maxSpringID; }
	double getAverageY() { return this->_averageYPos; }
	double getSpringConst() { return this->_springConst; }
	void setGlobalVelSumBeforeAcc(int d, unsigned cosetid, long double VelSum) {_globalVelSumBeforeAcc[d][cosetid] = VelSum; }
	void setGlobalVelSumAfterAcc(int d, unsigned cosetid, long double VelSum) {_globalVelSumAfterAcc[d][cosetid] = VelSum; }
	void setGlobalVelSumAfterThT(int d, unsigned cosetid, long double VelSum) {_globalVelSumAfterThT[d][cosetid] = VelSum; }
	double getGlobalVelSumBeforeAcc(int d, int cosetid) {return this->_globalVelSumBeforeAcc[d][cosetid]; }
	double getGlobalVelSumAfterAcc(int d, int cosetid) {return this->_globalVelSumAfterAcc[d][cosetid]; }
	double getGlobalVelSumAfterThT(int d, int cosetid) {return this->_globalVelSumAfterThT[d][cosetid]; }
	double getGlobalVelSum(int d, unsigned cosetid){return this->_globalVelocitySum[d][cosetid]; }
	unsigned long getGlobalN(int cosetid){return this->_globalN[cosetid]; }
	double getGlobalForceSum(int d, int cosetid){return this->_globalForceSum[d][cosetid]; }
	double getGlobalSpringForceSum(int d, int cosetid){return this->_globalSpringForceSum[d][cosetid]; }
	//! @brief collects the total force vector of each component "cosetid" as well as the spring force vector
	void calculateForcesOnComponent(ParticleContainer* molCont, unsigned int cid);
	void collectForcesOnComponent(DomainDecompBase* domainDecomp, ParticleContainer* molCont, unsigned int cid);
	void resetForcesOnComponent(unsigned int cid);
	void calculateSpringForcesOnComponent(ParticleContainer* molCont, unsigned int cid);
	void collectSpringForcesOnComponent(DomainDecompBase* domainDecomp, ParticleContainer* molCont, unsigned int cid);
	void resetSpringForcesOnComponent(unsigned int cid);
	//! @brief counts the call of the method pressureGradient->collectSpringForcesOnComponent() to ensure that its just called once each timestep */
	void setCounter(unsigned long counter){_counter = counter; }
	unsigned long getCounter(){return _counter; }
	
	// for the measurements in the confinement
	double getGlobalTargetVelocity(int d, int coset) { return this->_globalTargetVelocity[d][coset]; }

private:
	unsigned int _localRank;

	/// calculate new value of the uniform acceleration each # timesteps
	unsigned int _universalConstantAccelerationTimesteps;
	/// assigns a component set ID to some of the components
	std::map<unsigned int, unsigned int> _universalComponentSetID;
	/// assigns the type of movement to the components
	std::map<unsigned, std::string> _movementCompID;
	/// local number of molecules that belong to a given component set ID
	std::map<unsigned int, unsigned long> _localN;
	/// global number of molecules that belong to a given component set ID
	std::map<unsigned int, unsigned long> _globalN;
	/// local number of molecules that belong to a given component set ID
	std::map<unsigned int, unsigned long> _localSpringN;
	/// global number of molecules that belong to a given component set ID
	std::map<unsigned int, unsigned long> _globalSpringN;
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
	
	//TEST: Spring Force Influence
	unsigned long _minSpringID, _maxSpringID;
	double _averageYPos, _springConst;
	/// global sum of the velocity vectors before artificial acceleration corresponding to a given component set ID
	std::map<unsigned int, long double> _globalVelSumBeforeAcc[3];
	/// global sum of the velocity vectors after artificial acceleration corresponding to a given component set ID
	std::map<unsigned int, long double> _globalVelSumAfterAcc[3];
	/// global sum of the velocity vectors after artificial acceleration and after thermostat-scaling corresponding to a given component set ID
	std::map<unsigned int, long double> _globalVelSumAfterThT[3];
	/// local sum of the force vectors corresponding to a given component set ID
	std::map<unsigned int, long double> _localForceSum[3];
	/// global sum of the force vectors corresponding to a given component set ID
	std::map<unsigned int, long double> _globalForceSum[3];
	/// local sum of the spring force vectors corresponding to a given component set ID
	std::map<unsigned int, long double> _localSpringForceSum[3];
	/// global sum of the spring force vectors corresponding to a given component set ID
	std::map<unsigned int, long double> _globalSpringForceSum[3];
	/// counts the call of the method pressureGradient->collectSpringForcesOnComponent()
	unsigned long _counter;
	
};	
#endif /* PRESSUREGRADIENT_H_ */
