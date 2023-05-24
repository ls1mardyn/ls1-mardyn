
#ifndef PRESSUREGRADIENT_H_
#define PRESSUREGRADIENT_H_

#include <deque>
#include <map>

class DomainDecompBase;
class ParticleContainer;

/*! @brief DEPRECATED CLASS
 *
 * PressureGradient functionality was lost at the latest with the move to the .XML input format.
 * Code includes and usage was removed for code clarity and performance with a pull request by Kruegener on 01/22/2019.
 * The removed code was attached to PressureGradient.h for reference.
 * If the functionality would be needed, it should be implemented via the standard PluginBase interfaces, as the new
 * plugin system is powerful enough to implement this functionality without direct intervention all over the code base.
 */
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
#endif /* PRESSUREGRADIENT_H_ */

/* PROCEDURE TO REMOVE UNUSED PRESSURE GRADIENT FROM CODEBASE:

 	####### REMOVED PG FROM DOMAIN CONSTRUCTOR -> REMOVED UNIVERSALPG -> REMOVED FORWARD DECL IN .H AND INCLUDE IN .CPP

		this->_universalPG = pg;
		// after: checkpointfilestream << _epsilonRF << std::endl;
				std::map<unsigned, unsigned> componentSets = this->_universalPG->getComponentSets();
				for( std::map<unsigned, unsigned>::const_iterator uCSIDit = componentSets.begin();
						uCSIDit != componentSets.end();
						uCSIDit++ )
				{
					if(uCSIDit->first > 100) continue;
					checkpointfilestream << " S\t" << 1+uCSIDit->first << "\t" << uCSIDit->second << "\n";
				}
				std::map<unsigned, double> tau = this->_universalPG->getTau();
				for( std::map<unsigned, double>::const_iterator gTit = tau.begin();
						gTit != tau.end();
						gTit++ )
				{
					unsigned cosetid = gTit->first;
					double* ttargetv = this->_universalPG->getTargetVelocity(cosetid);
					double* tacc = this->_universalPG->getAdditionalAcceleration(cosetid);
					checkpointfilestream << " A\t" << cosetid << "\t"
						<< ttargetv[0] << " " << ttargetv[1] << " " << ttargetv[2] << "\t"
						<< gTit->second << "\t"
						<< tacc[0] << " " << tacc[1] << " " << tacc[2] << "\n";
					delete ttargetv;
					delete tacc;
				}

		// Function in Domain
			PressureGradient* getPG() { return this->_universalPG; }

	####### REMOVED FORWARD DECLARATION AND INCLUDE FROM SIMULATION.H/.CPP

		PressureGradient* _pressureGradient;
 		//after: 	_longRangeCorrection->calculateLongRange(); in #######
			if (_pressureGradient->isAcceleratingUniformly()) {
			global_log->info() << "Initialising uniform acceleration." << std::endl;
			unsigned long uCAT = _pressureGradient->getUCAT();
			global_log->info() << "uCAT: " << uCAT << " steps." << std::endl;
			_pressureGradient->determineAdditionalAcceleration(
					_domainDecomposition, _moleculeContainer, uCAT
							* _integrator->getTimestepLength());
			global_log->info() << "Uniform acceleration initialised." << std::endl;
			}

		// first in simulate()
			// (universal) constant acceleration (number of) timesteps
			unsigned uCAT = _pressureGradient->getUCAT();

		// after: _domain->calculateThermostatDirectedVelocity(_moleculeContainer); in simulate()
			if (_pressureGradient->isAcceleratingUniformly()) {
				if (!(_simstep % uCAT)) {
					global_log->debug() << "Determine the additional acceleration" << std::endl;
					_pressureGradient->determineAdditionalAcceleration(
							_domainDecomposition, _moleculeContainer, uCAT
							* _integrator->getTimestepLength());
				}
				global_log->debug() << "Process the uniform acceleration" << std::endl;
				_integrator->accelerateUniformly(_moleculeContainer, _domain);
				_pressureGradient->adjustTau(this->_integrator->getTimestepLength());
			}

		// in initialize() before Domain()
			global_log->info() << "Creating PressureGradient ... " << std::endl;
			_pressureGradient = new PressureGradient(ownrank);

 	####### REMOVED FUNCTION ONLY CALLED BY PG FROM INTEGRATOR, LEAPFROG AND LEAPFROGRMM

	// Integrator
 		virtual void accelerateUniformly(
			ParticleContainer* molCont,
			Domain* domain
		) = 0;

		virtual void accelerateInstantaneously(
			ParticleContainer* molCont,
			Domain* domain
		) = 0;

 	// LeapfrogRMM
 		void accelerateUniformly(
			ParticleContainer* molCont,
			Domain* domain
		) {}

		void accelerateInstantaneously(
			ParticleContainer* molCont,
			Domain* domain
		) {}

	// Leapfrog
 		virtual void accelerateUniformly(
			ParticleContainer* molCont,
			Domain* domain
		);

		virtual void accelerateInstantaneously(
			ParticleContainer* molCont,
			Domain* domain
		);

		void Leapfrog::accelerateUniformly(ParticleContainer* molCont, Domain* domain) {
			std::map<unsigned, double>* additionalAcceleration = domain->getPG()->getUAA();
			std::vector<Component> comp = *(_simulation.getEnsemble()->getComponents());
			std::vector<Component>::iterator compit;
			std::map<unsigned, double> componentwiseVelocityDelta[3];
			for (compit = comp.begin(); compit != comp.end(); compit++) {
				unsigned cosetid = domain->getPG()->getComponentSet(compit->ID());
				if (cosetid != 0)
					for (unsigned d = 0; d < 3; d++)
						componentwiseVelocityDelta[d][compit->ID()] = _timestepLength * additionalAcceleration[d][cosetid];
				else
					for (unsigned d = 0; d < 3; d++)
						componentwiseVelocityDelta[d][compit->ID()] = 0;
			}

			#if defined(_OPENMP)
			#pragma omp parallel
			#endif
			{
				for (auto thismol = molCont->iterator(); thismol.isValid(); ++thismol) {
					unsigned cid = thismol->componentid();
					mardyn_assert(componentwiseVelocityDelta[0].find(cid) != componentwiseVelocityDelta[0].end());
					thismol->vadd(componentwiseVelocityDelta[0][cid],
								  componentwiseVelocityDelta[1][cid],
								  componentwiseVelocityDelta[2][cid]);
				}
			}
		}

		void Leapfrog::accelerateInstantaneously(ParticleContainer* molCont, Domain* domain) {
			std::vector<Component> comp = *(_simulation.getEnsemble()->getComponents());
			std::vector<Component>::iterator compit;
			std::map<unsigned, double> componentwiseVelocityDelta[3];
			for (compit = comp.begin(); compit != comp.end(); compit++) {
				unsigned cosetid = domain->getPG()->getComponentSet(compit->ID());
				if (cosetid != 0)
					for (unsigned d = 0; d < 3; d++)
						componentwiseVelocityDelta[d][compit->ID()] = domain->getPG()->getMissingVelocity(cosetid, d);
				else
					for (unsigned d = 0; d < 3; d++)
						componentwiseVelocityDelta[d][compit->ID()] = 0;
			}

			#if defined(_OPENMP)
			#pragma omp parallel
			#endif
			{

				for (auto thismol = molCont->iterator(); thismol.isValid(); ++thismol) {
					unsigned cid = thismol->componentid();
					mardyn_assert(componentwiseVelocityDelta[0].find(cid) != componentwiseVelocityDelta[0].end());
					thismol->vadd(componentwiseVelocityDelta[0][cid],
								  componentwiseVelocityDelta[1][cid],
								  componentwiseVelocityDelta[2][cid]);
				}
			}
		}

 	##### removed getPG() related functions

 	// ASCIIReader.cpp after: else if((token == "NumberOfMolecules") || (token == "N")) {

 		else if((token == "AssignCoset") || (token == "S")) {
			unsigned component_id, cosetid;
			_phaseSpaceHeaderFileStream >> component_id >> cosetid;
			component_id--; // FIXME: Component ID starting with 0 in program ...
			domain->getPG()->assignCoset( component_id, cosetid );
		}
		else if((token == "Accelerate") || (token == "A")) {
			unsigned cosetid;
			_phaseSpaceHeaderFileStream >> cosetid;
			double v[3];
			for(unsigned d = 0; d < 3; d++)
				_phaseSpaceHeaderFileStream >> v[d];
			double tau;
			_phaseSpaceHeaderFileStream >> tau;
			double ainit[3];
			for(unsigned d = 0; d < 3; d++)
				_phaseSpaceHeaderFileStream >> ainit[d];
			domain->getPG()->specifyComponentSet(cosetid, v, tau, ainit, timestep);
		}

 	// MPI_IOReader after: else if((token == "NumberOfMolecules") || (token == "N")) {

 		else if((token == "AssignCoset") || (token == "S")) {
			unsigned component_id, cosetid;
			_phaseSpaceHeaderFileStream >> component_id >> cosetid;
			component_id--; // FIXME: Component ID starting with 0 in program ...
			domain->getPG()->assignCoset( component_id, cosetid );
		}
		else if((token == "Accelerate") || (token == "A")) {
			unsigned cosetid;
			_phaseSpaceHeaderFileStream >> cosetid;
			double v[3];
			for(unsigned d = 0; d < 3; d++)
				_phaseSpaceHeaderFileStream >> v[d];
			double tau;
			_phaseSpaceHeaderFileStream >> tau;
			double ainit[3];
			for(unsigned d = 0; d < 3; d++)
				_phaseSpaceHeaderFileStream >> ainit[d];
			domain->getPG()->specifyComponentSet(cosetid, v, tau, ainit, timestep);
		}
	##### Removed unused includes
 		Simulation
 		Leapfrog
 		LeapfrogRMM
 		ASCIIReader
 		MPI_IOReader
 		BinaryReader
 		LinkedCellsTest
		MDGenerator
			- inhereting generators
*/