//
// Created by Kruegener
//

#include "GrandCanonicalEnsemble.h"
#include "DomainBase.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "ChemicalPotential.h"

GrandCanonicalEnsemble::GrandCanonicalEnsemble() :
		_N(0), _V(0), _T(0), _mu(0), _p(0), _E(0), _E_trans(0), _E_rot(0) {
	_type = muVT;
	_simulationDomain = global_simulation->getDomain();
}

void GrandCanonicalEnsemble::initConfigXML(ParticleContainer* moleculeContainer) {
	int ownrank = 0;
	unsigned long globalNumMolecules = _simulationDomain->getglobalNumMolecules(true, moleculeContainer, nullptr);
#ifdef ENABLE_MPI
	MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &ownrank) );
#endif

	unsigned idi = _lmu.size();
	unsigned j = 0;
	std::list<ChemicalPotential>::iterator cpit;
	for(cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
		cpit->setIncrement(idi);
		double tmp_molecularMass = global_simulation->getEnsemble()->getComponent(cpit->getComponentID())->m();
		cpit->setSystem(_simulationDomain->getGlobalLength(0),
						_simulationDomain->getGlobalLength(1), _simulationDomain->getGlobalLength(2),
						tmp_molecularMass);
		cpit->setGlobalN(global_simulation->getEnsemble()->getComponent(cpit->getComponentID())->getNumMolecules());
		cpit->setNextID(j + (int) (1.001 * (256 + globalNumMolecules)));

		cpit->setSubdomain(ownrank, moleculeContainer->getBoundingBoxMin(0),
						   moleculeContainer->getBoundingBoxMax(0),
						   moleculeContainer->getBoundingBoxMin(1),
						   moleculeContainer->getBoundingBoxMax(1),
						   moleculeContainer->getBoundingBoxMin(2),
						   moleculeContainer->getBoundingBoxMax(2));
		/* TODO: thermostat */
		double Tcur = _simulationDomain->getCurrentTemperature(0);
		/* FIXME: target temperature from thermostat ID 0 or 1?  */
		double Ttar = _simulationDomain->severalThermostats() ? _simulationDomain->getTargetTemperature(1)
															  : _simulationDomain->getTargetTemperature(0);
		if((Tcur < 0.85 * Ttar) || (Tcur > 1.15 * Ttar))
			Tcur = Ttar;
		cpit->submitTemperature(Tcur);
		if(global_simulation->getH() != 0.0)
			cpit->setPlanckConstant(global_simulation->getH());
		j++;
	}
}

void GrandCanonicalEnsemble::prepare_start() {
	if(_lmu.size() > 0) {
		/* TODO: thermostat */
		double Tcur = _simulationDomain->getGlobalCurrentTemperature();
		/* FIXME: target temperature from thermostat ID 0 or 1? */
		double
				Ttar = _simulationDomain->severalThermostats() ? _simulationDomain->getTargetTemperature(1)
															   : _simulationDomain->getTargetTemperature(0);
		if((Tcur < 0.85 * Ttar) || (Tcur > 1.15 * Ttar))
			Tcur = Ttar;

		list<ChemicalPotential>::iterator cpit;
		if(global_simulation->getH() == 0.0)
			global_simulation->setH(sqrt(6.2831853 * Ttar));
		for(cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
			cpit->submitTemperature(Tcur);
			cpit->setPlanckConstant(global_simulation->getH());
		}
	}
}

void GrandCanonicalEnsemble::beforeEventNewTimestep(ParticleContainer* moleculeContainer,
													DomainDecompBase* domainDecomposition,
													unsigned long simstep) {
	/** @todo What is this good for? Where come the numbers from? Needs documentation */
	if(simstep >= _initGrandCanonical) {
		unsigned j = 0;
		list<ChemicalPotential>::iterator cpit;
		for(cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
			if(!((simstep + 2 * j + 3) % cpit->getInterval())) {
				cpit->prepareTimestep(moleculeContainer, domainDecomposition);
			}
			j++;
		}
	}
}

void GrandCanonicalEnsemble::afterForces(ParticleContainer* moleculeContainer, DomainDecompBase* domainDecomposition,
										 CellProcessor* cellProcessor,
										 unsigned long simstep) {

	if(simstep >= _initGrandCanonical) {
		unsigned j = 0;
		list<ChemicalPotential>::iterator cpit;
		for(cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
			if(!((simstep + 2 * j + 3) % cpit->getInterval())) {
				global_log->debug() << "Grand canonical ensemble(" << j << "): test deletions and insertions"
									<< endl;
				this->_simulationDomain->setLambda(cpit->getLambda());
				this->_simulationDomain->setDensityCoefficient(cpit->getDensityCoefficient());
				double localUpotBackup = _simulationDomain->getLocalUpot();
				double localVirialBackup = _simulationDomain->getLocalVirial();
				cpit->grandcanonicalStep(moleculeContainer, _simulationDomain->getGlobalCurrentTemperature(),
										 this->_simulationDomain,
										 cellProcessor);
				_simulationDomain->setLocalUpot(localUpotBackup);
				_simulationDomain->setLocalVirial(localVirialBackup);

#ifndef NDEBUG
				/* check if random numbers inserted are the same for all processes... */
				cpit->assertSynchronization(domainDecomposition);
#endif

				int localBalance = cpit->getLocalGrandcanonicalBalance();
				int balance = cpit->grandcanonicalBalance(domainDecomposition);
				global_log->debug() << "   b[" << ((balance > 0) ? "+" : "") << balance << "("
									<< ((localBalance > 0) ? "+" : "") << localBalance << ")" << " / c = "
									<< cpit->getComponentID() << "]   " << endl;
				_simulationDomain->Nadd(cpit->getComponentID(), balance, localBalance);
			}

			j++;
		}

		// TODO: Originally this comes after the deleteOuterParticles call

		_simulationDomain->evaluateRho(moleculeContainer->getNumberOfParticles(), domainDecomposition);


	}

}

void GrandCanonicalEnsemble::storeSample(Molecule* m, uint32_t componentid) {
	std::list<ChemicalPotential>::iterator cpit;
	for(cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
		if(!cpit->hasSample() && (componentid == cpit->getComponentID())) {
			cpit->storeMolecule(*m);
		}
	}
}