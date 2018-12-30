//
// Created by Kruegener on 12/29/2018.
//

#include "Simulation.h"
#include "GrandCanonicalEnsemble.h"
#include "DomainBase.h"
#include "ChemicalPotential.h"
#include "Domain.h"

GrandCanonicalEnsemble::GrandCanonicalEnsemble() :
        _N(0), _V(0), _T(0), _mu(0), _p(0), _E(0), _E_trans(0), _E_rot(0) {
    _type = "muVT";
    _simulationDomain = global_simulation->getDomain();
}

void GrandCanonicalEnsemble::initConfigXML(ParticleContainer *moleculeContainer, double h) {
    int ownrank = 0;
#ifdef ENABLE_MPI
    MPI_CHECK( MPI_Comm_rank(MPI_COMM_WORLD, &ownrank) );
#endif

    unsigned idi = _lmu.size();
    unsigned j = 0;
    std::list<ChemicalPotential>::iterator cpit;
    for (cpit = _lmu.begin(); cpit != _lmu.end(); cpit++) {
        cpit->setIncrement(idi);
        double tmp_molecularMass = global_simulation->getEnsemble()->getComponent(cpit->getComponentID())->m();
        cpit->setSystem(_simulationDomain->getGlobalLength(0),
                        _simulationDomain->getGlobalLength(1), _simulationDomain->getGlobalLength(2),
                        tmp_molecularMass);
        cpit->setGlobalN(global_simulation->getEnsemble()->getComponent(cpit->getComponentID())->getNumMolecules());
        cpit->setNextID(j + (int) (1.001 * (256 + _simulationDomain->getglobalNumMolecules())));

        cpit->setSubdomain(ownrank, moleculeContainer->getBoundingBoxMin(0),
                           moleculeContainer->getBoundingBoxMax(0),
                           moleculeContainer->getBoundingBoxMin(1),
                           moleculeContainer->getBoundingBoxMax(1),
                           moleculeContainer->getBoundingBoxMin(2),
                           moleculeContainer->getBoundingBoxMax(2));
        /* TODO: thermostat */
        double Tcur = _simulationDomain->getCurrentTemperature(0);
        /* FIXME: target temperature from thermostat ID 0 or 1?  */
        double Ttar = _simulationDomain->severalThermostats() ? _simulationDomain->getTargetTemperature(1) : _simulationDomain->getTargetTemperature(0);
        if ((Tcur < 0.85 * Ttar) || (Tcur > 1.15 * Ttar))
            Tcur = Ttar;
        cpit->submitTemperature(Tcur);
        if (h != 0.0)
            cpit->setPlanckConstant(h);
        j++;
    }
}
