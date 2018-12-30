//
// Created by Kruegener
//

#ifndef MARDYN_GRANDCANONICAL_H
#define MARDYN_GRANDCANONICAL_H


#include "EnsembleBase.h"

class ChemicalPotential;
class Domain;
class DomainDecompBase;
class CellProcessor;

class GrandCanonicalEnsemble : public Ensemble{

public:
    GrandCanonicalEnsemble();

    virtual ~GrandCanonicalEnsemble() override {
    }

    // TODO: Implement STUB
    void readXML(XMLfileUnits& xmlconfig) override {
        global_log->info() << "[GrandCanonicalEnsemble] readXML not implemented!" << std::endl;
    };

    unsigned long N() override {
        return _N;
    }
    double V() override {
        return _V;
    }
    double T() override {
        return _T;
    }

    double mu() override {
        return _mu;
    }
    double p() override {
        return _p;
    }
    double E() override {
        return _E;
    }

    // TODO: Implement
    void updateGlobalVariable(ParticleContainer *particleContainer, GlobalVariable variable) override {};

    /*! Returns _lmu pointer for processing by external plugins */
    std::list<ChemicalPotential>* getLmu() override {return &_lmu;}

    /*! Runs steps formerly in initConfigXML in simulation.cpp */
    void initConfigXML(ParticleContainer *moleculeContainer) override;

    /*! Runs steps formerly in prepare_start in simulation.cpp */
    void prepare_start() override;

    /*! Runs steps formerly in simulate in simulation.cpp */
    void beforeEventNewTimestep(ParticleContainer *moleculeContainer, DomainDecompBase *domainDecomposition,
                                unsigned long simstep) override;

    /*! Runs steps formerly in afterForces(simulate) in simulation.cpp */
    void afterForces(ParticleContainer *moleculeContainer, DomainDecompBase *domainDecomposition, CellProcessor *cellProcessor,
                         unsigned long simstep) override;

private:

    /** List of ChemicalPotential objects, each of which describes a
	 * particular control volume for the grand canonical ensemble with
	 * respect to one of the simulated components.
	 *
	 * It may at first be unclear why one could want to specify
	 * several grand canonical ensembles, which are then stored in a
	 * list. However, note that for every component a distinct
	 * chemical potential can be specified, and this is of course
	 * essential in certain cases. Also, different chemical potentials
	 * can be specified for different control volumes to induce a
	 * gradient of the chemical potential.
	 */
    std::list<ChemicalPotential> _lmu;

    unsigned long _N;
    double _V;
    double _T;

    double _mu;
    double _p;
    double _E;

    double _E_trans;
    double _E_rot;

    // Taken from simulation.cpp defaults. usually too large to have ever been used
    // Functionality of GrandCanonical not proven, probably lost during move to new input format
    unsigned long _initGrandCanonical = 10000000;

    Domain* _simulationDomain;
};


#endif //MARDYN_GRANDCANONICAL_H
