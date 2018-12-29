//
// Created by Kruegener on 12/29/2018.
//

#ifndef MARDYN_GRANDCANONICAL_H
#define MARDYN_GRANDCANONICAL_H


#include "EnsembleBase.h"

class ChemicalPotential;

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

    // Returns _lmu pointer for processing by external plugins
    std::list<ChemicalPotential>* getLmu() override {return &_lmu;}

    void initConfigXML(ParticleContainer *moleculeContainer, double h) override;

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

    Domain* _simulationDomain;
};


#endif //MARDYN_GRANDCANONICAL_H
