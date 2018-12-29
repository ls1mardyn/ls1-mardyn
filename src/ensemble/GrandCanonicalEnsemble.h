//
// Created by Kruegener on 12/29/2018.
//

#ifndef MARDYN_GRANDCANONICAL_H
#define MARDYN_GRANDCANONICAL_H


#include "EnsembleBase.h"
#include "ChemicalPotential.h"

class GrandCanonicalEnsemble : public Ensemble{

    GrandCanonicalEnsemble() {}

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
};


#endif //MARDYN_GRANDCANONICAL_H
