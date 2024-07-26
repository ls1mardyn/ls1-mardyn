//
// Created by Kruegener
//

#pragma once


#include "EnsembleBase.h"
#include "utils/mardyn_assert.h"
#include "ChemicalPotential.h"

class Domain;

class DomainDecompBase;

class CellProcessor;

/** @brief GrandCanonicalEnsemble was moved here from Simulation.cpp and partially adapted to Ensemble base class<br>
 * <br>
 * The GrandCanonicalEnsemble probably used to be a functional muVT ensemble class, that was somewhat dropped with
 * the move to XML inputs. The routines were all over simulation.cpp but never ran.<br>
 * This is a STUB implementation, mainly lacking readXML and updateGlobalVariable.<br>
 * <br>
 * Work done till now:<br>
 *      - GrandCanonical inherits from Ensemble base class like the NVT Ensemble<br>
 *      - Ensemble Base was extended with additional calls at the appropriate points in the simulation loop<br>
 *      - Simulation.h/.cpp is rid of _lmu list and only uses general _ensemble calls<br>
 * <br>
 * To Do:<br>
 *      - Make GrandCanonical functional again<br>
 *      - Adapt oldInput style to XML<br>
 *      - Test against a known simulation result<br>
 */
class GrandCanonicalEnsemble : public Ensemble {

public:
	GrandCanonicalEnsemble();

	~GrandCanonicalEnsemble() override = default;

	// TODO: Implement STUB
	void readXML(XMLfileUnits& xmlconfig) override {
		Log::global_log->info() << "[GrandCanonicalEnsemble] readXML not implemented!" << std::endl;
		mardyn_exit(-1);
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
	void updateGlobalVariable(ParticleContainer* particleContainer, GlobalVariable variable) override {
		Log::global_log->info() << "[GrandCanonicalEnsemble] updateGlobalVariable not implemented!" << std::endl;
		mardyn_exit(-1);
	};

	/*! Runs steps formerly in initConfigXML in simulation.cpp */
	void initConfigXML(ParticleContainer* moleculeContainer) override;

	/*! Runs steps formerly in prepare_start in simulation.cpp */
	void prepare_start() override;

	/*! Runs steps formerly in simulate in simulation.cpp */
	void beforeEventNewTimestep(ParticleContainer* moleculeContainer, DomainDecompBase* domainDecomposition,
								unsigned long simstep) override;

	/*! Runs steps formerly in afterForces(simulate) in simulation.cpp */
	void afterForces(ParticleContainer* moleculeContainer, DomainDecompBase* domainDecomposition,
					 CellProcessor* cellProcessor,
					 unsigned long simstep) override;

	/*! stores a molecule as a sample for a given component */
	void storeSample(Molecule* m, uint32_t componentid) override;

	/*! Returns _lmu pointer for processing by external plugins */
	std::list<ChemicalPotential>* getLmu() override { return &_lmu; }

private:

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

	/** List of ChemicalPotential objects needed for GrandCanonical only, each of which describes a
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
