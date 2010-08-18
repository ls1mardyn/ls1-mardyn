#include "ensemble/CanonicalEnsemble.h"
#include "utils/Logger.h"
#include "particleContainer/ParticleContainer.h"
#include "molecules/Molecule.h"
#include "molecules/Component.h"
#ifdef PARALLEL
#include "parallel/CollectiveCommunication.h"
#endif
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"

#include <map>

using namespace std;
using Log::global_log;

void CanonicalEnsemble::updateGlobalVariable( GlobalVariable variable ) {

	const int numComponents = this->numComponents();


	Molecule *tM;


	/* calculate local variables */

	/* "fixed" variables of this ensemble */
	if ( variable & NUM_PARTICLES ) {
		global_log->info() << "Updating particle counts" << endl;
		/* initializes the number of molecules present in each component! */
		unsigned long numMolecules[numComponents];
		for( int cid = 0; cid < numComponents; cid++) 
			numMolecules[cid] = 0;
		for( tM = _particles->begin(); tM != _particles->end(); tM = _particles->next() ) {
			Molecule& molecule = *tM;
			const int cid = molecule.componentid();
			numMolecules[cid]++;
		}
#ifdef PARALLEL
		_simulation.domainDecomposition().collCommInit(numComponents);
		for( int cid = 0; cid < numComponents; cid++)
			_simulation.domainDecomposition().collCommAppendUnsLong(numMolecules[cid]);
		_simulation.domainDecomposition().collCommAllreduceSum();
#endif
		_N = 0;
		for( int cid = 0; cid < numComponents; cid++) {
#ifdef PARALLEL
			numMolecules[cid] =  _simulation.domainDecomposition().collCommGetUnsLong();
#endif
			global_log->debug() << "Number of molecules in component " << cid << ": " << numMolecules[cid] << endl;
			_N += numMolecules[cid];
			(*_components)[cid].setNumMolecules(numMolecules[cid]);
		}
	}

	if ( variable & VOLUME ) {
		global_log->info() << "Updating volume" << endl;
	  /* TODO: calculate actual volume or return specified volume as 
	   * the canonical ensemble should have a fixed volume? */
	}

	/* variable variables of this ensemble */
	if ( variable & CHEMICAL_POTENTIAL ) {
		global_log->info() << "Updating chemical potential" << endl;
	}

	if ( variable & PRESSURE ) {
		global_log->info() << "Updating pressure" << endl;
	}

	if ( variable & ENERGY | variable & TEMPERATURE) {
		global_log->info() << "Updating energy" << endl;
	  double E_trans[numComponents];
	  double E_rot[numComponents];
	  for( int cid = 0; cid < numComponents; cid++)
		  E_trans[cid] = E_rot[cid] = 0.0;
	  for( tM = _particles->begin(); tM != _particles->end(); tM = _particles->next() ) {
		  Molecule& molecule = *tM;
		  const int cid = molecule.componentid();
		  double E_trans_loc = 0.0;
		  double E_rot_loc = 0.0;
		  molecule.calculate_mv2_Iw2( E_trans_loc, E_rot_loc );
		  E_trans[cid] += E_trans_loc;  // 2*k_{B} * E_{trans}
		  E_rot[cid]   += E_rot_loc;  // 2*k_{B} * E_{rot}
	  }
#ifdef PARALLEL
	  _simulation.domainDecomposition().collCommInit(2*numComponents);
	  for( int cid = 0; cid < numComponents; cid++ ) {
		  _simulation.domainDecomposition().collCommAppendDouble(E_trans[cid]);
		  _simulation.domainDecomposition().collCommAppendDouble(E_rot[cid]);
	  }
	  _simulation.domainDecomposition().collCommAllreduceSum();
#endif
	  _E = _E_trans = _E_rot = 0.0;
	  for( int cid = 0; cid < numComponents; cid++) {
#ifdef PARALLEL
		  E_trans[cid] =  _simulation.domainDecomposition().collCommGetDouble();
#endif
		  global_log->debug() << "Kinetic energy in component " << cid << ": " << 
			  "E_trans = " << E_trans[cid] << ", E_rot = " << E_rot[cid] << endl;
		  (*_components)[cid].setE_trans(E_trans[cid]);
		  (*_components)[cid].setE_rot(E_rot[cid]);
		  _E_trans += E_trans[cid];
		  _E_rot   += E_rot[cid];
	  }

	  global_log->debug() << "Total Kinetic energy: 2*E_{trans} = " << _E_trans 
		  << ", 2*E_{rot} = " << _E_rot << endl;
	  _E = _E_trans + _E_rot;


	}

	if ( variable & TEMPERATURE ) {
		global_log->info() << "Updating temperature" << endl;
		/* TODO: calculate actual temperature or return specified temperature as 
		 * the canonical ensemble should have a fixed temperature? */
		/* TODO: rotDOF missing */
		_T = _E / (double)(3 * _N /*+ rotDOF*/);
	}

	/* now calculate all local variables */
	//molecule.calculate_mv2_Iw2( mv2[cid], Iw2[cid]);

	/* calculate global variables from local variables */

	/* save variables to components and ensemble */
}

