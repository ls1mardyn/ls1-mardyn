#include "CanonicalEnsemble.h"

#include <map>

#include "BoxDomain.h"
#include "DomainBase.h"
#include "molecules/Component.h"
#include "molecules/Molecule.h"
#ifdef ENABLE_MPI
#include "parallel/CollectiveCommunication.h"
#endif
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"


using namespace std;
using Log::global_log;

void CanonicalEnsemble::updateGlobalVariable( GlobalVariable variable ) {

	const int numComponents = this->numComponents();


	Molecule *tM;


	/* calculate local variables */

	/* "fixed" variables of this ensemble */
	if ( (variable & NUM_PARTICLES) | (variable & TEMPERATURE) ) {
		global_log->debug() << "Updating particle counts" << endl;
		/* initializes the number of molecules present in each component! */
		unsigned long *numMolecules = new unsigned long[numComponents];
		for( int cid = 0; cid < numComponents; cid++) 
			numMolecules[cid] = 0;
		for( tM = _particles->begin(); tM != _particles->end(); tM = _particles->next() ) {
			Molecule& molecule = *tM;
			const int cid = molecule.componentid();
			numMolecules[cid]++;
		}
#ifdef ENABLE_MPI
		_simulation.domainDecomposition().collCommInit(numComponents);
		for( int cid = 0; cid < numComponents; cid++)
			_simulation.domainDecomposition().collCommAppendUnsLong(numMolecules[cid]);
		_simulation.domainDecomposition().collCommAllreduceSum();
#endif
		_N = 0;
		for( int cid = 0; cid < numComponents; cid++) {
#ifdef ENABLE_MPI
			numMolecules[cid] =  _simulation.domainDecomposition().collCommGetUnsLong();
#endif
			global_log->debug() << "Number of molecules in component " << cid << ": " << numMolecules[cid] << endl;
			_N += numMolecules[cid];
			_components[cid].setNumMolecules(numMolecules[cid]);
		}
#ifdef ENABLE_MPI
		_simulation.domainDecomposition().collCommFinalize();
#endif
		delete [] numMolecules;
	}

	if ( variable & VOLUME ) {
		global_log->debug() << "Updating volume" << endl;
	  /* TODO: calculate actual volume or return specified volume as 
	   * the canonical ensemble should have a fixed volume? */
	}

	/* variable variables of this ensemble */
	if ( variable & CHEMICAL_POTENTIAL ) {
		global_log->debug() << "Updating chemical potential" << endl;
	}

	if ( variable & PRESSURE ) {
		global_log->info() << "Updating pressure" << endl;
	}

	if ( (variable & ENERGY) | (variable & TEMPERATURE) ) {
		global_log->debug() << "Updating energy" << endl;
	  double *E_trans = new double [numComponents];
	  double *E_rot = new double[numComponents];
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
#ifdef ENABLE_MPI
	  _simulation.domainDecomposition().collCommInit(2*numComponents);
	  for( int cid = 0; cid < numComponents; cid++ ) {
		  _simulation.domainDecomposition().collCommAppendDouble(E_trans[cid]);
		  _simulation.domainDecomposition().collCommAppendDouble(E_rot[cid]);
	  }
	  _simulation.domainDecomposition().collCommAllreduceSum();
#endif
	  _E = _E_trans = _E_rot = 0.0;
	  for( int cid = 0; cid < numComponents; cid++) {
#ifdef ENABLE_MPI
		  E_trans[cid] =  _simulation.domainDecomposition().collCommGetDouble();
		  E_rot[cid]   =  _simulation.domainDecomposition().collCommGetDouble();
#endif
		  global_log->debug() << "Kinetic energy in component " << cid << ": " << 
			  "E_trans = " << E_trans[cid] << ", E_rot = " << E_rot[cid] << endl;
		  _components[cid].setE_trans(E_trans[cid]);
		  _components[cid].setE_rot(E_rot[cid]);
		  _E_trans += E_trans[cid];
		  _E_rot   += E_rot[cid];
	  }
#ifdef ENABLE_MPI
	  _simulation.domainDecomposition().collCommFinalize();
#endif

	  global_log->debug() << "Total Kinetic energy: 2*E_trans = " << _E_trans 
		  << ", 2*E_rot = " << _E_rot << endl;
	  _E = _E_trans + _E_rot;
		delete [] E_trans;
		delete [] E_rot; 
	}

	if ( variable & TEMPERATURE ) {
		global_log->debug() << "Updating temperature" << endl;
		/* TODO: calculate actual temperature or return specified temperature as 
		 * the canonical ensemble should have a fixed temperature? */
		long long totalDegreesOfFreedom = 0;
		for( int cid = 0; cid < numComponents; cid++) {
			unsigned int rdf = _components[cid].getRotationalDegreesOfFreedom();
			long long N      = _components[cid].getNumMolecules();
			long long degreesOfFreedom = (3 + rdf) * N;
			totalDegreesOfFreedom += degreesOfFreedom;

			double E_kin = _components[cid].E(); 
			double T = E_kin / degreesOfFreedom;
			global_log->debug() << "Temprature of component " << cid << ": " <<
				              "T = " << T << endl;
			_components[cid].setT( T );
		}
		_T = _E / totalDegreesOfFreedom;
	}

	/* now calculate all local variables */
	//molecule.calculate_mv2_Iw2( mv2[cid], Iw2[cid]);

	/* calculate global variables from local variables */

	/* save variables to components and ensemble */
}

void CanonicalEnsemble::readXML(XMLfileUnits& xmlconfig) {
	Ensemble::readXML(xmlconfig);

	xmlconfig.getNodeValueReduced("temperature", _T);
	global_log->info() << "Temperature: " << _T << endl;
	string domaintype;
	xmlconfig.getNodeValue("domain@type", domaintype);
	global_log->info() << "Domain type: " << domaintype << endl;
	if( "box" == domaintype) {
		_domain = new BoxDomain();
	}
	else {
		global_log->error() << "Volume type not supported." << endl;
		exit(1);
	}
	xmlconfig.changecurrentnode("domain");
	_domain->readXML(xmlconfig);
	xmlconfig.changecurrentnode("..");
	_V = _domain->V();
	global_log->info() << "Volume: " << _V << endl;
	global_log->warning() << "Box dimensions not set yet in domain class" << endl;
}
