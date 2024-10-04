#include "CanonicalEnsemble.h"

#include <map>

#include "BoxDomain.h"
#include "DomainBase.h"
#include "Domain.h"
#include "molecules/Component.h"
#include "molecules/Molecule.h"

#ifdef ENABLE_MPI
#include "parallel/CollectiveCommunication.h"
#endif

#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "utils/mardyn_assert.h"
#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"


CanonicalEnsemble::CanonicalEnsemble() :
		_N(0), _V(0), _T(0), _mu(0), _p(0), _E(0), _E_trans(0), _E_rot(0) {
	_type = NVT;
}


void CanonicalEnsemble::updateGlobalVariable(ParticleContainer* particleContainer, GlobalVariable variable) {

	const int numComponents = this->numComponents();

	/* calculate local variables */

	/* "fixed" variables of this ensemble */
	if((variable & NUM_PARTICLES) | (variable & TEMPERATURE)) {
		Log::global_log->debug() << "Updating particle counts" << std::endl;
		/* initializes the number of molecules present in each component! */
		std::vector<unsigned long> numMolecules(numComponents, 0ul);

#if defined(_OPENMP)
#pragma omp parallel
#endif
		{
			std::vector<unsigned long> numMolecules_private(numComponents, 0ul);

			for (auto molecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
				 molecule.isValid(); ++molecule) {
				numMolecules_private[molecule->componentid()]++;
			}

#if defined(_OPENMP)
#pragma omp critical
#endif
			{
				for(unsigned int i = 0; i < numMolecules.size(); i++) {
					numMolecules[i] += numMolecules_private[i];
				}

			}
		}
#ifdef ENABLE_MPI
		_simulation.domainDecomposition().collCommInit(numComponents);
		for( int cid = 0; cid < numComponents; cid++)
			_simulation.domainDecomposition().collCommAppendUnsLong(numMolecules[cid]);
		_simulation.domainDecomposition().collCommAllreduceSum();
#endif
		_N = 0;
		for(int cid = 0; cid < numComponents; cid++) {
#ifdef ENABLE_MPI
			numMolecules[cid] =  _simulation.domainDecomposition().collCommGetUnsLong();
#endif
			Log::global_log->debug() << "Number of molecules in component " << cid << ": " << numMolecules[cid] << std::endl;
			_N += numMolecules[cid];
			_components[cid].setNumMolecules(numMolecules[cid]);
		}
#ifdef ENABLE_MPI
		_simulation.domainDecomposition().collCommFinalize();
#endif
	}

	if(variable & VOLUME) {
		Log::global_log->debug() << "Updating volume" << std::endl;
		/* TODO: calculate actual volume or return specified volume as
		 * the canonical ensemble should have a fixed volume? */
	}

	/* variable variables of this ensemble */
	if(variable & CHEMICAL_POTENTIAL) {
		Log::global_log->debug() << "Updating chemical potential" << std::endl;
	}

	if(variable & PRESSURE) {
		Log::global_log->info() << "Updating pressure" << std::endl;
	}

	if((variable & ENERGY) | (variable & TEMPERATURE)) {
		Log::global_log->debug() << "Updating energy" << std::endl;
		std::vector<double> E_trans(numComponents, 0.);
		std::vector<double> E_rot(numComponents, 0.);


#if defined(_OPENMP)
#pragma omp parallel
#endif
		{
			std::vector<unsigned long> E_trans_priv(numComponents, 0.);
			std::vector<unsigned long> E_rot_priv(numComponents, 0.);

			for(auto molecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); molecule.isValid(); ++molecule) {
				const int cid = molecule->componentid();
				double E_trans_loc = 0.0;
				double E_rot_loc = 0.0;
				molecule->calculate_mv2_Iw2(E_trans_loc, E_rot_loc);
				E_trans_priv[cid] += E_trans_loc;  // 2*k_{B} * E_{trans}
				E_rot_priv[cid] += E_rot_loc;  // 2*k_{B} * E_{rot}
			}
#if defined(_OPENMP)
#pragma omp critical
#endif
			{
				for(unsigned int i = 0; i < E_trans_priv.size(); i++) {
					E_trans[i] += E_trans_priv[i];
					E_rot[i] += E_rot[i];
				}
			}
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
		for(int cid = 0; cid < numComponents; cid++) {
#ifdef ENABLE_MPI
			E_trans[cid] =  _simulation.domainDecomposition().collCommGetDouble();
		  E_rot[cid]   =  _simulation.domainDecomposition().collCommGetDouble();
#endif
			Log::global_log->debug() << "Kinetic energy in component " << cid << ": " <<
								"E_trans = " << E_trans[cid] << ", E_rot = " << E_rot[cid] << std::endl;
			_components[cid].setE_trans(E_trans[cid]);
			_components[cid].setE_rot(E_rot[cid]);
			_E_trans += E_trans[cid];
			_E_rot += E_rot[cid];
		}
#ifdef ENABLE_MPI
		_simulation.domainDecomposition().collCommFinalize();
#endif

		Log::global_log->debug() << "Total Kinetic energy: 2*E_trans = " << _E_trans
							<< ", 2*E_rot = " << _E_rot << std::endl;
		_E = _E_trans + _E_rot;
	}

	if(variable & TEMPERATURE) {
		Log::global_log->debug() << "Updating temperature" << std::endl;
		/* TODO: calculate actual temperature or return specified temperature as
		 * the canonical ensemble should have a fixed temperature? */
		long long totalDegreesOfFreedom = 0;
		for(int cid = 0; cid < numComponents; cid++) {
			unsigned int rdf = _components[cid].getRotationalDegreesOfFreedom();
			long long N = _components[cid].getNumMolecules();
			long long degreesOfFreedom = (3 + rdf) * N;
			totalDegreesOfFreedom += degreesOfFreedom;

			double E_kin = _components[cid].E();
			double T = E_kin / degreesOfFreedom;
			Log::global_log->debug() << "Temperature of component " << cid << ": " <<
								"T = " << T << std::endl;
			_components[cid].setT(T);
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
	Log::global_log->info() << "Temperature: " << _T << std::endl;
	std::string domaintype;
	xmlconfig.getNodeValue("domain@type", domaintype);
	Log::global_log->info() << "Domain type: " << domaintype << std::endl;
	if("box" == domaintype) {
		_domain = new BoxDomain();
	} else {
		std::ostringstream error_message;
		error_message << "Volume type not supported." << std::endl;
		MARDYN_EXIT(error_message);
	}
	xmlconfig.changecurrentnode("domain");
	_domain->readXML(xmlconfig);
	xmlconfig.changecurrentnode("..");
	_V = _domain->V();
	Log::global_log->info() << "Volume: " << _V << std::endl;
	Log::global_log->warning() << "Box dimensions not set yet in domain class" << std::endl;
}

void CanonicalEnsemble::beforeThermostat(unsigned long simstep, unsigned long initStatistics) {
	if(simstep >= initStatistics) {
		global_simulation->getDomain()->record_cv();
	}
}
