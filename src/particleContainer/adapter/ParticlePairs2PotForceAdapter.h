
#ifndef PARTICLEPAIRS2POTFORCEADAPTER_H_
#define PARTICLEPAIRS2POTFORCEADAPTER_H_

#include "molecules/potforce.h"
#include "particleContainer/handlerInterfaces/ParticlePairsHandler.h"
#include "io/RDF.h"
#include "Domain.h"
#include "WrapOpenMP.h"

//! @brief calculate pair forces and collect macroscopic values
//! @author Martin Bernreuther <bernreuther@hlrs.de> et al. (2010)
//!
//! used to calculate the force between all pairs and sum up macroscopic values (e.g. Upot)
//! The idea is, that after the call of init(), processPair(...) is called for all
//! particle pairs in the datastructure. processPair(...) calculates the interaction
//! of the two particles and collects macroscopic values in local member variables.
//! At the end (all pairs have been processed), finish() is called, which stores
//! the macroscopic values in _domain.
class ParticlePairs2PotForceAdapter : public ParticlePairsHandler {
public:
	//! Constructor
	ParticlePairs2PotForceAdapter(Domain& domain) :
		_domain(domain), _virial(0.0), _upot6LJ(0.0), _upotXpoles(0.0), _myRF(0.0) {
		//this->_doRecordRDF = false;

		const int numThreads = mardyn_get_max_threads();
		Log::global_log->info() << "ParticlePairs2PotForceAdapter: allocate data for " << numThreads << " threads." << std::endl;
		_threadData.resize(numThreads);
		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			PP2PFAThreadData * myown = new PP2PFAThreadData();
			const int myid = mardyn_get_thread_num();
			_threadData[myid] = myown;
		} // end pragma omp parallel
	}

	//! Destructor
	~ParticlePairs2PotForceAdapter() {
		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			const int myid = mardyn_get_thread_num();
			delete _threadData[myid];
		} // end pragma omp parallel
	}

	//! @brief initialize macroscopic values
	//!
	//! each pair contributes to the macroscopic values (potential energy,...)
	//! All those values are initialized with zero, and then for each pair,
	//! they are increased by the pairs contribution
	void init() {
		_virial = 0;
		_upot6LJ = 0;
		_upotXpoles = 0;
		_myRF = 0;

		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			const int myid = mardyn_get_thread_num();
			_threadData[myid]->initComp2Param(_domain.getComp2Params());
			_threadData[myid]->clear();
		}
	}

	//! @brief calculate macroscopic values
	//!
	//! After all pairs have been processes, Upot and Virial can be calculated
	//! and stored in _domain
	void finish() {
		double glob_virial = 0.0;
		double glob_upot6LJ = 0.0;
		double glob_upotXpoles = 0.0;
		double glob_myRF = 0.0;

		#if defined(_OPENMP)
		#pragma omp parallel reduction(+:glob_virial, glob_upot6LJ, glob_upotXpoles, glob_myRF)
		#endif
		{
			const int myid = mardyn_get_thread_num();
			glob_virial += _threadData[myid]->_virial;
			glob_upot6LJ += _threadData[myid]->_upot6LJ;
			glob_upotXpoles += _threadData[myid]->_upotXpoles;
			glob_myRF += _threadData[myid]->_myRF;
		} // end pragma omp parallel reduction

		_virial = glob_virial;
		_upot6LJ = glob_upot6LJ;
		_upotXpoles = glob_upotXpoles;
		_myRF = glob_myRF;

		_domain.setLocalUpot(_upot6LJ / 6. + _upotXpoles + _myRF);
		_domain.setLocalVirial(_virial + 3.0 * _myRF);
	}

	struct PP2PFAThreadData {
		PP2PFAThreadData() : _virial(0.0), _upot6LJ(0.0), _upotXpoles(0.0), _myRF(0.0), _comp2Param(0) {}

		~PP2PFAThreadData() {
			if(_comp2Param != 0) {
				delete _comp2Param;
				_comp2Param = 0;
			}
		}

		void initComp2Param(Comp2Param& c2p) {
			if (_comp2Param != 0) {
				delete _comp2Param;
			}
			_comp2Param = new Comp2Param(c2p);
		}

		void clear() {
			_virial = 0.0;
			_upot6LJ = 0.0;
			_upotXpoles = 0.0;
			_myRF = 0.0;
		}

		double _virial;
		double _upot6LJ;
		double _upotXpoles;
		double _myRF;
		Comp2Param * _comp2Param;
	};

	std::vector<PP2PFAThreadData *> _threadData;

    /** calculate force between pairs and collect macroscopic contribution
     *
     * For all pairs, the force between the two Molecules has to be calculated
     * and stored in the molecules. For original pairs(pairType 0), the contributions
     * to the macroscopic values have to be collected
     *
     * @param molecule1       molecule 1
     * @param molecule2       molecule 2
     * @param distanceVector  distance vector from molecule 2 to molecule 1
     * @param pairType        molecule pair type (see PairType)
     * @param dd              square of the distance between the two molecules
     * @param calculateLJ     true if we shall calculate the LJ interaction, otherwise false (default true)
     *
     * @return                interaction energy
     */
	double processPair(Molecule& molecule1, Molecule& molecule2, double distanceVector[3], PairType pairType, double dd, bool calculateLJ = true) {
		const int tid = mardyn_get_thread_num();
		PP2PFAThreadData &my_threadData = *_threadData[tid];

		ParaStrm& params = (* my_threadData._comp2Param)(molecule1.componentid(), molecule2.componentid());
		params.reset_read();

		switch (pairType) {

            double dummy1, dummy2, dummy3, dummy4[3], Virial3[3];

            case MOLECULE_MOLECULE :
//                if ( _rdf != NULL ) _rdf->observeRDF(molecule1, molecule2, dd); // moved to RDFCellProcessor
                PotForce(molecule1, molecule2, params, distanceVector, my_threadData._upot6LJ, my_threadData._upotXpoles, my_threadData._myRF, Virial3, calculateLJ );
                my_threadData._virial += 2*(Virial3[0]+Virial3[1]+Virial3[2]);
                return my_threadData._upot6LJ + my_threadData._upotXpoles;
            case MOLECULE_HALOMOLECULE :

                PotForce(molecule1, molecule2, params, distanceVector, dummy1, dummy2, dummy3, dummy4, calculateLJ);
                return 0.0;
            case MOLECULE_MOLECULE_FLUID :
                dummy1 = 0.0; // 6*U_LJ
                dummy2 = 0.0; // U_polarity
                dummy3 = 0.0; // U_dipole_reaction_field

                FluidPot(molecule1, molecule2, params, distanceVector, dummy1, dummy2, dummy3, calculateLJ);
                return dummy1 / 6.0 + dummy2 + dummy3;
            default:
				std::ostringstream error_message;
				error_message << "[ParticlePairs2PotForceAdapter9] pairType is unknown" << std::endl;
                MARDYN_EXIT(error_message);
        }
        return 0.0;
	}


//	void recordRDF() {
//		this->_doRecordRDF = true;
//	}

private:
	//! @brief reference to the domain is needed to store the calculated macroscopic values
	Domain& _domain;

	//! @brief variable used to sum the virial contribution of all pairs
	double _virial;
	//! @brief variable used to sum the Upot6LJ contribution of all pairs
	double _upot6LJ;
	//! @brief variable used to sum the UpotXpoles contribution of all pairs
	double _upotXpoles;
	//! @brief variable used to sum the MyRF contribution of all pairs
	double _myRF;

//	bool _doRecordRDF;
};

#endif /*PARTICLEPAIRS2POTFORCEADAPTER_H_*/
