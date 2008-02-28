#ifndef PARTICLEPAIRS2POTFORCEADAPTER_H_
#define PARTICLEPAIRS2POTFORCEADAPTER_H_

namespace datastructures {
  class ParticlePairs2PotForceAdapter; 
}

//! @brief calculate pair forces and collect macroscopic values
//! @author Martin Buchholz
//!
//! used to calculate the force between all pairs and sum up macroscopic values (e.g. Upot)
//! The idea is, that after the call of init(), processPair(...) is called for all
//! particle pairs in the datastructure. processPair(...) calculates the interaction
//! of the two particles and collects macroscopic values in local member variables.
//! At the end (all pairs have been processed), finish() is called, which stores
//! the macroscopic values in _domain.
class datastructures::ParticlePairs2PotForceAdapter: public datastructures::ParticlePairsHandler<Molecule>{
  public:
    //! Constructor
    ParticlePairs2PotForceAdapter(Domain& domain): _domain(domain){
    }
    
    //! Destructor
    ~ParticlePairs2PotForceAdapter(){
    }
    
    //! @brief initialize macroscopic values
    //!
    //! each pair contributes to the macroscopic values (potential energy,...)
    //! All those values are initialized with zero, and then for each pair, 
    //! they are increased by the pairs contribution
    void init(){
      _virial = 0;
      _upot6LJ = 0;
      _upotXpoles = 0;
      _myRF = 0;
    }
    
    //! @brief calculate macroscopic values
    //!
    //! After all pairs have been processes, Upot and Virial can be calculated
    //! and stored in _domain
    void finish(){
      _domain.setLocalUpot(_upot6LJ/6.+_upotXpoles+_myRF);
      _domain.setLocalVirial(_virial+3.0*_myRF);
    }
    
    //! @brief calculate force between pairs and collect macroscopic contribution
    //!
    //! For all pairs, the force between the two Molecules has to be calculated
    //! and stored in the molecules. For original pairs(pairType 0), the contributions
    //! to the macroscopic values have to be collected
    void processPair(Molecule& particle1, Molecule& particle2, double distanceVector[3], int pairType){
      ParaStrm& params=_domain.getComp2Params()(particle1.componentid(),particle2.componentid());
      params.reset_read();
      if(pairType == 0){
        PotForce(particle1,particle2,params,distanceVector,_upot6LJ,_upotXpoles,_myRF,_virial);
      }
      else if(pairType == 1){
        PotForce(particle1,particle2,params,distanceVector,_dummy1,_dummy2,_dummy3,_dummy4);
      }
    }
    
    
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
    
    //! @brief dummy variable used for pairs which don't contribute to the macroscopic values
    double _dummy1;
    //! @brief dummy variable used for pairs which don't contribute to the macroscopic values
    double _dummy2;
    //! @brief dummy variable used for pairs which don't contribute to the macroscopic values
    double _dummy3;
    //! @brief dummy variable used for pairs which don't contribute to the macroscopic values
    double _dummy4;
};

#endif /*PARTICLEPAIRS2POTFORCEADAPTER_H_*/
