/***************************************************************************
 *   Copyright (C) 2010 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef PARTICLEPAIRS2POTFORCEADAPTER_H_
#define PARTICLEPAIRS2POTFORCEADAPTER_H_

#include "molecules/potforce.h"
#include "datastructures/handlerInterfaces/ParticlePairsHandler.h"

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
  ParticlePairs2PotForceAdapter(Domain& domain) : _domain(domain) {
    this->_doRecordRDF = false;
  }

  //! Destructor
  ~ParticlePairs2PotForceAdapter() {
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
    _upotTersoff = 0;
  }

  //! @brief calculate macroscopic values
  //!
  //! After all pairs have been processes, Upot and Virial can be calculated
  //! and stored in _domain
  void finish() {
    _domain.setLocalUpot(_upot6LJ / 6. + _upotXpoles + _upotTersoff + _myRF);
    _domain.setLocalVirial(_virial + 3.0 * _myRF);
  }

  //! @brief calculate force between pairs and collect macroscopic contribution
  //!
  //! For all pairs, the force between the two Molecules has to be calculated
  //! and stored in the molecules. For original pairs(pairType 0), the contributions
  //! to the macroscopic values have to be collected
  double processPair(Molecule& particle1, Molecule& particle2, double distanceVector[3], int pairType, double dd, bool calculateLJ) {
    ParaStrm& params = _domain.getComp2Params()(particle1.componentid(), particle2.componentid());
    params.reset_read();
    if (pairType == 0) {
      if (this->_doRecordRDF)
        this->_domain.observeRDF(dd, particle1.componentid(), particle2.componentid());

      PotForce(particle1, particle2, params, distanceVector, _upot6LJ, _upotXpoles, _myRF, _virial, calculateLJ);
      return _upot6LJ + _upotXpoles;
    }
    else if (pairType == 1) {
      PotForce(particle1, particle2, params, distanceVector, _dummy1, _dummy2, _dummy3, _dummy4, calculateLJ);
      return 0.0;
    }
    else if (pairType == 2) {
      _dummy1 = 0.0; // 6*U_LJ
      _dummy2 = 0.0; // U_polarity
      _dummy3 = 0.0; // U_dipole_reaction_field

      FluidPot(particle1, particle2, params, distanceVector, _dummy1, _dummy2, _dummy3, calculateLJ);
      return _dummy1 / 6.0 + _dummy2 + _dummy3;
    }
    else
      exit(666);
  }

  //! Only for so-called original pairs (pair type 0) the contributions
  //! to the macroscopic values have to be collected
  //!
  //! @brief register Tersoff neighbours
  void preprocessTersoffPair(Molecule& particle1, Molecule& particle2, bool pairType) {
    particle1.addTersoffNeighbour(&particle2, pairType);
    particle2.addTersoffNeighbour(&particle1, pairType);
  }

  //! @brief process Tersoff interaction
  //!
  void processTersoffAtom(Molecule& particle1, double params[15], double delta_r) {
    TersoffPotForce(&particle1, params, _upotTersoff, delta_r);
  }

  void recordRDF() {
    this->_doRecordRDF = true;
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
  //! @brief variable used to sum the Tersoff internal energy contribution of all pairs
  double _upotTersoff;
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

  bool _doRecordRDF;
};

#endif /*PARTICLEPAIRS2POTFORCEADAPTER_H_*/
