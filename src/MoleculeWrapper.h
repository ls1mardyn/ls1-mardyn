#ifndef _MOLECULARDYNAMICS_COUPLING_INTERFACE_MOLECULEWRAPPER_H_
#define _MOLECULARDYNAMICS_COUPLING_INTERFACE_MOLECULEWRAPPER_H_

#include "tarch/la/Vector.h"
namespace moleculardynamics {
  namespace coupling {
    namespace interface {
      template<class Molecule,unsigned int dim>
      class MoleculeWrapper;
      template<class Molecule,unsigned int dim>
      class ConstMoleculeWrapper;
    }
  }
}


/** wrapper around a molecule in order to access its data.
 *
 *  @author Philipp Neumann
 */
template<class Molecule,unsigned int dim>
class moleculardynamics::coupling::interface::MoleculeWrapper {
  public:
    MoleculeWrapper(Molecule molecule);
    ~MoleculeWrapper();

    /** returns/ sets the velocity of the molecule */
    tarch::la::Vector<dim,double> getVelocity() const;
    void setVelocity(const tarch::la::Vector<dim,double>& velocity);

    /** returns/ sets the position of the molecule */
    tarch::la::Vector<dim,double> getPosition() const;
    void setPosition(const tarch::la::Vector<dim,double>& position);

    /** sets the force acting on this molecule. This function is only called in the USHER
     *  scheme so far If you want to set the force of a newly created molecule,
     *  you need to implement this function.
     */
    void setForce(const tarch::la::Vector<dim,double>& force);
    tarch::la::Vector<dim,double> getForce() const;

    /** returns/ sets the potential energy of the molecule */
    double getPotentialEnergy() const;
    void setPotentialEnergy(const double& potentialEnergy);
};


#endif // _MOLECULARDYNAMICS_COUPLING_INTERFACE_MOLECULEWRAPPER_H_
