
namespace moleculardynamics {
  namespace coupling {
    namespace interface {
      template<class Molecule,unsigned int dim>
      class MardynMoleculeWrapper;
    }
  }
}



template<class Molecule, unsigned int dim>
class moleculardynamics::coupling::interface::MardynMoleculeWrapper {
public:
	MardynMoleculeWrapper(Molecule* molecule) : mol(molecule){
	}


	/** returns/ sets the velocity of the molecule */
	tarch::la::Vector<dim,double> getVelocity() const {
		tarch::la::Vector<3, double>(mol->v(0), mol->v(1), mol->v(2));
	}
	void setVelocity(const tarch::la::Vector<dim,double>& velocity) {
		mol->setv(0, velocity(0));
		mol->setv(1, velocity(1));
		mol->setv(2, velocity(2));
	}

	/** returns/ sets the position of the molecule */
	tarch::la::Vector<dim,double> getPosition() const {
		return tarch::la::Vector<3, double>(mol->r(0), mol->r(1), mol->r(2));
	}
	void setPosition(const tarch::la::Vector<dim,double>& position) {
		mol->setr(0, position(0));
		mol->setr(1, position(1));
		mol->setr(2, position(2));
	}

	/** sets the force acting on this molecule. This function is only called in the USHER
	 *  scheme so far If you want to set the force of a newly created molecule,
	 *  you need to implement this function.
	 */
	void setForce(const tarch::la::Vector<dim,double>& force) {
		double f[3];
		for (int i = 0; i < 3; i++) f[i] = force(i);
		mol->setF(f);
	}
	tarch::la::Vector<dim,double> getForce() {
		tarch::la::Vector<3, double> (mol->F(0), mol->F(1), mol->F(2));
	}


private:
	Molecule* mol;
};
