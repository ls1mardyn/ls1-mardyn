// Wrapper for the molecule to be used with the coupling tool
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
	double* getVelocity() const {
		double* vel = new double[3];
		vel[0] = mol->v(0);
		vel[1] = mol->v(1);
		vel[2] = mol->v(2);
		return vel;
	}
	void setVelocity(double* velocity) {
		mol->setv(0, velocity[0]);
		mol->setv(1, velocity[1]);
		mol->setv(2, velocity[2]);
	}

	/** returns/ sets the position of the molecule */
	double* getPosition() const {
		double* pos = new double[3];
		pos[0] = mol->r(0);
		pos[1] = mol->r(1);
		pos[2] = mol->r(2);
		return pos;
	}
	void setPosition(double* position) {
		mol->setr(0, position[0]);
		mol->setr(1, position[1]);
		mol->setr(2, position[2]);
	}

	/** sets the force acting on this molecule. This function is only called in the USHER
	 *  scheme so far If you want to set the force of a newly created molecule,
	 *  you need to implement this function.
	 */
	void setForce(double* force) {
		mol->setF(force);
	}

	double* getForce() {
		double* f = new double[3];
		f[0] = mol->F(0);
		f[1] = mol->F(1);
		f[2] = mol->F(2);
		return f;
	}


private:
	Molecule* mol;
};
