#ifndef SIMPLEMOLECULE_H_
#define SIMPLEMOLECULE_H_
//#include "testsuite/molecule_cu_suite.h"

//namespace molecules {
//  class SimpleMolecule;
//}

//! @brief representation of a real Molecule
//!
//! The class SimleMolecule is used to represent molecules and operations
//! on them, e.g. methods to calculate the forces between molecules
//! or update the position of a molecule,...
//! This class is mainly for test purpuses, as it is only a very basic molecule model
class SimpleMolecule{
  public:
  	//! construct a new molecule from the given coordinates
  	SimpleMolecule(int id, int type,
  	         double xPos, double yPos, double zPos,
  	         double xVel, double yVel, double zVel);

  	//! returns the id of the molecule
  	int getId();
    //! returns the type of the molecule
  	int getType();
    //! returns the element of the position vector for the given dimension 
    double getPosition(int dimension);
    //! returns the element of the velocity vector for the given dimension
    double getVelocity(int dimension);
    //! returns the element of the force vector for the given dimension
	  double getForce(int dimension);
	
    //! sets the id to the given value
    void setId(int id);
    //! sets the type to the given value
    void setType(int type);
    //! sets int the position vector the element at the given dimension to the given value
    void setPosition(double position, int dimension);
    //! sets int the velocity vector the element at the given dimension to the given value
    void setVelocity(double velocity, int dimension);
    //! sets int the force vector the element at the given dimension to the given value
    void setForce(double force, int dimension);

	  //! calculate and return the squared distance between the two molecules
	  static double calcDistanceSquare(SimpleMolecule& molecule1, SimpleMolecule& molecule2);

	  //! calculate the force between molecule1 and molecule2 and store it in both molecules
	  static void calcForce(SimpleMolecule& molecule1, SimpleMolecule& molecule2);

    
  private:
  	//! Id of the molecule
  	int id;
  	//! type of the molecule
  	int type;	
   	//! Array for the Position (x, y and z-value)
	  double position[3];
  	//! Array for the Velocity (x, y and z-value)
  	double velocity[3];
  	//! Array for the Force (x, y and z-value)
  	double force[3];

};
	
#endif /*SIMPLEMOLECULE_H_*/
