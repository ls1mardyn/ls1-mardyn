/*
 * Usher algorithm for inserting particles on the position
 * where energy is equal to the mean energy of the system
 */

#ifndef _MOLECULARDYNAMICS_COUPLING_PARTICLEINSERTION_H_
#define _MOLECULARDYNAMICS_COUPLING_PARTICLEINSERTION_H_

#include "particleContainer/LinkedCells.h"
#include "molecules/Molecule.h"
#include <vector>
#include "MardynMoleculeWrapper.h"
#include <cstdlib>
#include <cmath>
#include <string>
#include "Simulation.h"
using namespace std;

#define PI 3.1415926535

namespace moleculardynamics {
namespace coupling {
template<class Molecule, class LinkedCells, unsigned int dim>
class ParticleInsertion;
}
}


/** handles particle insertion (via Usher algorithm) and random particle deletion.
 *
 *  @author Philipp Neumann
 */
template<class Molecule, class LinkedCells, unsigned int dim>
class moleculardynamics::coupling::ParticleInsertion {
public:
	ParticleInsertion();
	~ParticleInsertion() {
	}

	/** this state is returned by the insertDeleteMass() function and tells the user, if mass was inserted/deleted
	 *  or if nothing happened at all.
	 */
	enum Action {
		NoAction = 0, Insertion = 1, Deletion = 2
	};




	/*
	 * Determines position where the energy is within the xiMax from the
	 * target energy U_0
	 * Position should be in the allowed region defined by
	 * allowedLow and allowedHigh
	 * Returns -1 if failed, otherwise number of performed steps
	 * At the end, molecule has the correct position
	 * Please view cpph file for parameter explanations
	 */
	int findParticlePosition(
			LinkedCells* linkedCells,
			Molecule* molecule,
			double U_0, double* energy, double* old_energy,
			bool largeStepsizeOnOverlap, bool restartIfIncreases, int seed,
			int intIterMax, int restartMax,
			int rotationsMax, double maxAngle, double maxAllowedAngle, double minAngle, double xiMax,
			vector<double>* vec_energy, vector<double>* vec_angle,
			vector<double*>* vec_lj, vector<double*>* vec_center,
			string name_energy, string name_angle, string name_lj, string name_center,
			double* allowed_low, double* allowed_high) const;

	/*
	 * Molecule rotations towards teh target energy
	 * Please see cpph file for parameter explanations
	 */
	int rotateMolecule(Molecule* molecule,
			moleculardynamics::coupling::interface::MardynMoleculeWrapper<Molecule,
			dim> wrapper, LinkedCells* linkedCells,
			int rotationsMax, double maxAngle, double maxAllowedAngle, double minAngle,
			double U_0, double U_overlap, int* timestep,
			double xiMax, double* energy, double* energyOld, double* force, double* absForce,
			Quaternion* q, bool stopIfOverlapping, bool log,
			vector<double>* vec_energy, vector<double>* vec_angle,
			vector<double*>* vec_lj, vector<double*>* vec_center,
			string name_energy, string name_angle, string name_lj, string name_center) const;
	void writeMoleculeVtk(Molecule* molecule, LinkedCells* linkedCells,
			Simulation* simulation, int* timestep) const;

	void doLogging(vector<double>* vec_energy, vector<double>* vec_angle,
			vector<double*>* vec_lj, vector<double*>* vec_center,
			string name_energy, string name_angle, string name_lj,
			string name_center) const;

	void writeToVectors(vector<double>* vec_energy, vector<double>* vec_angle,
			vector<double*>* vec_lj, vector<double*>* vec_center, double energy, double angle,
			Molecule* molecule) const;

	void write_vector(vector<double> vec, string file_name) const;
	void write_vector_3d(vector<double*> vec, string file_name) const;

};
#include "ParticleInsertion.cpph"
#endif // _MOLECULARDYNAMICS_COUPLING_PARTICLEINSERTION_H_
