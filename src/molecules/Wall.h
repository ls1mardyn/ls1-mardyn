//Calculation of the Fluid Wall interaction by a function

#ifndef WALL_H_
#define WALL_H_

#include "molecules/Molecule.h"
#include "molecules/Component.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "Domain.h"

#include <vector>
#include <cmath>
#include <string>
#include <map>

class XMLfileUnits;
class Wall{
	public:
		// constructor and destructor
		Wall();
		~Wall();

		void readXML(XMLfileUnits& xmlconfig);
		void initializeLJ93(const std::vector<Component>* components, double in_rhoWall, double in_sigWall, double in_epsWall, std::vector<double> in_xi, std::vector<double> in_eta, double in_yOffWall, double in_yWallCut);
		void calcTSLJ_9_3( ParticleContainer* partContainer, Domain* domain );
		void initializeLJ104(const std::vector<Component>* components, double in_rhoWall, double in_sigWall, double in_epsWall, std::vector<double> in_xi, std::vector<double> in_eta, double in_yOffWall, double y_cut,  double Delta);
		void calcTSLJ_10_4( ParticleContainer* partContainer, Domain* domain );

	private:
		double _rhoW, _yc, _yOff, Delta;
		double* _eps_wi;
		double* _sig3_wi;
		double* _sig2_wi;
		double* _sig_wi;
		double* _uShift_9_3;
		double* _uPot_9_3;
		double* _uShift_10_4;
		double* _uPot_10_4;
		unsigned _nc;
		double _dWidth;
		double _dWidthHalf;
};

#endif /*WALL_H_*/
