//Calculation of the Fluid Wall interaction by a function

#ifndef MIRROR_H_
#define MIRROR_H_

#include "molecules/Molecule.h"
#include "molecules/Component.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "Domain.h"

#include <vector>
#include <cmath>
#include <string>
#include <map>

class Mirror{
	public:
		// constructor and destructor
		Mirror();
		~Mirror();

		void initialize(const std::vector<Component>* components, double in_yMirr, double in_forceConstant);
		void VelocityChange( ParticleContainer* partContainer, Domain* domain );

	private:
		double  _yMirr;
		double  _forceConstant;
};

#endif /*MIRROR_H_*/
