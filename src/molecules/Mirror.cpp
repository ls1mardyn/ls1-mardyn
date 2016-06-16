
#include "Mirror.h"

using namespace std;
using Log::global_log;


Mirror::Mirror(){}

Mirror::~Mirror(){}

void Mirror::initialize(const std::vector<Component>* /*components*/, double in_yMirr, double in_forceConstant)
{
  global_log->info() << "Initializing the mirror function.\n";
  this->_yMirr = in_yMirr;
  this->_forceConstant = in_forceConstant;
       
}

void Mirror::VelocityChange( ParticleContainer* partContainer, Domain* /*domain*/)
{

  double regionLowCorner[3], regionHighCorner[3];
  vector<Molecule*> particlePtrsForRegion;
  
  //global_log->info() << "Mirror is called\n";

  particlePtrsForRegion.clear();
  
  /*!*** Mirror boundary in y-direction on top of the simulation box ****/
  if(partContainer->getBoundingBoxMax(1) > _yMirr){ // if linked cell in the region of the mirror boundary
    for(unsigned d = 0; d < 3; d++){
     regionLowCorner[d] = partContainer->getBoundingBoxMin(d);
     regionHighCorner[d] = partContainer->getBoundingBoxMax(d);
    }
    // regionLowCorner[1] = (partContainer->getBoundingBoxMin(1) > _yMirr) ? (partContainer->getBoundingBoxMin(1)) : _yMirr;
    partContainer->getRegion(regionLowCorner, regionHighCorner, particlePtrsForRegion);

    std::vector<Molecule*>::iterator particlePtrIter;
    
    for(particlePtrIter = particlePtrsForRegion.begin(); particlePtrIter != particlePtrsForRegion.end(); particlePtrIter++){
		double additionalForce[3];
		additionalForce[0]=0;
		additionalForce[2]=0;
		//global_log->info() << "Mirror\n";
		if ((*particlePtrIter)->r(1) > _yMirr){
			//global_log->info()<<"particle groesser als mirror\n";
			//double _forceConstant = 1*10^(20); 
			double distance = (*particlePtrIter)->r(1) - _yMirr;
			additionalForce[1] = -_forceConstant*distance;
			//global_log->info() << additionalForce[0] << "\n";
			//global_log->info() << additionalForce[1] << "\n";
			//global_log->info() << additionalForce[2] << "\n";
			(*particlePtrIter)->Fljcenteradd(0,additionalForce);
		}
 //     if((*particlePtrIter)->v(1) > 0){
	//(*particlePtrIter)->setv(1, -(*particlePtrIter)->v(1)); 
	// }   
   }    
  } // end Mirror boundary

} // end mthod calcTSLJ_9_3(...)

