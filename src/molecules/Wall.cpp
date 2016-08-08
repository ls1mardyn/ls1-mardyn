
#include "Wall.h"


using namespace std;
using Log::global_log;


Wall::Wall(){}

Wall::~Wall(){}

void Wall::initializeLJ93(const std::vector<Component>* components, 
		 double in_rhoWall, double in_sigWall, double in_epsWall, double* in_xi, double* in_eta, 
		 double in_yOffWall, double in_yWallCut)
{
  global_log->info() << "Initializing the wall function LJ-9-3.\n";
  this->_rhoW = in_rhoWall;
  this->_yc = in_yWallCut; 
  this->_yOff = in_yOffWall;
  

  /*!*** So far: only 1CLJ components allowed ****/
  _nc = components->size();
  _eps_wi = new double[_nc];
  _sig3_wi = new double[_nc];
  _uShift_9_3 = new double[_nc];
  _uPot_9_3 = new double[_nc];

  for (unsigned i = 0; i < _nc; i++)
  {
    _eps_wi[i] = in_xi[i]*sqrt(in_epsWall * components->at(i).ljcenter(0).eps());
    double sig_wi;
    sig_wi = 0.5*in_eta[i]*(in_sigWall + components->at(i).ljcenter(0).sigma());
    _sig3_wi[i] = sig_wi *sig_wi *sig_wi;
    double sig9_sf = _sig3_wi[i] *_sig3_wi[i] *_sig3_wi[i];
    double y = _yc - _yOff;
    double y3 = y*y*y;
    double y9 = y3*y3*y3;
    _uShift_9_3[i] = 4.0/3.0*M_PI*_rhoW*_eps_wi[i]*_sig3_wi[i]*(sig9_sf/15.0/y9 - _sig3_wi[i]/2.0/y3);
    _uPot_9_3[i] = 0.0;

  }
  
    
}

void Wall::initializeLJ104(const std::vector<Component>* components, 
		 double in_rhoWall, double in_sigWall, double in_epsWall, double* in_xi, double* in_eta, 
		 double in_yOffWall, double in_yWallCut, double Delta)
{
  global_log->info() << "Initializing the wall function LJ-10-4.\n";
  this->_rhoW = in_rhoWall;
  this->_yc = in_yWallCut;
  this->_yOff = in_yOffWall;
  

  /*!*** So far: only 1CLJ components allowed ****/
  _nc = components->size();
  _eps_wi = new double[_nc];
  _sig3_wi = new double[_nc];
  _sig2_wi = new double[_nc];
  _sig_wi = new double[_nc];
  _uShift_10_4 = new double[_nc];
  _uPot_10_4 = new double[_nc];


  for (unsigned i = 0; i < _nc; i++)
  {
    _eps_wi[i] = in_xi[i]*sqrt(in_epsWall * components->at(i).ljcenter(0).eps());
    double sig_wi;
    sig_wi = 0.5*in_eta[i]*(in_sigWall + components->at(i).ljcenter(0).sigma());
	_sig_wi[i] = sig_wi;
    _sig3_wi[i] = sig_wi *sig_wi *sig_wi;
	_sig2_wi[i] = sig_wi *sig_wi;
    double sig10_sf = _sig3_wi[i] *_sig3_wi[i] *_sig3_wi[i]*_sig_wi[i];
	double sig4_sf = _sig2_wi[i]*_sig2_wi[i];
	double y = _yc-_yOff;
	double y2 = y*y;
	double y4 = y2*y2;
	double y10 = y4*y4*y2;
	_uShift_10_4[i] = 2*M_PI*_eps_wi[i]*_rhoW*_sig2_wi[i]*Delta*(2/5*sig10_sf/y10-sig4_sf/y4-sig4_sf/(3*Delta*((y + 0.61*Delta)*(y + 0.61*Delta)*(y + 0.61*Delta))));
	_uPot_10_4[i] = 0.0;

  }
  
    
}

void Wall::calcTSLJ_9_3( ParticleContainer* partContainer, Domain* domain)
{

  double regionLowCorner[3], regionHighCorner[3];
  vector<Molecule*> particlePtrsForRegion;
  
  /*! LJ-9-3 potential applied in y-direction */
  if(partContainer->getBoundingBoxMin(1) < _yc){ // if linked cell within the potential range (inside the potential's cutoff)
    for(unsigned d = 0; d < 3; d++){
     regionLowCorner[d] = partContainer->getBoundingBoxMin(d);
     regionHighCorner[d] = partContainer->getBoundingBoxMax(d);
    }
    //regionHighCorner[1] = (partContainer->getBoundingBoxMax(1) < _yc) ? (partContainer->getBoundingBoxMax(1)) : ;
    partContainer->getRegion(regionLowCorner, regionHighCorner, particlePtrsForRegion);
    
    std::vector<Molecule*>::iterator particlePtrIter;

    for(particlePtrIter = particlePtrsForRegion.begin(); particlePtrIter != particlePtrsForRegion.end(); particlePtrIter++){
      //! so far for 1CLJ only, several 1CLJ-components possible
      double y, y3, y9;
      unsigned cid = (*particlePtrIter) -> componentid();
      y = (*particlePtrIter)->r(1) - _yOff;
      if(y < _yc){
	y3 = y*y*y;
	y9 = y3*y3*y3;
	double f[3];
	for(unsigned d = 0; d < 3; d++)
	  f[d] = 0.0;
	
	double sig9_wi;
	sig9_wi = _sig3_wi[cid]*_sig3_wi[cid]*_sig3_wi[cid]; 
	f[cid] = 4.0*M_PI* _rhoW * _eps_wi[cid] * _sig3_wi[cid] * ( sig9_wi/5.0/y9 - _sig3_wi[cid]/2.0/y3 ) / y;
	_uPot_9_3[cid] += 4.0*M_PI* _rhoW * _eps_wi[cid] * _sig3_wi[cid] * ( sig9_wi/45.0/y9 - _sig3_wi[cid]/6.0/y3 ) - _uShift_9_3[cid];
	f[0] = 0;
	f[2] = 0;
	(*particlePtrIter)->Fljcenteradd(0, f);
	
      } // end if()
    }
  }
  
  double u_pot;
  u_pot = _uPot_9_3[0] + domain -> getLocalUpotCompSpecific();
  domain->setLocalUpotCompSpecific(u_pot);
  for(unsigned cid = 0; cid < _nc; cid++){
    _uPot_9_3[cid] = 0.0;
  }
  
  particlePtrsForRegion.clear();
  

} // end mthod calcTSLJ_9_3(...)

void Wall::calcTSLJ_10_4( ParticleContainer* partContainer, Domain* domain)
{

  double regionLowCorner[3], regionHighCorner[3];
  vector<Molecule*> particlePtrsForRegion;
  
  /*! LJ-10-4 potential applied in y-direction */
  if(partContainer->getBoundingBoxMin(1)){ // if linked cell within the potential range (inside the potential's cutoff)
    for(unsigned d = 0; d < 3; d++){
     regionLowCorner[d] = partContainer->getBoundingBoxMin(d);
     regionHighCorner[d] = partContainer->getBoundingBoxMax(d);
    }
    //regionHighCorner[1] = (partContainer->getBoundingBoxMax(1) < _yc) ? (partContainer->getBoundingBoxMax(1)) : ;
    partContainer->getRegion(regionLowCorner, regionHighCorner, particlePtrsForRegion);
    
    std::vector<Molecule*>::iterator particlePtrIter;

    for(particlePtrIter = particlePtrsForRegion.begin(); particlePtrIter != particlePtrsForRegion.end(); particlePtrIter++){
      //! so far for 1CLJ only, several 1CLJ-components possible
      double y, y2, y4, y5, y10, y11;
      unsigned cid = (*particlePtrIter) -> componentid();
      y = (*particlePtrIter)->r(1) - _yOff;
	  if(y < _yc){
	y2 = y*y;
	y4 = y2*y2;
	y5 = y4*y;
	y10 = y5*y5;
	y11 = y10*y;
	double f[3];
	for(unsigned d = 0; d < 3; d++)
	  f[d] = 0.0;
	
	double sig4_wi, sig2_wi, sig5_wi, sig10_wi, bracket, bracket3, term3, preFactor, force[3];
	//double term1, term2;
	//double Potential;
	//double pressure;
	sig2_wi = _sig2_wi[cid];
	sig4_wi = _sig2_wi[cid] * _sig2_wi[cid];
	sig5_wi = sig4_wi * _sig_wi[cid];
	sig10_wi = sig5_wi * sig5_wi;
	bracket = y + 0.61*Delta;
	bracket3 = bracket*bracket*bracket;
	double term1 = sig10_wi/y10;
	double term2 = sig4_wi/y4;
	term3 = sig4_wi/(3*Delta*bracket3);
	preFactor = 2*M_PI*_eps_wi[cid]*_rhoW*sig2_wi*Delta;
	//Potential[cid] = preFactor*(2/5*term1-term2-term3);
	_uPot_10_4[cid] += preFactor*(2/5*term1-term2-term3) - _uShift_10_4[cid];
	force[cid] = preFactor*(4*sig10_wi/y11-4*sig4_wi/y5-term3*3/(bracket));
	//pressure[cid] = 0.5*force[cid]*y;
	f[0] = 0;
	f[2] = 0;
	f[1] = force[1];
	(*particlePtrIter)->Fljcenteradd(0, f);

	  } //end if()
    }
  }
  
  double u_pot;
  u_pot = _uPot_10_4[0] + domain -> getLocalUpotCompSpecific();
  domain->setLocalUpotCompSpecific(u_pot);
  for(unsigned cid = 0; cid < _nc; cid++){
    _uPot_10_4[cid] = 0.0;
  }
  
  particlePtrsForRegion.clear();
  

} // end mthod calcTSLJ_10_4(...)