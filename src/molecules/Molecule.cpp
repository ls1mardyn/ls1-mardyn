#include "Molecule.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iterator>

#include "utils/Logger.h"

#include "Domain.h"

using namespace std;
using Log::global_log;


Molecule::Molecule(unsigned long id, Component *component,
	                 double rx, double ry, double rz,
	                 double vx, double vy, double vz,
	                 double q0, double q1, double q2, double q3,
	                 double Dx, double Dy, double Dz)
		: _q(q0, q1, q2, q3) {
	_id = id;
	_component = component;
	_r[0] = rx;
	_r[1] = ry;
	_r[2] = rz;
	_rOld[0] = rx;
	_rOld[1] = ry;
	_rOld[2] = rz;
	_v[0] = vx;
	_v[1] = vy;
	_v[2] = vz;
	_L[0] = Dx;
	_L[1] = Dy;
	_L[2] = Dz; 
	_sites_d = _sites_F = _osites_e = NULL;
	//_springSites_F = NULL;
	_numTersoffNeighbours = 0;
	fixedx = rx;
	fixedy = ry;
	
	for(int d = 0; d < 3; d++){
	  _directedVelocity[d] = 0.0;
	  for(int e = 0; e < 3; e++){
	    setVirialForce(d, e, 0.0);
	    setVirialKin(d, e, 0.0);
	    setVirialForceConfinement(d, e, 0.0);
	    setVirialKinConfinement(d, e, 0.0);
	  }
	  setPressureVirial(d, 0.0);
	  setPressureKin(d, 0.0);
	  setPressureVirialConfinement(d, 0.0);
	  setPressureKinConfinement(d, 0.0);
	  setPressureVirial_barostat(d, 0.0);
	  setPressureKin_barostat(d, 0.0);
	  setDiffusiveHeatflux(d, 0.0);
	  setConvectivePotHeatflux(d, 0.0);
	  setDirectedVelocity(d, 0.0);
	  setDirectedVelocitySlab(d, 0.0);
	  setDirectedVelocityStress(d, 0.0);
	  setDirectedVelocityConfinement(d, 0.0);
	}
	
	// VTK Molecule Data
	_T = 0.0;
// 	_rho = 0.0;
	_vAverage[0] = 0.0;
	_vAverage[1] = 0.0;
	_vAverage[2] = 0.0;
	_count = 0;
	
	// Hardy stresses
	_HardyStress = false;
	_HardyConfinement = false;
	_virialForceHardyStress.clear();
	_virialForceHardyConfinement.clear();
	_weightingFuncStress = string("Linear");
	_weightingFuncConfinement = string("Linear"); 
	if(_component != NULL) {
		setupCache();
	}
}

Molecule::Molecule(const Molecule& m) {
	_id = m._id;
	_component = m._component;
	_r[0] = m._r[0];
	_r[1] = m._r[1];
	_r[2] = m._r[2];
	_rOld[0] = m._rOld[0];
	_rOld[1] = m._rOld[1];
	_rOld[2] = m._rOld[2];
	_v[0] = m._v[0];
	_v[1] = m._v[1];
	_v[2] = m._v[2];
	_q = m._q;
	_L[0] = m._L[0];
	_L[1] = m._L[1];
	_L[2] = m._L[2];
	_F[0] = m._F[0];
	_F[1] = m._F[1];
	_F[2] = m._F[2];
	_M[0] = m._M[0];
	_M[1] = m._M[1];
	_M[2] = m._M[2];
	_sites_d = _sites_F = _osites_e = NULL;
	//_springSites_F = NULL;
	fixedx = m.fixedx;
	fixedy = m.fixedy;

	for(int d = 0; d < 3; d++){
	  _directedVelocity[d] = 0.0;
	  for(int e = 0; e < 3; e++){
	    setVirialForce(d, e, 0.0);
	    setVirialKin(d, e, 0.0);
	    setVirialForceConfinement(d, e, 0.0);
	    setVirialKinConfinement(d, e, 0.0);
	  }
	  setPressureVirial(d, 0.0);
	  setPressureKin(d, 0.0);
	  setPressureVirialConfinement(d, 0.0);
	  setPressureKinConfinement(d, 0.0);
	  setPressureVirial_barostat(d, 0.0);
	  setPressureKin_barostat(d, 0.0);
	  setDiffusiveHeatflux(d, 0.0);
	  setConvectivePotHeatflux(d, 0.0);
	  setDirectedVelocity(d, 0.0);
	  setDirectedVelocitySlab(d, 0.0);
	  setDirectedVelocityStress(d, 0.0);
	  setDirectedVelocityConfinement(d, 0.0);
	}
	
	// VTK Molecule Data
	_T = m._T;
// 	_rho = m._rho;
	_vAverage[0] = m._vAverage[0];
	_vAverage[1] = m._vAverage[1];
	_vAverage[2] = m._vAverage[2];
	_count = m._count;
	
	// Hardy stresses
	_HardyStress = m._HardyStress;
	_HardyConfinement = m._HardyConfinement;
	_virialForceHardyStress.clear();
	_virialForceHardyConfinement.clear();
	_weightingFuncStress = m._weightingFuncStress;
	_weightingFuncConfinement = m._weightingFuncConfinement;
	if(_component != NULL) {
		setupCache();
	}
}

void Molecule::upd_preF(double dt, double vcorr, double Dcorr, Domain *dom) {
	assert(_m > 0);
	double dt_halve = .5 * dt;
	double dtInv2m = dt_halve / _m;
	int thermostat = dom->getThermostat(this->componentid());
	int dim;
	if(dom->isScaling1Dim(thermostat)  && dom->getAlphaTransCorrection(thermostat) == false){
	    map<int, int> Dim = dom->getDim();
	    dim = Dim[thermostat];
		for (unsigned short d = 0; d < 3; ++d) {
		   if(d == dim)
		      _v[d] = vcorr * _v[d] + dtInv2m * _F[d];
		   else
		      _v[d] = _v[d] + dtInv2m * _F[d];
		   
		   //rOld necessary for the calculation of the diffusion coefficient
		    _rOld[d] = _r[d];  
		   
		   _r[d] += dt * _v[d];
		}
	}
	else{
	    for (unsigned short d = 0; d < 3; ++d) {
		_v[d] = vcorr * _v[d] + dtInv2m * _F[d];
		
		//rOld necessary for the calculation of the diffusion coefficient
		 _rOld[d] = _r[d];
		
		_r[d] += dt * _v[d];
	    }
	}

	double w[3];
	_q.rotate(_L, w);
	for (unsigned short d = 0; d < 3; ++d)
		w[d] *= _invI[d];
	Quaternion qhalfstep;
	_q.differentiate(w, qhalfstep);
	qhalfstep.scale(dt_halve);
	qhalfstep.add(_q);
	double qcorr = 1. / sqrt(qhalfstep.magnitude2());
	qhalfstep.scale(qcorr);
	for (unsigned short d = 0; d < 3; ++d)
		_L[d] = Dcorr * _L[d] + dt_halve * _M[d];
	qhalfstep.rotate(_L, w);
	for (unsigned short d = 0; d < 3; ++d)
		w[d] *= _invI[d];
	Quaternion qincr;
	qhalfstep.differentiate(w, qincr);
	qincr.scale(dt);
	_q.add(qincr);
	qcorr = 1. / sqrt(_q.magnitude2());
	_q.scale(qcorr);

}

void Molecule::upd_cache() {
	unsigned int i;
	unsigned int ns;

	_q.normalize();

	ns = numLJcenters();
	for (i = 0; i < ns; ++i)
		_q.rotateinv(_component->ljcenter(i).r(), &(_ljcenters_d[i*3]));
	ns = numCharges();
	for (i = 0; i < ns; ++i)
		_q.rotateinv(_component->charge(i).r(), &(_charges_d[i*3]));
	ns = numDipoles();
	for (i = 0; i < ns; ++i) {
		const Dipole& di = _component->dipole(i);
		_q.rotateinv(di.r(), &(_dipoles_d[i*3]));
		_q.rotateinv(di.e(), &(_dipoles_e[i*3]));
	}
	ns = numQuadrupoles();
	for (i = 0; i < ns; ++i) {
		const Quadrupole& qi = _component->quadrupole(i);
		_q.rotateinv(qi.r(), &(_quadrupoles_d[i*3]));
		_q.rotateinv(qi.e(), &(_quadrupoles_e[i*3]));
	}
	ns = numTersoff();
	for (i = 0; i < ns; i++)
		_q.rotateinv(_component->tersoff(i).r(), &(_tersoff_d[i*3]));
}

void Molecule::upd_postF(double dt_halve, double& summv2, double& summv2_1Dim, double& sumIw2, Domain *dom) {
	double dtInv2m = dt_halve / _m;
	double v2 = 0.0;
	double v_aux = 0.0;
	for (unsigned short d = 0; d < 3; d++) {
		_v[d] += dtInv2m * _F[d];
		// in the case of directed velocities
		v_aux = _v[d] - _directedVelocity[d];
		v2 += v_aux * v_aux;
		_L[d] += dt_halve * _M[d];
	}
    assert(!isnan(v2)); // catches NaN
    summv2 += _m * v2;
 
	double w[3];
	_q.rotate(_L, w); // L = D = Iw
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
    assert(!isnan(Iw2)); // catches NaN
	sumIw2 += Iw2;
	
  // one-dimensional thermostat
	if(dom->isScaling1Dim(dom->getThermostat(this->componentid()))){
	    map<int, int> Dim = dom->getDim();
	    for( map<int, int>::const_iterator gDit = Dim.begin();
				gDit != Dim.end();
				gDit++ )
		{
		  if(dom->getThermostat(this->componentid()) == gDit->first){
		    v2 = 0.0;
		    v_aux = 0.0;
		    // in the case of directed velocities
		    v_aux = _v[gDit->second] - _directedVelocity[gDit->second];
		    v2 += v_aux * v_aux;	
		    assert(!isnan(v2)); // catches NaN
		    summv2_1Dim += _m * v2;
		  }
		}	  
	}
}

double Molecule::U_rot() {
	double w[3];
	_q.rotate(_L, w);
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
	return 0.5 * Iw2;
}

void Molecule::calculate_mv2_Iw2(double& summv2, double& summv2_1Dim, double& sumIw2, int dimToThermostat, Domain *dom) {
  
	// with  correction due to directed movement
	double v2 = (_v[0]-_directedVelocity[0])*(_v[0]-_directedVelocity[0]) + (_v[1]-_directedVelocity[1])*(_v[1]-_directedVelocity[1]) + (_v[2]-_directedVelocity[2])*(_v[2]-_directedVelocity[2]);
	summv2 += _m * v2;
	
	// without correction due to directed movement
	//summv2 += _m * v2();
	
	double v2_1Dim = 0.0;
	if(dom->isScaling1Dim(dom->getThermostat(this->componentid()))){ 
	  v2_1Dim = (v(dimToThermostat)-_directedVelocity[dimToThermostat])*(v(dimToThermostat)-_directedVelocity[dimToThermostat]);
	  summv2_1Dim += _m * v2_1Dim;
	}
	double w[3];
	_q.rotate(_L, w);
	double Iw2 = 0.;
	
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
	sumIw2 += Iw2;
}

void Molecule::calculate_mv2_Iw2(double& summv2, double& sumIw2, string directedVelocityType) {
	  
	string Default ("Default");
	string Slab ("Slab");
	string Stress ("Stress");
	string Confinement ("Confinement");
	
	// with  correction due to directed movement
	double v2;
	if(directedVelocityType == Slab)
	  v2 = (_v[0]-_directedVelocitySlab[0])*(_v[0]-_directedVelocitySlab[0]) + (_v[1]-_directedVelocitySlab[1])*(_v[1]-_directedVelocitySlab[1]) + (_v[2]-_directedVelocitySlab[2])*(_v[2]-_directedVelocitySlab[2]);
	else if(directedVelocityType == Stress)
	  v2 = (_v[0]-_directedVelocityStress[0])*(_v[0]-_directedVelocityStress[0]) + (_v[1]-_directedVelocityStress[1])*(_v[1]-_directedVelocityStress[1]) + (_v[2]-_directedVelocityStress[2])*(_v[2]-_directedVelocityStress[2]);
	else if(directedVelocityType == Confinement)
	  v2 = (_v[0]-_directedVelocityConfinement[0])*(_v[0]-_directedVelocityConfinement[0]) + (_v[1]-_directedVelocityConfinement[1])*(_v[1]-_directedVelocityConfinement[1]) + (_v[2]-_directedVelocityConfinement[2])*(_v[2]-_directedVelocityConfinement[2]);
	else if(directedVelocityType == Default)
	  v2 = (_v[0]-_directedVelocity[0])*(_v[0]-_directedVelocity[0]) + (_v[1]-_directedVelocity[1])*(_v[1]-_directedVelocity[1]) + (_v[2]-_directedVelocity[2])*(_v[2]-_directedVelocity[2]);
	else
	  v2 = (_v[0]-_directedVelocity[0])*(_v[0]-_directedVelocity[0]) + (_v[1]-_directedVelocity[1])*(_v[1]-_directedVelocity[1]) + (_v[2]-_directedVelocity[2])*(_v[2]-_directedVelocity[2]);
	summv2 += _m * v2;
	
	// without correction due to directed movement
	//summv2 += _m * v2();
	
	double w[3];
	_q.rotate(_L, w);
	double Iw2 = 0.;
	
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
	sumIw2 += Iw2;
}

void Molecule::calculate_mv2_Iw2(double& summv2, double& sumIw2, double offx, double offy, double offz) {
	double vcx = _v[0] - offx;
	double vcy = _v[1] - offy;
	double vcz = _v[2] - offz;

	summv2 += _m * (vcx*vcx + vcy*vcy + vcz*vcz);

	double w[3];
	_q.rotate(_L, w);
	double Iw2 = 0.;
	for (unsigned short d = 0; d < 3; ++d) {
		w[d] *= _invI[d];
		Iw2 += _I[d] * w[d] * w[d];
	}
	sumIw2 += Iw2;
}

void Molecule::scale_v(double s, double offx, double offy, double offz) {
	this->vsub(offx, offy, offz);
	this->scale_v(s);
	this->vadd(offx, offy, offz);
}

void Molecule::calculateHardyIntersection(const double drmEingang[3], double mjx, double mjy, double mjz, Domain *dom, string stress, string weightingFunc) {
	long int xun=0, yun=0, xun_tot=0, yun_tot=0;
	long int xun2=0, yun2=0;
	unsigned unID=0, unID_tot;
	double xMax, yMax;
	long double deltaX=0.0, deltaY=0.0;
	long double xPlane_min, xPlane_max, yPlane_min, yPlane_max;
	long double lambda;
	long double HardyIntersectionX, HardyIntersectionY, HardyIntersectionZ;
	bool hardyIntersection = true;
	// the input string "stress" is like a flag, it decides whether the output is prepared for the confinement or the whole system
	string Stress ("Stress");
	string Confinement ("Confinement");
	string Linear ("Linear");
	string Pyramide ("Pyramide");

	// Vorzeichen bei drm genau verkehrt....
	double drm[3];
	for (int i = 0; i < 3; i++)
	  drm[i] = (-1)*drmEingang[i];
	
	long double fracX, fracY, fracZ;
	double xPlaneShift = 0.0, yPlaneShift = 0.0;
	xMax = dom->getGlobalLength(0);
	yMax = dom->getGlobalLength(1);
	
	// clearing of all maps
	_HardyIntersectionX.clear();
	_HardyIntersectionY.clear();
	_HardyIntersectionZ.clear();
	_bondFractionUNID.clear();
	
	if (stress == Stress){
	  // data transfer for determination of the control volumes
	  xun_tot = dom->getUniversalNProfileUnits_Stress(0);
	  yun_tot = dom->getUniversalNProfileUnits_Stress(1);
	
	  deltaX = 1/dom->getUniversalInvProfileUnit_Stress(0);
	  deltaY = 1/dom->getUniversalInvProfileUnit_Stress(1);
	
	  // assignement of both interacting molecules to their control volumes
	  xun = floor(_r[0] / deltaX);
	  yun = floor(_r[1] / deltaY);
	  xun2 = floor(mjx / deltaX);
	  yun2 = floor(mjy / deltaY);
	}else if (stress == Confinement){
	  // data transfer for determination of the control volumes
	  xun_tot = dom->getUniversalNProfileUnitsStressConfinement(0);
	  yun_tot = dom->getUniversalNProfileUnitsStressConfinement(1);
	
	  deltaX = 1/dom->getUniversalInvProfileUnitStressConfinement(0);
	  deltaY = 1/dom->getUniversalInvProfileUnitStressConfinement(1);
	  
	  
	  // assignement of both interacting molecules to their control volumes
	  xun = floor((_r[0]-dom->getConfinementEdge(0)) / deltaX);
	  yun = floor((_r[1]-(dom->get_confinementMidPoint(3))) / deltaY);
	  xun2 = floor((mjx-dom->getConfinementEdge(0)) / deltaX);
	  yun2 = floor((mjy-(dom->get_confinementMidPoint(3))) / deltaY);
	  
	  if ((xun < 0 || xun >= xun_tot || yun < 0 || yun >= yun_tot) && (xun2 < 0 || xun2 >= xun_tot || yun2 < 0 || yun2 >= yun_tot)){
	    hardyIntersection = false;
	  }
	  	  
	  // smallest value that is a multiple of (_r[0]-dom->getConfinementEdge(0)) / deltaX    
	  xPlaneShift = dom->getConfinementEdge(0);
	  yPlaneShift = dom->get_confinementMidPoint(3);
	  
	  //cout << " x1 " << _r[0] << " y1 " << _r[1] << " x2 " << mjx << " y2 " << mjy << " xun " << xun << " yun " << yun << " xun2 " << xun2 << " yun2 " << yun2 << " dx " << deltaX << " dy " << deltaY << " unX " << xun_tot << " unY " << yun_tot << " shX " << xPlaneShift << " shY " << yPlaneShift << endl;
	}
	
	unID_tot = unsigned(xun_tot*yun_tot);
	
	// calculation of the coordinates of the intersection points of the bond between molecules with the edges of the control volumes
	// differentiation in four cases
	// ------- case 1 --------- /> 
	if(hardyIntersection && drm[0] >= 0. && drm[1] >= 0.){
	  xPlane_min = xun * deltaX + xPlaneShift;
	  yPlane_min = yun * deltaY + yPlaneShift;
	  xPlane_max = (xun2+1) * deltaX + xPlaneShift;
	  yPlane_max = (yun2+1) * deltaY + yPlaneShift;
	  if(xPlane_min + deltaX != xPlane_max){
	    for (long double i = xPlane_min + deltaX; i < xPlane_max; i += deltaX){
	      lambda = (i - _r[0])/(mjx - _r[0]);
	      fracX = lambda * (mjx - _r[0]);
	      fracY = lambda * (mjy - _r[1]);
	      fracZ = lambda * (mjz - _r[2]);
	      HardyIntersectionX = _r[0] + fracX;
	      HardyIntersectionY = _r[1] + fracY;
	      HardyIntersectionZ = _r[2] + fracZ;	      
	      // checks the existence of the key = lambda
	      if(_HardyIntersectionX.count(lambda) == 0){
		_HardyIntersectionX.insert(std::pair<long double, long double>(lambda, HardyIntersectionX));
		_HardyIntersectionY.insert(std::pair<long double, long double>(lambda, HardyIntersectionY));
		_HardyIntersectionZ.insert(std::pair<long double, long double>(lambda, HardyIntersectionZ));
	      }else{
		lambda += 1.0E-300;
		_HardyIntersectionX.insert(std::pair<long double, long double>(lambda, HardyIntersectionX));
		_HardyIntersectionY.insert(std::pair<long double, long double>(lambda, HardyIntersectionY));
		_HardyIntersectionZ.insert(std::pair<long double, long double>(lambda, HardyIntersectionZ));
	      }
	    }
	  }
	  if(yPlane_min + deltaY != yPlane_max){
	    for (long double i = yPlane_min + deltaY; i < yPlane_max; i += deltaY){
	      lambda = (i - _r[1])/(mjy - _r[1]);
	      fracX = lambda * (mjx - _r[0]);
	      fracY = lambda * (mjy - _r[1]);
	      fracZ = lambda * (mjz - _r[2]);
	      HardyIntersectionX = _r[0] + fracX;
	      HardyIntersectionY = _r[1] + fracY;
	      HardyIntersectionZ = _r[2] + fracZ;
	      // checks the existence of the key = lambda
	      if(_HardyIntersectionX.count(lambda) == 0){
		_HardyIntersectionX.insert(std::pair<long double, long double>(lambda, HardyIntersectionX));
		_HardyIntersectionY.insert(std::pair<long double, long double>(lambda, HardyIntersectionY));
		_HardyIntersectionZ.insert(std::pair<long double, long double>(lambda, HardyIntersectionZ));
	      }else{
		lambda += 1.0E-300;
		_HardyIntersectionX.insert(std::pair<long double, long double>(lambda, HardyIntersectionX));
		_HardyIntersectionY.insert(std::pair<long double, long double>(lambda, HardyIntersectionY));
		_HardyIntersectionZ.insert(std::pair<long double, long double>(lambda, HardyIntersectionZ));
	      }
	    }
	  }
	// ------- case 2 --------- \> 
	}else if(hardyIntersection && drm[0] >= 0. && drm[1] <= 0.){
	  xPlane_min = xun * deltaX + xPlaneShift;
	  yPlane_min = yun2 * deltaY + yPlaneShift;
	  xPlane_max = (xun2+1) * deltaX + xPlaneShift;
	  yPlane_max = (yun+1) * deltaY + yPlaneShift;
	  if(xPlane_min + deltaX != xPlane_max){
	    for (long double i = xPlane_min + deltaX; i < xPlane_max; i += deltaX){
	      lambda = (i - _r[0])/(mjx - _r[0]);
	      fracX = lambda * (mjx - _r[0]);
	      fracY = lambda * (mjy - _r[1]);
	      fracZ = lambda * (mjz - _r[2]);
	      HardyIntersectionX = _r[0] + fracX;
	      HardyIntersectionY = _r[1] + fracY;
	      HardyIntersectionZ = _r[2] + fracZ;
	      // checks the existence of the key = lambda
	      if(_HardyIntersectionX.count(lambda) == 0){
		_HardyIntersectionX.insert(std::pair<long double, long double>(lambda, HardyIntersectionX));
		_HardyIntersectionY.insert(std::pair<long double, long double>(lambda, HardyIntersectionY));
		_HardyIntersectionZ.insert(std::pair<long double, long double>(lambda, HardyIntersectionZ));
	      }else{
		lambda += 1.0E-300;
		_HardyIntersectionX.insert(std::pair<long double, long double>(lambda, HardyIntersectionX));
		_HardyIntersectionY.insert(std::pair<long double, long double>(lambda, HardyIntersectionY));
		_HardyIntersectionZ.insert(std::pair<long double, long double>(lambda, HardyIntersectionZ));
	      }
	    }
	  }
	  if(yPlane_min + deltaY != yPlane_max){
	    for (long double i = yPlane_max - deltaY; i > yPlane_min; i -= deltaY){
	      lambda = (i - _r[1])/(mjy - _r[1]);
	      fracX = lambda * (mjx - _r[0]);
	      fracY = lambda * (mjy - _r[1]);
	      fracZ = lambda * (mjz - _r[2]);
	      HardyIntersectionX = _r[0] + fracX;
	      HardyIntersectionY = _r[1] + fracY;
	      HardyIntersectionZ = _r[2] + fracZ;
	      // checks the existence of the key = lambda
	      if(_HardyIntersectionX.count(lambda) == 0){
		_HardyIntersectionX.insert(std::pair<long double, long double>(lambda, HardyIntersectionX));
		_HardyIntersectionY.insert(std::pair<long double, long double>(lambda, HardyIntersectionY));
		_HardyIntersectionZ.insert(std::pair<long double, long double>(lambda, HardyIntersectionZ));
	      }else{
		lambda += 1.0E-300;
		_HardyIntersectionX.insert(std::pair<long double, long double>(lambda, HardyIntersectionX));
		_HardyIntersectionY.insert(std::pair<long double, long double>(lambda, HardyIntersectionY));
		_HardyIntersectionZ.insert(std::pair<long double, long double>(lambda, HardyIntersectionZ));
	      }
	    }
	  }
	// ------- case 3 --------- </ 
	}else if(hardyIntersection && drm[0] <= 0. && drm[1] <= 0.){
	  xPlane_min = xun2 * deltaX + xPlaneShift;
	  yPlane_min = yun2 * deltaY + yPlaneShift;
	  xPlane_max = (xun+1) * deltaX + xPlaneShift;
	  yPlane_max = (yun+1) * deltaY + yPlaneShift;
	  if(xPlane_min + deltaX != xPlane_max){
	    for (long double i = xPlane_min + deltaX; i < xPlane_max; i += deltaX){
	      lambda = (i - _r[0])/(mjx - _r[0]);
	      fracX = lambda * (mjx - _r[0]);
	      fracY = lambda * (mjy - _r[1]);
	      fracZ = lambda * (mjz - _r[2]);
	      HardyIntersectionX = _r[0] + fracX;
	      HardyIntersectionY = _r[1] + fracY;
	      HardyIntersectionZ = _r[2] + fracZ;
	      // checks the existence of the key = lambda
	      if(_HardyIntersectionX.count(lambda) == 0){
		_HardyIntersectionX.insert(std::pair<long double, long double>(lambda, HardyIntersectionX));
		_HardyIntersectionY.insert(std::pair<long double, long double>(lambda, HardyIntersectionY));
		_HardyIntersectionZ.insert(std::pair<long double, long double>(lambda, HardyIntersectionZ));
	      }else{
		lambda += 1.0E-300;
		_HardyIntersectionX.insert(std::pair<long double, long double>(lambda, HardyIntersectionX));
		_HardyIntersectionY.insert(std::pair<long double, long double>(lambda, HardyIntersectionY));
		_HardyIntersectionZ.insert(std::pair<long double, long double>(lambda, HardyIntersectionZ));
	      }
	    }
	  }
	  if(yPlane_min + deltaY != yPlane_max){
	    for (long double i = yPlane_min + deltaY; i < yPlane_max; i += deltaY){
	      lambda = (i - _r[1])/(mjy - _r[1]);
	      fracX = lambda * (mjx - _r[0]);
	      fracY = lambda * (mjy - _r[1]);
	      fracZ = lambda * (mjz - _r[2]);
	      HardyIntersectionX = _r[0] + fracX;
	      HardyIntersectionY = _r[1] + fracY;
	      HardyIntersectionZ = _r[2] + fracZ;
	      // checks the existence of the key = lambda
	      if(_HardyIntersectionX.count(lambda) == 0){
		_HardyIntersectionX.insert(std::pair<long double, long double>(lambda, HardyIntersectionX));
		_HardyIntersectionY.insert(std::pair<long double, long double>(lambda, HardyIntersectionY));
		_HardyIntersectionZ.insert(std::pair<long double, long double>(lambda, HardyIntersectionZ));
	      }else{
		lambda += 1.0E-300;
		_HardyIntersectionX.insert(std::pair<long double, long double>(lambda, HardyIntersectionX));
		_HardyIntersectionY.insert(std::pair<long double, long double>(lambda, HardyIntersectionY));
		_HardyIntersectionZ.insert(std::pair<long double, long double>(lambda, HardyIntersectionZ));
	      }
	    }
	  }
	// ------- case 4 --------- <\ // 
	}else if(hardyIntersection && drm[0] <= 0. && drm[1] >= 0.){
	  xPlane_min = xun2 * deltaX + xPlaneShift;
	  yPlane_min = yun * deltaY + yPlaneShift;
	  xPlane_max = (xun+1) * deltaX + xPlaneShift;
	  yPlane_max = (yun2+1) * deltaY + yPlaneShift;
	  if(xPlane_min + deltaX != xPlane_max){
	    for (long double i = xPlane_min + deltaX; i < xPlane_max; i += deltaX){
	      lambda = (i - _r[0])/(mjx - _r[0]);
	      fracX = lambda * (mjx - _r[0]);
	      fracY = lambda * (mjy - _r[1]);
	      fracZ = lambda * (mjz - _r[2]);
	      HardyIntersectionX = _r[0] + fracX;
	      HardyIntersectionY = _r[1] + fracY;
	      HardyIntersectionZ = _r[2] + fracZ;
	      // checks the existence of the key = lambda
	      if(_HardyIntersectionX.count(lambda) == 0){
		_HardyIntersectionX.insert(std::pair<long double, long double>(lambda, HardyIntersectionX));
		_HardyIntersectionY.insert(std::pair<long double, long double>(lambda, HardyIntersectionY));
		_HardyIntersectionZ.insert(std::pair<long double, long double>(lambda, HardyIntersectionZ));
	      }else{
		lambda += 1.0E-300;
		_HardyIntersectionX.insert(std::pair<long double, long double>(lambda, HardyIntersectionX));
		_HardyIntersectionY.insert(std::pair<long double, long double>(lambda, HardyIntersectionY));
		_HardyIntersectionZ.insert(std::pair<long double, long double>(lambda, HardyIntersectionZ));
	      }
	    }
	  }
	  if(yPlane_min + deltaY != yPlane_max){
	    for (long double i = yPlane_max - deltaY; i > yPlane_min; i -= deltaY){
	      lambda = (i - _r[1])/(mjy - _r[1]);
	      fracX = lambda * (mjx - _r[0]);
	      fracY = lambda * (mjy - _r[1]);
	      fracZ = lambda * (mjz - _r[2]);
	      HardyIntersectionX = _r[0] + fracX;
	      HardyIntersectionY = _r[1] + fracY;
	      HardyIntersectionZ = _r[2] + fracZ;
	      // checks the existence of the key = lambda
	      if(_HardyIntersectionX.count(lambda) == 0){
		_HardyIntersectionX.insert(std::pair<long double, long double>(lambda, HardyIntersectionX));
		_HardyIntersectionY.insert(std::pair<long double, long double>(lambda, HardyIntersectionY));
		_HardyIntersectionZ.insert(std::pair<long double, long double>(lambda, HardyIntersectionZ));
	      }else{
		lambda += 1.0E-300;
		_HardyIntersectionX.insert(std::pair<long double, long double>(lambda, HardyIntersectionX));
		_HardyIntersectionY.insert(std::pair<long double, long double>(lambda, HardyIntersectionY));
		_HardyIntersectionZ.insert(std::pair<long double, long double>(lambda, HardyIntersectionZ));
	      }
	    }
	  }
	}
		
	// first and last insert are the particle centre of mass themselves
	_HardyIntersectionX[0.0] = _r[0];
	_HardyIntersectionY[0.0] = _r[1];
	_HardyIntersectionZ[0.0] = _r[2];
	
	_HardyIntersectionX[1.0] = mjx;
	_HardyIntersectionY[1.0] = mjy;
	_HardyIntersectionZ[1.0] = mjz;
	
	
	long double halfDeltaX = 0, halfDeltaY = 0, deltaXMol = 0, deltaYMol = 0, deltaXCentre = 0, deltaYCentre = 0;
	int sign[2]= {1,1};
	// calculate contribution of 2D weighting function being pyramidal
	if(hardyIntersection && weightingFunc == Pyramide){
	  halfDeltaX = 0.5*deltaX;
	  halfDeltaY = 0.5*deltaY;
	  deltaXMol = mjx-_r[0];
	  deltaYMol = mjy-_r[1];
	  sign[0] = signumFunction(deltaXMol);  
	  sign[1] = signumFunction(deltaYMol);
	}
	
	long double previousValue = 0.0;
	// assign the _bondFraction[a] to a control volume [b]
      if(hardyIntersection){
	for(std::map<long double, long double>::iterator it=_HardyIntersectionX.begin(); it!=_HardyIntersectionX.end(); ++it){
	  std::map<long double, int> _root;
	  
	  if(it->first == 0.0){
	    previousValue = it->first;
	    continue;
	  }
	  long double xMean, yMean;
	  xMean = 0.5*(_HardyIntersectionX[previousValue]+_HardyIntersectionX[it->first]);
	  yMean = 0.5*(_HardyIntersectionY[previousValue]+_HardyIntersectionY[it->first]); 
	  
	  // transformation of coordinates from halo cells to in-box-cells
	  if(xMean < 0.0)
	    xMean += xMax;
	  if(xMean > xMax)
	    xMean -= xMax;
	  if(yMean < 0.0)
	    yMean += yMax;
	  if(yMean > yMax)
	    yMean -= yMax;
	  	  
	  // the bond fraction is assigned to the unID that is between the points _HardyIntersectionX[i] and _HardyIntersectionX[i+1]
	  if (stress == Stress){
	    xun = floor(xMean / deltaX);
	    yun = floor(yMean / deltaY);
	  }else if (stress == Confinement){
	    xun = floor((xMean-dom->getConfinementEdge(0)) / deltaX);
	    yun = floor((yMean-(dom->get_confinementMidPoint(3))) / deltaY);
	  }
	  
	  if (xun >= 0 && yun >= 0 && xun < xun_tot && yun < yun_tot){
	    unID = xun * yun_tot  + yun;
	    
	    if (weightingFunc == Linear && unID >= 0 && unID < unID_tot)
	      _bondFractionUNID[unID] = (it->first) - previousValue;
	     
	    // calculate contribution of 2D weighting function being pyramidal
	    if(weightingFunc == Pyramide  && unID >= 0 && unID < unID_tot){
	      _bondFractionUNID[unID] = 0.0;
	      if (stress == Stress){
		  deltaXCentre = _r[0]-(xun*deltaX+halfDeltaX);
		  deltaYCentre = _r[1]-(yun*deltaY+halfDeltaY);
		  if(abs(deltaXCentre) > abs(xMax/2)){
		    if(_r[0] > mjx){
		      if(deltaXCentre < 0) 
			deltaXCentre = xMax+deltaXCentre;
		      else if(deltaXCentre > 0)
			deltaXCentre = xMax-deltaXCentre;
		    }else{
		      if(deltaXCentre < 0) 
			deltaXCentre = (-1)*(xMax+deltaXCentre);
		      else if(deltaXCentre > 0)
			deltaXCentre = (-1)*(xMax-deltaXCentre);
		    } 
		  }
		  if(abs(deltaYCentre) > abs(yMax/2)){
		    if(_r[1] > mjy){
		      if(deltaYCentre < 0) 
			deltaYCentre = yMax+deltaYCentre;
		      else if(deltaYCentre > 0)
			deltaYCentre = yMax-deltaYCentre;
		    }else{
		      if(deltaYCentre < 0) 
			deltaYCentre = (-1)*(yMax+deltaYCentre);
		      else if(deltaYCentre > 0)
			deltaYCentre = (-1)*(yMax-deltaYCentre);
		    }
		  }
	      }else if (stress == Confinement){
		  deltaXCentre = _r[0]-(xun*deltaX+halfDeltaX+dom->getConfinementEdge(0));
		  deltaYCentre = _r[1]-(yun*deltaY+halfDeltaY+dom->get_confinementMidPoint(3));
		  if(abs(deltaXCentre) > abs(xMax/2)){
		    if(_r[0] > mjx){
		      if(deltaXCentre < 0) 
			deltaXCentre = xMax+deltaXCentre;
		      else if(deltaXCentre > 0)
			deltaXCentre = xMax-deltaXCentre;
		    }else{
		      if(deltaXCentre < 0) 
			deltaXCentre = (-1)*(xMax+deltaXCentre);
		      else if(deltaXCentre > 0)
			deltaXCentre = (-1)*(xMax-deltaXCentre);
		    } 
		  }
		  if(abs(deltaYCentre) > abs(yMax/2)){
		    if(_r[1] > mjy){
		      if(deltaYCentre < 0) 
			deltaYCentre = yMax+deltaYCentre;
		      else if(deltaYCentre > 0)
			deltaYCentre = yMax-deltaYCentre;
		    }else{
		      if(deltaYCentre < 0) 
			deltaYCentre = (-1)*(yMax+deltaYCentre);
		      else if(deltaYCentre > 0)
			deltaYCentre = (-1)*(yMax-deltaYCentre);
		    }
		  }
	      }
	        
	      _root.insert(std::pair<long double, int>(previousValue,1));
	      _root.insert(std::pair<long double, int>(it->first,1));
	      
	      // calculate roots of the function f = lambda*deltaXMol + deltaXCentre and g = lambda*deltaYMol + deltaYCentre due to integration over absolute value of both
	      if((-1)*deltaXCentre/deltaXMol > previousValue && (-1)*deltaXCentre/deltaXMol < it->first){
		  _root.insert(std::pair<long double, int>((-1)*deltaXCentre/deltaXMol,sign[0]));}
	      if((-1)*deltaYCentre/deltaYMol > previousValue && (-1)*deltaYCentre/deltaYMol < it->first){
		  _root.insert(std::pair<long double, int>((-1)*deltaYCentre/deltaYMol,sign[1]));}  
	      
	      double signumF = 1.0, signumG = 1.0;
	      
	      long double previousValueRoot = 0;
	      for(std::map<long double, int>::iterator it2=_root.begin(); it2!=_root.end(); ++it2){
		if(it2->first == previousValue){
		  previousValueRoot = previousValue;
		  continue;
		}
		
		signumF = (previousValueRoot+(it2->first))*0.5*deltaXMol+deltaXCentre;
		signumG = (previousValueRoot+(it2->first))*0.5*deltaYMol+deltaYCentre;
		sign[0] = signumFunction(signumF);
		sign[1] = signumFunction(signumG);

		  _bondFractionUNID[unID] += bondFunctionPyramide(previousValueRoot, (it2->first), sign[0], sign[1], halfDeltaX, halfDeltaY, deltaXMol, deltaYMol, deltaXCentre, deltaYCentre);

		previousValueRoot = it2->first;
	      }
	    }  
	  }
	  previousValue = it->first;
	}
      }
}
 
long double Molecule::bondFunctionPyramide(long double n0, long double n1, int signF, int signG, long double halfDeltaX, long double halfDeltaY, long double deltaXMol, long double deltaYMol, long double deltaXCentre, long double deltaYCentre){
     long double bondFunc, bondFunc0, bondFunc1;
     
     //bondFunc1 = 1/(6*halfDeltaX*halfDeltaX*halfDeltaY*halfDeltaY)*(n1*(3*halfDeltaX*(2*halfDeltaY-(2*deltaYCentre+deltaYMol*n1)*signG)+signF*((-3)*halfDeltaY*(2*deltaXCentre+deltaXMol*n1)+(6*deltaXCentre*deltaYCentre+3*deltaXCentre*deltaYMol*n1+3*deltaXMol*deltaYCentre*n1+2*deltaXMol*deltaYMol*n1*n1)*signG)));
     //bondFunc0 = 1/(6*halfDeltaX*halfDeltaX*halfDeltaY*halfDeltaY)*(n0*(3*halfDeltaX*(2*halfDeltaY-(2*deltaYCentre+deltaYMol*n0)*signG)+signF*((-3)*halfDeltaY*(2*deltaXCentre+deltaXMol*n0)+(6*deltaXCentre*deltaYCentre+3*deltaXCentre*deltaYMol*n0+3*deltaXMol*deltaYCentre*n0+2*deltaXMol*deltaYMol*n0*n0)*signG)));
     bondFunc1 = 1/(halfDeltaX*halfDeltaY)*(n1-1/(2*halfDeltaX)*n1*(deltaXMol*n1+2*deltaXCentre)*signF-1/(2*halfDeltaY)*n1*(deltaYMol*n1+2*deltaYCentre)*signG
		    +1/(halfDeltaX*halfDeltaY*6)*n1*signF*signG*(n1*deltaXMol*(2*deltaYMol*n1+3*deltaYCentre)+3*deltaXCentre*(deltaYMol*n1+2*deltaYCentre)));
     bondFunc0 = 1/(halfDeltaX*halfDeltaY)*(n0-1/(2*halfDeltaX)*n0*(deltaXMol*n0+2*deltaXCentre)*signF-1/(2*halfDeltaY)*n0*(deltaYMol*n0+2*deltaYCentre)*signG
		    +1/(halfDeltaX*halfDeltaY*6)*n0*signF*signG*(n0*deltaXMol*(2*deltaYMol*n0+3*deltaYCentre)+3*deltaXCentre*(deltaYMol*n0+2*deltaYCentre)));
     
     bondFunc = bondFunc1 - bondFunc0;

     return bondFunc;
}

double Molecule::weightingFunctionPyramide(unsigned xun, unsigned yun, double deltaX, double deltaY, double xStart, double yStart){
     double CentreX, CentreY, weighting;
     CentreX = xun*deltaX + 0.5*deltaX + xStart;
     CentreY = yun*deltaY + 0.5*deltaY + yStart;
     if(abs(CentreX-_r[0]) < 0.5*deltaX && abs(CentreY-_r[1]) < 0.5*deltaY)
	weighting = 4/(deltaX*deltaY)*(1-2/deltaX*abs(CentreX-_r[0]))*(1-2/deltaY*abs(CentreY-_r[1]));
     else
        weighting = 0.0;
     return weighting;
}

void Molecule::write(ostream& ostrm) const {
	ostrm << _id << "\t" << (_component->ID() + 1) << "\t"
	      << _r[0] << " " << _r[1] << " " << _r[2] << "\t"
	      << _v[0] << " " << _v[1] << " " << _v[2] << "\t"
	      << _q.qw() << " " << _q.qx() << " " << _q.qy() << " " << _q.qz() << "\t"
	      << _L[0] << " " << _L[1] << " " << _L[2] << "\t"
	      << endl;
}

void Molecule::addTersoffNeighbour(Molecule* m, bool pairType) {
	// this->_Tersoff_neighbours.insert(pair<Molecule*, bool>(m, (pairType > 0)));
	for (int j = 0; j < _numTersoffNeighbours; j++) {
		if (m->_id == _Tersoff_neighbours_first[j]->id()) {
			this->_Tersoff_neighbours_first[j] = m;
			this->_Tersoff_neighbours_second[j] = pairType;
			return;
		}
	}

	this->_Tersoff_neighbours_first[_numTersoffNeighbours] = m;
	this->_Tersoff_neighbours_second[_numTersoffNeighbours] = pairType;
	this->_numTersoffNeighbours++;
	if (_numTersoffNeighbours > MAX_TERSOFF_NEIGHBOURS) {
		global_log->error() << "Tersoff neighbour list overflow: Molecule " << m->_id << " has more than " << MAX_TERSOFF_NEIGHBOURS << " Tersoff neighbours." << endl;
		exit(1);
	}
}

double Molecule::tersoffParameters(double params[15]) //returns delta_r
{
	const Tersoff* t = &_component->tersoff()[0];
	params[ 0] = t->R();
	params[ 1] = t->S();
	params[ 2] = t->h();
	params[ 3] = t->cSquare();
	params[ 4] = t->dSquare();
	params[ 5] = t->A();
	params[ 6] = t->minusLambda();
	params[ 7] = t->minusMu();
	params[ 8] = t->beta();
	params[ 9] = t->n();
	params[10] = M_PI / (t->S() - t->R());
	params[11] = 1.0 + t->cSquare() / t->dSquare();
	params[12] = t->S() * t->S();
	params[13] = -(t->B());
	params[14] = -0.5 / t->n();

	return 0.000001 * (t->S() - t->R());
}

// private functions
// these are only used when compiling molecule.cpp and therefore might be inlined without any problems

inline void Molecule::setupCache() {
	assert(_component != NULL);

	_m = _component->m();
	_I[0] = _component->I11();
	_I[1] = _component->I22();
	_I[2] = _component->I33();
	for (unsigned short d = 0; d < 3; ++d) {
		if (_I[d] != 0.)
			_invI[d] = 1. / _I[d];
		else
			_invI[d] = 0.;
	}

	int numsites = _component->numSites();
	_sites_d = new double[3*numsites];
	assert(_sites_d);
	_ljcenters_d = &(_sites_d[0]);
	_charges_d = &(_ljcenters_d[3*numLJcenters()]);
	_dipoles_d = &(_charges_d[3*numCharges()]);
	_quadrupoles_d = &(_dipoles_d[3*numDipoles()]);
	_tersoff_d = &(_quadrupoles_d[3*numQuadrupoles()]);

	int numorientedsites = _component->numOrientedSites();
	_osites_e = new double[3*numorientedsites];
	assert(_osites_e);
	_dipoles_e = &(_osites_e[0]);
	_quadrupoles_e = &(_dipoles_e[3*numDipoles()]);

	_sites_F = new double[3*numsites];
	//_springSites_F = new double[3*numsites];
	assert(_sites_F);
	_ljcenters_F = &(_sites_F[0]);
	_charges_F = &(_ljcenters_F[3*numLJcenters()]);
	_dipoles_F = &(_charges_F[3*numCharges()]);
	_quadrupoles_F = &(_dipoles_F[3*numDipoles()]);
	_tersoff_F = &(_quadrupoles_F[3*numQuadrupoles()]);
	//_spring_F = &(_springSites_F[0]);

	this->clearFM();
}

void Molecule::clearFM() {
	int numSites = _component->numSites();
	for (int i = 0; i < 3*numSites; i++) {
		_sites_F[i] = 0.;
		//_springSites_F[i] = 0.;
	}
	_F[0] = _F[1] = _F[2] = 0.;
	_M[0] = _M[1] = _M[2] = 0.;
}

void Molecule::calcFM() {
	//_M[0] = _M[1] = _M[2] = 0.;
	unsigned int ns = numSites();
	for (unsigned int si = 0; si < ns; ++si) {
		const double* Fsite = site_F(si);
		const double* dsite = site_d(si);
#ifndef NDEBUG
		/*
		 * catches NaN assignments
		 */
		for (int d = 0; d < 3; d++) {
			if (isnan(dsite[d])) {
				global_log->error() << "Severe dsite[" << d << "] error for site " << si << " of m" << _id << endl;
				assert(false);
			}
			if (isnan(Fsite[d])) {
				global_log->error() << "Severe Fsite[" << d << "] error for site " << si << " of m" << _id << endl;
				assert(false);
			}
		}
#endif

		Fadd(Fsite);
		_M[0] += dsite[1] * Fsite[2] - dsite[2] * Fsite[1];
		_M[1] += dsite[2] * Fsite[0] - dsite[0] * Fsite[2];
		_M[2] += dsite[0] * Fsite[1] - dsite[1] * Fsite[0];
	}
}


/**
 * catches NaN values and missing data
 *
 * @note Use isnan from cmath to check for nan.
 * If that's not available (C99), compare the value with itself. If the value
 * is NaN, the comparison will evaluate to false (according to IEEE754 spec.)
 */
void Molecule::check(unsigned long id) {
#ifndef NDEBUG
	assert(_id == id);
	assert(_m > 0.0);
	for (int d = 0; d < 3; d++) {
		assert(!isnan(_r[d]));
		assert(!isnan(_v[d]));
		assert(!isnan(_L[d]));
		assert(!isnan(_F[d]));
		assert(!isnan(_M[d]));
		assert(!isnan(_I[d]));
		assert(!isnan(_invI[d]));
	}
#endif
}

bool Molecule::isLessThan(const Molecule& m2) const {
	if (_r[2] < m2.r(2))
		return true;
	else if (_r[2] > m2.r(2))
		return false;
	else {
		if (_r[1] < m2.r(1))
			return true;
		else if (_r[1] > m2.r(1))
			return false;
		else {
			if (_r[0] < m2.r(0))
				return true;
			else if (_r[0] > m2.r(0))
				return false;
			else {
				global_log->error() << "LinkedCells::isFirstParticle: both Particles have the same position" << endl;
				exit(1);
			}
		}
	}
	return false; /* Silence warnings about missing return statement */
}


std::ostream& operator<<( std::ostream& os, const Molecule& m ) {
	os << "ID: " << m.id() << "\n";
	os << "r:  (" << m.r(0) << ", " << m.r(1) << ", " << m.r(2) << ")\n" ;
	os << "v:  (" << m.v(0) << ", " << m.v(1) << ", " << m.v(2) << ")\n" ;
	os << "F:  (" << m.F(0) << ", " << m.F(1) << ", " << m.F(2) << ")\n" ;
	os << "q:  [[" << m.q().qw() << ", " << m.q().qx() << "], [" << m.q().qy() << ", " << m.q().qz()<< "]]\n" ;
	os << "w:  (" << m.D(0) << ", " << m.D(1) << ", " << m.D(2) << ")" ;
	return os;
}


unsigned long Molecule::totalMemsize() const {
	unsigned long size = sizeof (*this);

	//_sites_d
	size += sizeof(double) * _component->numSites()* 3;
	// site orientation _osites_e
	size += sizeof(double) * _component->numOrientedSites() * 4;
	// site Forces _sites_F
	// row order: Fx1,Fy1,Fz1,Fx2,Fy2,Fz2,...
	size += sizeof(double) * _component->numSites() * 3;

	return size;
}
