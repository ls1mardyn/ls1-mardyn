/*
 * M. T. Horsch, LS1/Mardyn project moderated by Martin Bernreuther
 * (C)2010 GNU General Public License
 */
#include "PressureGradient.h"

#include <cmath>
#include <iostream>

#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"

using namespace Log;
using namespace std;

PressureGradient::PressureGradient(int rank) {
	this->_localRank = rank;
	this->_universalConstantAccelerationTimesteps = 0;
	if(!rank){
		for(unsigned short int d=0; d < 3; d++)
			this->_globalVelocitySum[d] = map<unsigned int, long double>();
	}
	this->_universalConstantTau = true;
	this->_universalZetaFlow = 0.0;
	this->_universalTauPrime = 0.0;
	this->_shearRampTime = 1;
	this->_doApplyShearForce = false;
}

void PressureGradient::specifyComponentSet(unsigned int cosetid, double v[3], double tau, double ainit[3], double timestep)
{
	this->_localN[cosetid] = 0;
	for(unsigned short int d = 0; d < 3; d++)
	{
		this->_universalAdditionalAcceleration[d][cosetid] = ainit[d];
		this->_localVelocitySum[d][cosetid] = 0.0;
	}
	this->_universalTau[cosetid] = tau;
	if(this->_localRank == 0)
	{
		this->_globalN[cosetid] = 0;
		for(unsigned short int d = 0; d < 3; d++)
		{
			this->_globalTargetVelocity[d][cosetid] = v[d];
			this->_globalVelocitySum[d][cosetid] = 0.0;
		}
		this->_globalVelocityQueuelength[cosetid] = (unsigned int)ceil(
			(this->_universalTauPrime == 0.0)?
				sqrt(this->_universalTau[cosetid] / (timestep*this->_universalConstantAccelerationTimesteps))
				: this->_universalTauPrime / (timestep*this->_universalConstantAccelerationTimesteps)
		);
	}
	// _globalTargetVelocity is distributed to all processes
	if(!isAcceleratingUniformly())
	{
	 for(unsigned short int d = 0; d < 3; d++)
		{
			this->_globalTargetVelocity[d][cosetid] = v[d];
		}
	}
}

bool PressureGradient::isAcceleratingInstantaneously(unsigned numberOfComp)
{
	  double vel = 0.0;
	  bool acc = false;
	  for (unsigned cid = 0; cid < numberOfComp; cid++)
	    for (int d = 0; d < 3; d++)
	      vel += abs(this->getGlobalTargetVelocity(d, cid));
	  
	  if (vel != 0)
	    acc = true;
	    
	  return acc;
}

/*
 * diese Version beschleunigt in alle Raumrichtungen
 */
void PressureGradient::determineAdditionalAcceleration
(
 DomainDecompBase* domainDecomp,
 ParticleContainer* molCont, double dtConstantAcc )
{
	for( map<unsigned int, double>::iterator uAAit = _universalAdditionalAcceleration[0].begin();
			uAAit != _universalAdditionalAcceleration[0].end();
			uAAit++ )
	{
		this->_localN[uAAit->first] = 0;
		for(unsigned short int d = 0; d < 3; d++)
			this->_localVelocitySum[d][uAAit->first] = 0.0;
	}
	for(Molecule* thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
	{
		unsigned int cid = thismol->componentid();
		map<unsigned int, unsigned int>::iterator uCSIDit = this->_universalComponentSetID.find(cid);
		if(uCSIDit == _universalComponentSetID.end()) continue;
		unsigned cosetid = uCSIDit->second;
		this->_localN[cosetid]++;
		for(unsigned short int d = 0; d < 3; d++)
			this->_localVelocitySum[d][cosetid] += thismol->v(d);
	}

	// domainDecomp->collectCosetVelocity(&_localN, _localVelocitySum, &_globalN, _globalVelocitySum);
	//
	domainDecomp->collCommInit( 4 * _localN.size() );
	for( map<unsigned int, unsigned int long>::iterator lNit = _localN.begin(); lNit != _localN.end(); lNit++ )
	{
		domainDecomp->collCommAppendUnsLong(lNit->second);
		for(int d = 0; d < 3; d++)
			domainDecomp->collCommAppendDouble( this->_localVelocitySum[d][lNit->first]);
	}
	domainDecomp->collCommAllreduceSum();
	for( map<unsigned int, unsigned long>::iterator lNit = _localN.begin(); lNit != _localN.end(); lNit++ )
	{
		_globalN[lNit->first] = domainDecomp->collCommGetUnsLong();
		for(int d = 0; d < 3; d++)
			_globalVelocitySum[d][lNit->first] = domainDecomp->collCommGetDouble();
	}
	domainDecomp->collCommFinalize();

	map<unsigned int, long double>::iterator gVSit;
	if(!this->_localRank)
	{
		for(gVSit = _globalVelocitySum[0].begin(); gVSit != _globalVelocitySum[0].end(); gVSit++)
		{
#ifndef NDEBUG
			global_log->debug() << "required entries in velocity queue: " << _globalVelocityQueuelength[gVSit->first] << endl;
			global_log->debug() << "entries in velocity queue: " << _globalPriorVelocitySums[0][gVSit->first].size() << endl;
#endif
			for(unsigned short int d = 0; d < 3; d++)
			{
				while(_globalPriorVelocitySums[d][gVSit->first].size() < _globalVelocityQueuelength[gVSit->first])
					_globalPriorVelocitySums[d][gVSit->first].push_back(_globalVelocitySum[d][gVSit->first]);
			}
		}
	}

	if(!this->_localRank)
	{
		for(gVSit = _globalVelocitySum[0].begin(); gVSit != _globalVelocitySum[0].end(); gVSit++)
		{
			double invgN = 1.0 / this->_globalN[gVSit->first];
			double invgtau = 1.0 / this->_universalTau[gVSit->first];
			double invgtau2 = invgtau * invgtau;
			double previousVelocity[3];
			for(unsigned short int d = 0; d < 3; d++)
			{
				previousVelocity[d] = invgN * _globalPriorVelocitySums[d][gVSit->first].front();
				this->_globalPriorVelocitySums[d][gVSit->first].pop_front();
				this->_universalAdditionalAcceleration[d][gVSit->first]
					+= dtConstantAcc * invgtau2 *
					( this->_globalTargetVelocity[d][gVSit->first]
					  - 2.0*this->_globalVelocitySum[d][gVSit->first]*invgN
					  + previousVelocity[d] );
			}
#ifndef NDEBUG
			global_log->debug() << "accelerator no. " << gVSit->first 
				<< "previous vz: " << previousVelocity[2] 
				<< "current vz: " << _globalVelocitySum[2][gVSit->first]*invgN << endl;
#endif
		}
	}

	domainDecomp->collCommInit(7*this->_localN.size());
	for( map<unsigned int, unsigned long>::iterator lNit = _localN.begin();
			lNit != this->_localN.end();
			lNit++ )
	{
		domainDecomp->collCommAppendUnsLong(_globalN[lNit->first]);
		for(int d=0; d < 3; d++)
		{
			domainDecomp->collCommAppendDouble(
					this->_universalAdditionalAcceleration[d][lNit->first]
					);
			domainDecomp->collCommAppendDouble(
					this->_globalVelocitySum[d][lNit->first]
					);
		}
	}
	domainDecomp->collCommBroadcast();
	for( map<unsigned int, unsigned long>::iterator lNit = _localN.begin();
			lNit != this->_localN.end();
			lNit++ )
	{
		_globalN[lNit->first] = domainDecomp->collCommGetUnsLong();
		for(int d=0; d < 3; d++)
		{
			this->_universalAdditionalAcceleration[d][lNit->first]
				= domainDecomp->collCommGetDouble();
			this->_globalVelocitySum[d][lNit->first]
				= domainDecomp->collCommGetDouble();
		}
	}
	domainDecomp->collCommFinalize();
}

double PressureGradient::getDirectedVelocity(unsigned int cosetid, unsigned short int d)
{
	if(!this->_localRank) 
		return this->_globalVelocitySum[d][cosetid] / this->_globalN[cosetid];
	else return 0.0;
}

double PressureGradient::getUniformAcceleration(unsigned int cosetid)
{
	double aa = 0.0;
	for(unsigned short int d = 0; d < 3; d++)
		aa += _universalAdditionalAcceleration[d][cosetid] * _universalAdditionalAcceleration[d][cosetid];
	return sqrt(aa);
}
double PressureGradient::getUniformAcceleration(unsigned int cosetid, unsigned short int d)
{
	return this->_universalAdditionalAcceleration[d][cosetid];
}

void PressureGradient::prepare_getMissingVelocity(DomainDecompBase* domainDecomp, ParticleContainer* molCont, unsigned int cid, unsigned numberOfComp, unsigned directedVelTime)
{
        
	  this->_localN[cid] = 0;
          this->_globalN[cid] = 0;  
	  for(unsigned short int d = 0; d < 3; d++){
		this->_localVelocitySum[d][cid] = 0.0;
		this->_globalVelocitySum[d][cid] = 0.0;
          }

	//TEST this->_globalPriorAccVelocitySums[0][cid].size() OR this->_globalPriorAccVelocitySums[0].size() ???
        if(!this->_localRank && _globalPriorAccVelocitySums[0][cid].size() == 0){
            for(unsigned short int d = 0; d < 3; d++)
                this->_averagedAccVelocitySum[d][cid] = 0.0;
	    this->_averagedAccN[cid] = 0;
	}

	for(Molecule* thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
	{
	  if(thismol->componentid() == cid){
		this->_localN[cid]++;
		for(unsigned short int d = 0; d < 3; d++)
			this->_localVelocitySum[d][cid] += thismol->v(d);
	  }
	}

	domainDecomp->collCommInit(4);
        domainDecomp->collCommAppendUnsLong(this->_localN[cid]);
	for(int d = 0; d < 3; d++)
		domainDecomp->collCommAppendDouble( this->_localVelocitySum[d][cid]);
	domainDecomp->collCommAllreduceSum();
        _globalN[cid] = domainDecomp->collCommGetUnsLong();
        for(int d = 0; d < 3; d++)
                _globalVelocitySum[d][cid] = domainDecomp->collCommGetDouble();
	
	domainDecomp->collCommFinalize();

        for(int d = 0; d < 3; d++){
	  if(this->_globalN[cid] != 0)
	    this->_directedAccVel[d][cid] = this->_globalVelocitySum[d][cid]/this->_globalN[cid];
	  else
	    this->_directedAccVel[d][cid] = 0.0;
        }

	for(int d = 0; d < 3; d++){
                this->_globalPriorAccVelocitySums[d][cid].push_back(this->_globalVelocitySum[d][cid]);
                this->_averagedAccVelocitySum[d][cid] += this->_globalVelocitySum[d][cid]; 
        }
        this->_globalPriorAccN[cid].push_back(this->_globalN[cid]);
        this->_averagedAccN[cid] += this->_globalN[cid];
	
        for(int d = 0; d < 3; d++){
                if(this->_averagedAccN[cid] != 0)
                    this->_directedAccVelAverage[d][cid] = this->_averagedAccVelocitySum[d][cid]/this->_averagedAccN[cid];
                else
                    this->_directedAccVelAverage[d][cid] = 0.0;
        }
        //TEST this->_globalPriorAccVelocitySums[0][cid].size() OR this->_globalPriorAccVelocitySums[cid].size() ???
        if(this->_globalPriorAccVelocitySums[0][cid].size() == directedVelTime){
            for(int d = 0; d < 3; d++){
                this->_averagedAccVelocitySum[d][cid] -= this->_globalPriorAccVelocitySums[d][cid].front();
                this->_globalPriorAccVelocitySums[d][cid].pop_front();
            }
            this->_averagedAccN[cid] -= this->_globalPriorAccN[cid].front();
            this->_globalPriorAccN[cid].pop_front();
	}

	domainDecomp->collCommInit(3);
        for(int d = 0; d < 3; d++)
		domainDecomp->collCommAppendDouble(this->_directedAccVelAverage[d][cid]);
	domainDecomp->collCommBroadcast();
        for(int d = 0; d < 3; d++)
		this->_directedAccVelAverage[d][cid] = domainDecomp->collCommGetDouble();
	domainDecomp->collCommFinalize();
}

double PressureGradient::getMissingVelocity(unsigned int cid, unsigned short int d, unsigned long simstep, unsigned long init)
{
	double v_directed = this->_directedAccVel[d][cid];
	double v_missing = 0.0;
	double rampTime = 50000*_globalTargetVelocity[d][cid];
	
	double slowAcceleration = 1.0;
	// slowAcceleration = [0;1]
	if((simstep-init) > 0 && (simstep-init) < rampTime)
			slowAcceleration = (double)(simstep-init)/rampTime;
	if(slowAcceleration > 1.0)
			slowAcceleration = 1.0;
	
	if (this->_globalTargetVelocity[d][cid] != 0.0 && this->_globalTargetVelocity[d][cid] > 1e-5)
	  v_missing = slowAcceleration * this->_globalTargetVelocity[d][cid] - 2*v_directed + this->_directedAccVelAverage[d][cid];
	else if (this->_globalTargetVelocity[d][cid] != 0.0 && this->_globalTargetVelocity[d][cid] <= 1e-5)
	  v_missing = (-1)*v_directed; 
	
	this->_currentGlobalTargetVelocity[d][cid] = slowAcceleration * this->_globalTargetVelocity[d][cid];
	
	return v_missing;
}

double* PressureGradient::getTargetVelocity(unsigned int set)
{
	double* retv = new double[3];
	for(int d=0; d < 3; d++)
		retv[d] = this->_globalTargetVelocity[d][set];
	return retv;
}

double* PressureGradient::getAdditionalAcceleration(unsigned int set)
{
	double* retv = new double[3];
	for(int d=0; d < 3; d++)
		retv[d] = this->_universalAdditionalAcceleration[d][set];
	return retv;
}

void PressureGradient::setupShearRate(double xmin, double xmax, double ymin, double ymax, unsigned cid, double shearRate, double shearWidth, bool shearForce)
{	
	this->_shearRateBox[0] = xmin;
	this->_shearRateBox[1] = xmax;
	this->_shearRateBox[2] = ymin;
	this->_shearRateBox[3] = ymax;
	this->_shearRate = shearRate;
	this->_shearWidth = shearWidth;
	this->_shearComp = cid;
	this->_doApplyShearForce = shearForce;
	this->_shearRampTime = ceil((_shearRateBox[3]-_shearRateBox[2])*_shearRate*10000);
	cout << "Shear info: the final shear velocity is firstly reached " << _shearRampTime << " timesteps after the initStatistics!";
}

void PressureGradient::prepareShearRate(ParticleContainer* molCont, DomainDecompBase* domainDecomp, unsigned directedVelTime)
{	
	
	
	unsigned yuns;
	unsigned yun;
	
	// if shearWidth > 0.0 -->  the shear velocity is controlled in 2 stripes at the box margins and 2 stripes in the middle of the box
	// if shearWidth == 0.0 --> the shear velocity is controlled in n stripes with a width of 0.1 sigma; more accurate v-profile 
	if(this->_shearWidth > 0.0)
		yuns = 4;
	else
		yuns = ceil((_shearRateBox[3] - _shearRateBox[2])*10);
	
	for(yun = 0; yun < yuns; yun++){
	    this->_localShearN[yun] = 0;
	    this->_globalShearN[yun] = 0;
	    this->_localShearVelocitySum[yun] = 0.0;
	    this->_globalShearVelocitySum[yun] = 0.0;
	}
	
	if(!this->_localRank && _globalPriorShearVelocitySums[0].size() == 0){
	  for(yun = 0; yun < yuns; yun++){
	    this->_averagedShearVelocitySum[yun] = 0.0;
	    this->_averagedShearN[yun] = 0;
	  }
	}
	  
	  
	for(Molecule* thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
	{
	 if(this->_shearWidth > 0.0){	
	  if(thismol->componentid() == _shearComp && thismol->r(0) >= _shearRateBox[0] && thismol->r(0) <= _shearRateBox[1] && thismol->r(1) >= _shearRateBox[2] && thismol->r(1) <= _shearRateBox[2] + _shearWidth){
		this->_localShearN[0]++;
		this->_localShearVelocitySum[0] += thismol->v(0);
	  }else if(thismol->componentid() == _shearComp && thismol->r(0) >= _shearRateBox[0] && thismol->r(0) <= _shearRateBox[1] && thismol->r(1) >= _shearRateBox[2] + 0.5 * (_shearRateBox[3]-_shearRateBox[2]) - _shearWidth && thismol->r(1) < _shearRateBox[2] + 0.5 * (_shearRateBox[3]-_shearRateBox[2])){
		this->_localShearN[1]++;
		this->_localShearVelocitySum[1] += thismol->v(0);
	  }else if(thismol->componentid() == _shearComp && thismol->r(0) >= _shearRateBox[0] && thismol->r(0) <= _shearRateBox[1] && thismol->r(1) >= _shearRateBox[2] + 0.5 * (_shearRateBox[3]-_shearRateBox[2]) && thismol->r(1) <= _shearRateBox[2] + 0.5 * (_shearRateBox[3]-_shearRateBox[2]) + _shearWidth){
		this->_localShearN[2]++;
		this->_localShearVelocitySum[2] += thismol->v(0);	
	  }else if(thismol->componentid() == _shearComp && thismol->r(0) >= _shearRateBox[0] && thismol->r(0) <= _shearRateBox[1] && thismol->r(1) >= _shearRateBox[3] - _shearWidth && thismol->r(1) <= _shearRateBox[3]){
		this->_localShearN[3]++;
		this->_localShearVelocitySum[3] += thismol->v(0);
	  }
	 }else{
		if(thismol->componentid() == _shearComp && thismol->r(0) >= _shearRateBox[0] && thismol->r(0) <= _shearRateBox[1] && thismol->r(1) >= _shearRateBox[2] && thismol->r(1) <= _shearRateBox[3]){ 
			yun = floor((thismol->r(1) - this->_shearRateBox[2])*10);
			this->_localShearN[yun]++;
			this->_localShearVelocitySum[yun] += thismol->v(0);
		}
	 }
	}
	
	domainDecomp->collCommInit( 2*yuns );
	for(yun = 0; yun < yuns; yun++){
	    domainDecomp->collCommAppendUnsLong(this->_localShearN[yun]);
	    domainDecomp->collCommAppendDouble( this->_localShearVelocitySum[yun]);
	}
	domainDecomp->collCommAllreduceSum();
	for(yun = 0; yun < yuns; yun++){
	    this->_globalShearN[yun] = domainDecomp->collCommGetUnsLong();
	    this->_globalShearVelocitySum[yun] = domainDecomp->collCommGetDouble();
	}
	domainDecomp->collCommFinalize();
	for(yun = 0; yun < yuns; yun++){
	  if(this->_globalShearN[yun] != 0)
	    this->_directedShearVel[yun] = this->_globalShearVelocitySum[yun]/this->_globalShearN[yun];
	  else
	    this->_directedShearVel[yun] = 0.0;
	}
	
	map<unsigned int, long double>::iterator gVSit;
	if(!this->_localRank)
	{
		for(gVSit = this->_globalShearVelocitySum.begin(); gVSit != this->_globalShearVelocitySum.end(); gVSit++)
		{
				this->_globalPriorShearVelocitySums[gVSit->first].push_back(this->_globalShearVelocitySum[gVSit->first]);
				this->_globalPriorShearN[gVSit->first].push_back(this->_globalShearN[gVSit->first]);
				this->_averagedShearVelocitySum[gVSit->first] += this->_globalShearVelocitySum[gVSit->first]; 
				this->_averagedShearN[gVSit->first] += this->_globalShearN[gVSit->first];
		}
	}

	if(!this->_localRank)
	{
		for(gVSit = this->_globalShearVelocitySum.begin(); gVSit != this->_globalShearVelocitySum.end(); gVSit++)
		{
			this->_directedShearVelAverage[gVSit->first] = this->_averagedShearVelocitySum[gVSit->first]/this->_averagedShearN[gVSit->first];
			if(this->_globalPriorShearVelocitySums[gVSit->first].size() == directedVelTime){
			  this->_averagedShearVelocitySum[gVSit->first] -= this->_globalPriorShearVelocitySums[gVSit->first].front();
			  this->_averagedShearN[gVSit->first] -= this->_globalPriorShearN[gVSit->first].front();
			  this->_globalPriorShearVelocitySums[gVSit->first].pop_front();
			  this->_globalPriorShearN[gVSit->first].pop_front();
			}
		}
	}

	domainDecomp->collCommInit(yuns);
	for( map<unsigned int, long double>::iterator lNit = this->_globalShearVelocitySum.begin();
			lNit != this->_globalShearVelocitySum.end();
			lNit++ )
		domainDecomp->collCommAppendDouble(this->_directedShearVelAverage[lNit->first]);
	domainDecomp->collCommBroadcast();
	for( map<unsigned int, long double>::iterator lNit = this->_globalShearVelocitySum.begin();
			lNit != this->_globalShearVelocitySum.end();
			lNit++ ){
		this->_directedShearVelAverage[lNit->first] = domainDecomp->collCommGetDouble();
	}
	domainDecomp->collCommFinalize();
}

void PressureGradient::calculateForcesOnComponent(ParticleContainer* molCont, unsigned int cid)
{
	for(Molecule* thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
	{	  
	  if(thismol->componentid() == cid){
		this->_localN[cid]++;
		for(unsigned short int d = 0; d < 3; d++)
			this->_localForceSum[d][cid] += thismol->F(d);
	  }
	}
}

void PressureGradient::collectForcesOnComponent(DomainDecompBase* domainDecomp, unsigned int cid)
{
	domainDecomp->collCommInit( 4 * _localN.size() );
	for( map<unsigned int, unsigned int long>::iterator lNit = _localN.begin(); lNit != _localN.end(); lNit++ )
	{
		domainDecomp->collCommAppendUnsLong(lNit->second);
		for(int d = 0; d < 3; d++)
			domainDecomp->collCommAppendDouble( this->_localForceSum[d][lNit->first]);
	}
	domainDecomp->collCommAllreduceSum();
	for( map<unsigned int, unsigned long>::iterator lNit = _localN.begin(); lNit != _localN.end(); lNit++ )
	{
		_globalN[lNit->first] = domainDecomp->collCommGetUnsLong();
		for(int d = 0; d < 3; d++)
			_globalForceSum[d][lNit->first] = domainDecomp->collCommGetDouble();
	}
	
	domainDecomp->collCommFinalize();
}

void PressureGradient::resetForcesOnComponent(unsigned int cid)
{
	this->_localN[cid] = 0;
	this->_globalN[cid] = 0;
	for(unsigned short int d = 0; d < 3; d++){
		this->_localForceSum[d][cid] = 0.0;
		this->_globalForceSum[d][cid] = 0.0;
	}
}

void PressureGradient::calculateSpringForcesOnComponent(ParticleContainer* molCont, unsigned int cid)
{
	for(Molecule* thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
	{
	  if(thismol->componentid() == cid){
		    this->_localSpringN[cid]++;
		    for(unsigned short int d = 0; d < 3; d++)
			    this->_localSpringForceSum[d][cid] += thismol->F_Spring(d);
	  }
	}
}

void PressureGradient::collectSpringForcesOnComponent(DomainDecompBase* domainDecomp, unsigned int cid)
{
	    domainDecomp->collCommInit( 4 * _localSpringN.size() );
	    for( map<unsigned int, unsigned int long>::iterator lNit = _localSpringN.begin(); lNit != _localSpringN.end(); lNit++ )
	    {
		domainDecomp->collCommAppendUnsLong(lNit->second);
		for(int d = 0; d < 3; d++)
			domainDecomp->collCommAppendDouble( this->_localSpringForceSum[d][lNit->first]);
	    }
	    domainDecomp->collCommAllreduceSum();
	    for( map<unsigned int, unsigned long>::iterator lNit = _localSpringN.begin(); lNit != _localSpringN.end(); lNit++ )
	    {
		_globalSpringN[lNit->first] = domainDecomp->collCommGetUnsLong();
		for(int d = 0; d < 3; d++)
			_globalSpringForceSum[d][lNit->first] = domainDecomp->collCommGetDouble();
	    }
	domainDecomp->collCommFinalize();
}

void PressureGradient::resetSpringForcesOnComponent(unsigned int cid)
{
	this->_localSpringN[cid] = 0;
	this->_globalSpringN[cid] = 0;
	for(unsigned short int d = 0; d < 3; d++){
		this->_localSpringForceSum[d][cid] = 0.0;
		this->_globalSpringForceSum[d][cid] = 0.0;
	}
}

void PressureGradient::specifyTauPrime(double tauPrime, double dt)
{
	this->_universalTauPrime = tauPrime;
	if(this->_localRank != 0) return;
	if(this->_universalConstantAccelerationTimesteps == 0)
	{
		cout << "SEVERE ERROR: unknown UCAT!\n";
		exit(78);
	}
	unsigned int vql = (unsigned int)ceil(tauPrime / (dt*this->_universalConstantAccelerationTimesteps));
	map<unsigned int, unsigned int>::iterator vqlit;
	for(vqlit = _globalVelocityQueuelength.begin(); vqlit != _globalVelocityQueuelength.end(); vqlit++)
	{
		vqlit->second = vql;
		cout << "coset " << vqlit->first << " will receive "
		     << vqlit->second << " velocity queue entries.\n";
	}
}

/*
 * quadratische Variante
 */
void PressureGradient::adjustTau(double dt) {
	if(this->_universalConstantTau) return;
	map<unsigned int, double>::iterator tauit;
	double increment;
	for(tauit = _universalTau.begin(); tauit != _universalTau.end(); tauit++)
	{
		increment = dt * this->_universalZetaFlow * sqrt(tauit->second);
		tauit->second += increment;
	}
}
