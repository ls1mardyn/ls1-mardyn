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
	if(!rank)
		for(unsigned short int d=0; d < 3; d++)
			this->_globalVelocitySum[d] = map<unsigned int, long double>();
	this->_universalConstantTau = true;
	this->_universalZetaFlow = 0.0;
	this->_universalTauPrime = 0.0;
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

void PressureGradient::specifySpringInfluence(unsigned long minSpringID, unsigned long maxSpringID, double averageYPos, double springConst)
{
	this->_minSpringID = minSpringID;
	this->_maxSpringID = maxSpringID;
	this->_averageYPos = averageYPos;
	this->_springConst = springConst;
	
	if(!this->_localRank)
	    cout << "Spring Force Effect is prepared for MolIDs " << this->_minSpringID << " to " << this->_maxSpringID << " with minimum y-Position " << this->_averageYPos << endl;
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

void PressureGradient::prepare_getMissingVelocity(DomainDecompBase* domainDecomp, ParticleContainer* molCont, unsigned int cid, unsigned numberOfComp)
{
	for (unsigned compID = 0; compID < numberOfComp; compID++){
	  this->_localN[compID] = 0;
	  for(unsigned short int d = 0; d < 3; d++)
		this->_localVelocitySum[d][compID] = 0.0;
	}
	for(Molecule* thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
	{
	  if(thismol->componentid() == cid){
		this->_localN[cid]++;
		for(unsigned short int d = 0; d < 3; d++)
			this->_localVelocitySum[d][cid] += thismol->v(d);
	  }
	}

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
}

double PressureGradient::getMissingVelocity(unsigned int cid, unsigned short int d)
{
	double v_directed = this->_globalVelocitySum[d][cid] / this->_globalN[cid];
	double v_missing = 0.0;
	
	if (this->_globalTargetVelocity[d][cid] != 0.0)
	  v_missing = this->_globalTargetVelocity[d][cid] - v_directed;
	
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

void PressureGradient::setupShearRate(double xmin, double xmax, double ymin, double ymax, unsigned cid, double shearRate)
{	
	this->_shearRateBox[0] = xmin;
	this->_shearRateBox[1] = xmax;
	this->_shearRateBox[2] = ymin;
	this->_shearRateBox[3] = ymax;
	this->_shearRate = shearRate;
	this->_shearComp = cid;
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

