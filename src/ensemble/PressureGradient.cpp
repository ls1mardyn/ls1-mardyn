#include "PressureGradient.h"

#include <cmath>
#include <iostream>

#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "utils/mardyn_assert.h"


PressureGradient::PressureGradient(int rank) {
	this->_localRank = rank;
	this->_universalConstantAccelerationTimesteps = 0;
	if(!rank)
		for(unsigned short int d=0; d < 3; d++)
			this->_globalVelocitySum[d] = std::map<unsigned int, long double>();
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
		std::cout << "coset " << cosetid << " will receive "
			<< _globalVelocityQueuelength[cosetid] << " velocity queue entries." << std::endl;
	}
}

/*
 * diese Version beschleunigt in alle Raumrichtungen
 */
void PressureGradient::determineAdditionalAcceleration
(
 DomainDecompBase* domainDecomp,
 ParticleContainer* molCont, double dtConstantAcc )
{
	for( auto uAAit = _universalAdditionalAcceleration[0].begin();
			uAAit != _universalAdditionalAcceleration[0].end();
			uAAit++ )
	{
		this->_localN[uAAit->first] = 0;
		for(unsigned short int d = 0; d < 3; d++)
			this->_localVelocitySum[d][uAAit->first] = 0.0;
	}

	// TODO: consider parallelization, but when we have input

	for(auto thismol = molCont->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); thismol.isValid(); ++thismol)
	{
		unsigned int cid = thismol->componentid();
		auto uCSIDit = this->_universalComponentSetID.find(cid);
		if(uCSIDit == _universalComponentSetID.end()) continue;
		unsigned cosetid = uCSIDit->second;
		this->_localN[cosetid]++;
		for(unsigned short int d = 0; d < 3; d++)
			this->_localVelocitySum[d][cosetid] += thismol->v(d);
	}

	// domainDecomp->collectCosetVelocity(&_localN, _localVelocitySum, &_globalN, _globalVelocitySum);
	//
	domainDecomp->collCommInit( 4 * _localN.size() );
	for( auto lNit = _localN.begin(); lNit != _localN.end(); lNit++ )
	{
		domainDecomp->collCommAppendUnsLong(lNit->second);
		for(int d = 0; d < 3; d++)
			domainDecomp->collCommAppendDouble( this->_localVelocitySum[d][lNit->first]);
	}
	domainDecomp->collCommAllreduceSum();
	for(auto lNit = _localN.begin(); lNit != _localN.end(); lNit++ )
	{
		_globalN[lNit->first] = domainDecomp->collCommGetUnsLong();
		for(int d = 0; d < 3; d++)
			_globalVelocitySum[d][lNit->first] = domainDecomp->collCommGetDouble();
	}
	domainDecomp->collCommFinalize();

	std::map<unsigned int, long double>::iterator gVSit;
	if(!this->_localRank)
	{
		for(gVSit = _globalVelocitySum[0].begin(); gVSit != _globalVelocitySum[0].end(); gVSit++)
		{
#ifndef NDEBUG
			Log::global_log->debug() << "required entries in velocity queue: " << _globalVelocityQueuelength[gVSit->first] << std::endl;
			Log::global_log->debug() << "entries in velocity queue: " << _globalPriorVelocitySums[0][gVSit->first].size() << std::endl;
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
					+= static_cast<long double>(dtConstantAcc) * invgtau2 *
					( this->_globalTargetVelocity[d][gVSit->first]
					  - 2.0*this->_globalVelocitySum[d][gVSit->first]*invgN
					  + previousVelocity[d] );
			}
#ifndef NDEBUG
			Log::global_log->debug() << "accelerator no. " << gVSit->first
				<< "previous vz: " << previousVelocity[2]
				<< "current vz: " << _globalVelocitySum[2][gVSit->first]*invgN << std::endl;
#endif
		}
	}

	domainDecomp->collCommInit(7*this->_localN.size());
	for( auto lNit = _localN.begin();
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
	for( auto lNit = _localN.begin();
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

double PressureGradient::getMissingVelocity(unsigned int cosetid, unsigned short int d)
{
	double vd = this->_globalVelocitySum[d][cosetid] / this->_globalN[cosetid];
	return this->_globalTargetVelocity[d][cosetid] - vd;
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

void PressureGradient::specifyTauPrime(double tauPrime, double dt)
{
	this->_universalTauPrime = tauPrime;
	if(this->_localRank != 0) return;
	if(this->_universalConstantAccelerationTimesteps == 0)
	{
		std::ostringstream error_message;
		error_message << "SEVERE ERROR: unknown UCAT!\n";
		MARDYN_EXIT(error_message);
	}
	unsigned int vql = (unsigned int)ceil(tauPrime / (dt*this->_universalConstantAccelerationTimesteps));
	std::map<unsigned int, unsigned int>::iterator vqlit;
	for(vqlit = _globalVelocityQueuelength.begin(); vqlit != _globalVelocityQueuelength.end(); vqlit++)
	{
		vqlit->second = vql;
		std::cout << "coset " << vqlit->first << " will receive "
		     << vqlit->second << " velocity queue entries.\n";
	}
}

/*
 * quadratische Variante
 */
void PressureGradient::adjustTau(double dt) {
	if(this->_universalConstantTau) return;
	std::map<unsigned int, double>::iterator tauit;
	double increment;
	for(tauit = _universalTau.begin(); tauit != _universalTau.end(); tauit++)
	{
		increment = dt * this->_universalZetaFlow * sqrt(tauit->second);
		tauit->second += increment;
	}
}

