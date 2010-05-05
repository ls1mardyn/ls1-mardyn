/***************************************************************************
 *   Copyright (C) 2010 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "Domain.h"

#include "datastructures/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "ensemble/GrandCanonical.h"

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

#include "Logger.h"
using Log::global_log;

using namespace std;

#define VERSION 20100321
#define RDF_MINIMAL_OUTPUT_STEPS 1023
#define MIN_BETA 0.9
#define KINLIMIT_PER_T 10.0

/*
 * cutoff correction
 */
double TICCu(int n,double rc,double sigma2);
double TICSu(int n,double rc,double sigma2,double tau);
double TISSu(int n,double rc,double sigma2,double tau1,double tau2);
double TICCv(int n,double rc,double sigma2);
double TICSv(int n,double rc,double sigma2,double tau);
double TISSv(int n,double rc,double sigma2,double tau1,double tau2);

double TICCu(int n,double rc,double sigma2)
{
	return -pow(rc,2*n+3) / (pow(sigma2,n)*(2*n+3));
}
double TICSu(int n,double rc,double sigma2,double tau)
{
	return -( pow(rc+tau,2*n+3) - pow(rc-tau,2*n+3) ) * rc / ( 4*pow(sigma2,n)*tau*(n+1)*(2*n+3) ) 
		+  ( pow(rc+tau,2*n+4) - pow(rc-tau,2*n+4) ) / ( 4*pow(sigma2,n)*tau*(n+1)*(2*n+3)*(2*n+4) );
}
double TISSu(int n,double rc,double sigma2,double tau1,double tau2)
{
	double tauMinus,tauPlus;
	tauPlus = tau1+tau2;
	tauMinus = tau1-tau2;
	return -(   pow(rc+tauPlus,2*n+4) - pow(rc+tauMinus,2*n+4) - pow(rc-tauMinus,2*n+4) + pow(rc-tauPlus,2*n+4) ) * rc / ( 8*pow(sigma2,n)*tau1*tau2*(n+1)*(2*n+3)*(2*n+4) ) +  (   pow(rc+tauPlus,2*n+5) - pow(rc+tauMinus,2*n+5) - pow(rc-tauMinus,2*n+5) + pow(rc-tauPlus,2*n+5) ) / ( 8*pow(sigma2,n)*tau1*tau2*(n+1)*(2*n+3)*(2*n+4)*(2*n+5) );
}
double TICCv(int n,double rc,double sigma2)
{
	return 2*n * TICCu(n,rc,sigma2);
}
double TICSv(int n,double rc,double sigma2,double tau)
{
	return -( pow(rc+tau,2*n+2) - pow(rc-tau,2*n+2) ) * rc*rc / ( 4*pow(sigma2,n)*tau*(n+1) ) - 3*TICSu(n,rc,sigma2,tau);
}

double TISSv(int n,double rc,double sigma2,double tau1,double tau2){
	double tauMinus,tauPlus;
	tauPlus = tau1+tau2;
	tauMinus = tau1-tau2;
	return -(   pow(rc+tauPlus,2*n+3) - pow(rc+tauMinus,2*n+3) - pow(rc-tauMinus,2*n+3) + pow(rc-tauPlus,2*n+3) ) * rc*rc / ( 8*pow(sigma2,n)*tau1*tau2*(n+1)*(2*n+3) ) - 3*TISSu(n,rc,sigma2,tau1,tau2);
}

/*
 * END cutoff correction
 */

Domain::Domain(int rank){
	_localRank = rank;
	_localUpot = 0;
	_localVirial = 0;   
	_globalUpot = 0;
	_globalVirial = 0; 
	_globalRho = 0;

	this->_componentToThermostatIdMap = map<int, int>();
	this->_localThermostatN = map<int, unsigned long>();
	this->_localThermostatN[-1] = 0;
	this->_localThermostatN[0] = 0;
	this->_universalThermostatN = map<int, unsigned long>();
	this->_universalThermostatN[-1] = 0;
	this->_universalThermostatN[0] = 0;
	this->_localRotationalDOF = map<int, unsigned long>();
	this->_localRotationalDOF[-1] = 0;
	this->_localRotationalDOF[0] = 0;
	this->_universalRotationalDOF = map<int, unsigned long>();
	this->_universalRotationalDOF[-1] = 0;
	this->_universalRotationalDOF[0] = 0;
	this->_globalLength[0] = 0;
	this->_globalLength[1] = 0;
	this->_globalLength[2] = 0;
	this->_universalBTrans = map<int, double>();
	this->_universalBTrans[0] = 1.0;
	this->_universalBRot = map<int, double>();
	this->_universalBRot[0] = 1.0;
	this->_universalTargetTemperature = map<int, double>();
	this->_universalTargetTemperature[0] = 1.0;
	this->_globalTemperatureMap = map<int, double>();
	this->_globalTemperatureMap[0] = 1.0;
	this->_local2KETrans[0] = 0.0;
	this->_local2KERot[0] = 0.0; 

	_currentTime = 0.0;

	this->_doCollectRDF = false;
	this->_universalRDFTimesteps = -1;
	this->_universalNVE = false;
	this->_globalUSteps = 0;
	this->_globalSigmaU = 0.0;
	this->_globalSigmaUU = 0.0;
#ifdef COMPLEX_POTENTIAL_SET
	this->_universalConstantAccelerationTimesteps = 30;
	if(!rank)
		for(int d = 0; d < 3; d++)
			this->_globalVelocitySum[d] = map<unsigned, long double>();
#endif
	this->_componentwiseThermostat = false;
#ifdef COMPLEX_POTENTIAL_SET
	this->_universalUndirectedThermostat = map<int, bool>();
	for(int d = 0; d < 3; d++)
	{
		this->_universalThermostatDirectedVelocity[d] = map<int, double>();
		this->_localThermostatDirectedVelocity[d] = map<int, double>();
	}
#endif
	this->_universalSelectiveThermostatCounter = 0;
	this->_universalSelectiveThermostatWarning = 0;
	this->_universalSelectiveThermostatError = 0;
}

void Domain::setLocalUpot(double Upot) {_localUpot = Upot;}

double Domain::getLocalUpot() const {return _localUpot; }

void Domain::setLocalVirial(double Virial) {_localVirial = Virial;}

double Domain::getLocalVirial() const {return _localVirial; }

/* methods accessing thermostat info */
double Domain::getGlobalBetaTrans() { return _universalBTrans[0]; }
double Domain::getGlobalBetaTrans(int thermostat) { return _universalBTrans[thermostat]; }
double Domain::getGlobalBetaRot() { return _universalBRot[0]; }
double Domain::getGlobalBetaRot(int thermostat) { return _universalBRot[thermostat]; }

void Domain::setLocalSummv2(double summv2, int thermostat)
{
#ifndef NDEBUG
	global_log->debug() << "* local thermostat " << thermostat << ":  mvv = " << summv2 << endl;
#endif
	this->_local2KETrans[thermostat] = summv2;
}

void Domain::setLocalSummv2(double summv2) { setLocalSummv2(summv2, 0); }

void Domain::setLocalSumIw2(double sumIw2, int thermostat)
{
	_local2KERot[thermostat] = sumIw2;
} 

void Domain::setLocalSumIw2(double sumIw2) { setLocalSumIw2(sumIw2, 0); }

double Domain::getGlobalPressure()
{
	double globalTemperature = _globalTemperatureMap[0];
	return globalTemperature * _globalRho + _globalRho * getAverageGlobalVirial()/3.;
}

double Domain::getAverageGlobalVirial() const { return _globalVirial/_globalNumMolecules; }

double Domain::getAverageGlobalUpot() const { return _globalUpot/_globalNumMolecules; }



double Domain::getCurrentTime(){ return _currentTime;}

void Domain::setCurrentTime(double curtime){ _currentTime = curtime;}

void Domain::advanceTime(double timestep){ _currentTime += timestep;}

vector<Component>& Domain::getComponents(){
	return _components; 
}

void Domain::addComponent(Component component){
	_components.push_back(component);
}

Comp2Param& Domain::getComp2Params(){
	return _comp2params; 
}

void Domain::calculateGlobalValues(
		DomainDecompBase* domainDecomp,
		ParticleContainer* particleContainer,
		bool collectThermostatVelocities,
		double Tfactor
		) {
	double Upot = _localUpot;
	double Virial = _localVirial;

	// To calculate Upot, Ukin and Pressure, intermediate values from all      
	// processes are needed. Here the         
	// intermediate values of all processes are summed up so that the root    
	// process can calculate the final values. to be able to calculate all     
	// values at this point, the calculation of the intermediate value sum_v2  
	// had to be moved from Thermostat to upd_postF and the final calculations  
	// of m_Ukin, m_Upot and Pressure had to be moved from Thermostat / upd_F  
	// to this point           
	domainDecomp->collCommInit(2);
	domainDecomp->collCommAppendDouble(Upot);
	domainDecomp->collCommAppendDouble(Virial);
	domainDecomp->collCommAllreduceSum();
	Upot = domainDecomp->collCommGetDouble();
	Virial = domainDecomp->collCommGetDouble();
	domainDecomp->collCommFinalize();

	/* FIXME: why should process 0 do this alone? 
	 * we should keep symmetry of all proccesses! */
	// Process 0 has to add the dipole correction:
	// m_UpotCorr and m_VirialCorr already contain constant (internal) dipole correction
	_globalUpot = Upot + _UpotCorr;
	_globalVirial = Virial + _VirialCorr;

	/*
	 * thermostat ID 0 represents the entire system
	 */
	map<int, unsigned long>::iterator thermit;
	if( _componentwiseThermostat )
	{
#ifndef NDEBUG
		global_log->debug() << "* applying a componentwise thermostat" << endl;
#endif
		this->_localThermostatN[0] = 0;
		this->_localRotationalDOF[0] = 0;
		this->_local2KETrans[0] = 0;
		this->_local2KERot[0] = 0;
		for(thermit = _localThermostatN.begin(); thermit != _localThermostatN.end(); thermit++)
		{
			if(thermit->first == 0) continue;
			this->_localThermostatN[0] += thermit->second;
			this->_localRotationalDOF[0] += this->_localRotationalDOF[thermit->first];
			this->_local2KETrans[0] += this->_local2KETrans[thermit->first];
			this->_local2KERot[0] += this->_local2KERot[thermit->first];
		}
	}

	for(thermit = _universalThermostatN.begin(); thermit != _universalThermostatN.end(); thermit++)
	{
		// number of molecules on the local process. After the reduce operation
		// num_molecules will contain the global number of molecules
		unsigned long numMolecules = _localThermostatN[thermit->first];
		double summv2 = _local2KETrans[thermit->first];
		assert(summv2 >= 0.0);
		unsigned long rotDOF = _localRotationalDOF[thermit->first];
		double sumIw2 = (rotDOF > 0)? _local2KERot[thermit->first]: 0.0;

		domainDecomp->collCommInit(4);
		domainDecomp->collCommAppendDouble(summv2);
		domainDecomp->collCommAppendDouble(sumIw2);
		domainDecomp->collCommAppendUnsLong(numMolecules);
		domainDecomp->collCommAppendUnsLong(rotDOF);
		domainDecomp->collCommAllreduceSum();
		summv2 = domainDecomp->collCommGetDouble();
		sumIw2 = domainDecomp->collCommGetDouble();
		numMolecules = domainDecomp->collCommGetUnsLong();
		rotDOF = domainDecomp->collCommGetUnsLong();
		domainDecomp->collCommFinalize();

		this->_universalThermostatN[thermit->first] = numMolecules;
		this->_universalRotationalDOF[thermit->first] = rotDOF;

		/* calculate the temperature of the entire system */
		if(numMolecules > 0)
			_globalTemperatureMap[thermit->first] =
				(summv2 + sumIw2) / (double)(3*numMolecules + rotDOF);
		else
			_globalTemperatureMap[thermit->first] = _universalTargetTemperature[thermit->first];

		double Ti = Tfactor * _universalTargetTemperature[thermit->first];
		if((Ti > 0.0) && !_universalNVE)
		{
			_universalBTrans[thermit->first] = pow(3.0*numMolecules*Ti / summv2, 0.4);
			if( sumIw2 == 0.0 ) 
				_universalBRot[thermit->first] = 1.0;
			else 
				_universalBRot[thermit->first] = pow(rotDOF*Ti / sumIw2, 0.4);
		}
		else
		{
			this->_universalBTrans[thermit->first] = 1.0;
			this->_universalBRot[thermit->first] = 1.0;
		}

		// heuristic handling of the unfortunate special case of an explosion in the system
		if( ( (_universalBTrans[thermit->first] < MIN_BETA) || (_universalBRot[thermit->first] < MIN_BETA) )
				&& (0 >= _universalSelectiveThermostatError) )
		{
			global_log->warning() << "Explosion warning (time t=" << _currentTime << ")." << endl;
			global_log->debug() << "Selective thermostat will be applied to set " << thermit->first
				<< " (beta_trans = " << this->_universalBTrans[thermit->first]
				<< ", beta_rot = " << this->_universalBRot[thermit->first] << "!)" << endl;
			int rot_dof;
			double Utrans, Urot;
			double limit_energy =  KINLIMIT_PER_T * Ti;
			double limit_rot_energy;
			double vcorr, Dcorr;
			Molecule* tM;
			for( tM = particleContainer->begin();
					tM != particleContainer->end();
					tM = particleContainer->next() )
			{
				Utrans = tM->Utrans();
				if(Utrans > limit_energy)
				{
					vcorr = sqrt(limit_energy / Utrans);
					global_log->debug() << ": v(m" << tM->id() << ") *= " << vcorr << endl;
					tM->scale_v(vcorr);
					tM->scale_F(vcorr);
				}

				rot_dof = _components[tM->componentid()].rot_dof();
				if(rot_dof > 0)
				{
					limit_rot_energy = 3.0*rot_dof * Ti;
					Urot = tM->Urot();
					if(Urot > limit_rot_energy)
					{
						Dcorr = sqrt(limit_rot_energy / Urot);
						global_log->debug() << "D(m" << tM->id() << ") *= " << Dcorr << endl;
						tM->scale_D(Dcorr);
						tM->scale_F(Dcorr);
					}
				}
			}

			/* FIXME: Unnamed constant 3960... */
			if(3960 >= _universalSelectiveThermostatCounter)
			{
				if( _universalSelectiveThermostatWarning > 0 )
					_universalSelectiveThermostatError = _universalSelectiveThermostatWarning;
				if( _universalSelectiveThermostatCounter > 0 )
					_universalSelectiveThermostatWarning = _universalSelectiveThermostatCounter;
				_universalSelectiveThermostatCounter = 4000;
			}
			_universalBTrans[thermit->first] = 1.0;
			_universalBRot[thermit->first] = pow(this->_universalBRot[thermit->first], 0.0091);
		}
#ifdef NDEBUG
		if( (_universalSelectiveThermostatCounter > 0) &&
				((_universalSelectiveThermostatCounter % 20) == 10) )
#endif
			/* FIXME: why difference counters? */
			global_log->debug() << "counter " << _universalSelectiveThermostatCounter
				<< ",\t warning " << _universalSelectiveThermostatWarning
				<< ",\t error " << _universalSelectiveThermostatError << endl;

		if(collectThermostatVelocities && _universalUndirectedThermostat[thermit->first])
		{
			double sigv[3];
			for(int d=0; d < 3; d++)
				sigv[d] = _localThermostatDirectedVelocity[d][thermit->first];

			domainDecomp->collCommInit(3);
			for(int d=0; d < 3; d++) domainDecomp->collCommAppendDouble(sigv[d]);
			domainDecomp->collCommAllreduceSum();
			for(int d=0; d < 3; d++) sigv[d] = domainDecomp->collCommGetDouble();
			domainDecomp->collCommFinalize();

			for(int d=0; d < 3; d++)
			{
				_localThermostatDirectedVelocity[d][thermit->first] = 0.0;
				if(numMolecules > 0)
					_universalThermostatDirectedVelocity[d][thermit->first] = sigv[d] / numMolecules;
				else 
					_universalThermostatDirectedVelocity[d][thermit->first] = 0.0;
			}

#ifndef NDEBUG
			global_log->debug() << "* thermostat " << thermit->first
				<< " directed velocity: ("
				<< _universalThermostatDirectedVelocity[0][thermit->first]
				<< " / " << _universalThermostatDirectedVelocity[1][thermit->first]
				<< " / " << _universalThermostatDirectedVelocity[2][thermit->first] 
				<< ")" << endl;
#endif
		}

#ifndef NDEBUG
		global_log->debug() << "* Th" << thermit->first << " N=" << numMolecules
			<< " DOF=" << rotDOF + 3.0*numMolecules
			<< " Tcur=" << _globalTemperatureMap[thermit->first]
			<< " Ttar=" << _universalTargetTemperature[thermit->first]
			<< " Tfactor=" << Tfactor
			<< " bt=" << _universalBTrans[thermit->first]
			<< " br=" << _universalBRot[thermit->first] << "\n";
#endif
	}

	if(this->_universalSelectiveThermostatCounter > 0)
		this->_universalSelectiveThermostatCounter--;
	if(this->_universalSelectiveThermostatWarning > 0)
		this->_universalSelectiveThermostatWarning--;
	if(this->_universalSelectiveThermostatError > 0)
		this->_universalSelectiveThermostatError--;
}

void Domain::calculateThermostatDirectedVelocity(ParticleContainer* partCont)
{
	Molecule* tM;
	if(this->_componentwiseThermostat)
	{
		for( map<int, bool>::iterator thit = _universalUndirectedThermostat.begin();
				thit != _universalUndirectedThermostat.end();
				thit ++ )
		{
			if(thit->second)
				for(int d=0; d < 3; d++) _localThermostatDirectedVelocity[d][thit->first] = 0.0;
		}
		for(tM = partCont->begin(); tM != partCont->end(); tM = partCont->next() )
		{
			int cid = tM->componentid();
			int thermostat = this->_componentToThermostatIdMap[cid];
			if(this->_universalUndirectedThermostat[thermostat])
			{
				for(int d=0; d < 3; d++)
					_localThermostatDirectedVelocity[d][thermostat] += tM->v(d);
			}
		}
	}
	else if(this->_universalUndirectedThermostat[0])
	{
		for(int d=0; d < 3; d++) _localThermostatDirectedVelocity[d][0] = 0.0;
		for(tM = partCont->begin(); tM != partCont->end(); tM = partCont->next() )
		{
			for(int d=0; d < 3; d++)
				_localThermostatDirectedVelocity[d][0] += tM->v(d);
		}
	}
}

void Domain::calculateVelocitySums(ParticleContainer* partCont)
{
	Molecule* tM;
	if(this->_componentwiseThermostat)
	{
		for(tM = partCont->begin(); tM != partCont->end(); tM = partCont->next() )
		{
			int cid = tM->componentid();
			int thermostat = this->_componentToThermostatIdMap[cid];
			this->_localThermostatN[thermostat]++;
			this->_localRotationalDOF[thermostat] += _components[cid].rot_dof();
			if(this->_universalUndirectedThermostat[thermostat])
			{
				tM->calculate_mv2_Iw2( this->_local2KETrans[thermostat],
						this->_local2KERot[thermostat],
						this->_universalThermostatDirectedVelocity[0][thermostat],
						this->_universalThermostatDirectedVelocity[1][thermostat],
						this->_universalThermostatDirectedVelocity[2][thermostat]  );
			}
			else
			{
				tM->calculate_mv2_Iw2(_local2KETrans[thermostat], _local2KERot[thermostat]);
			}
		}
	}
	else
	{
		for(tM = partCont->begin(); tM != partCont->end(); tM = partCont->next() )
		{
			this->_localThermostatN[0]++;
			this->_localRotationalDOF[0] += _components[ tM->componentid() ].rot_dof();
			if(this->_universalUndirectedThermostat[0])
			{
				tM->calculate_mv2_Iw2( this->_local2KETrans[0],
						this->_local2KERot[0],
						this->_universalThermostatDirectedVelocity[0][0],
						this->_universalThermostatDirectedVelocity[1][0],
						this->_universalThermostatDirectedVelocity[2][0]  );
			}
			else
			{
				tM->calculate_mv2_Iw2(_local2KETrans[0], _local2KERot[0]);
			}
		}
#ifndef NDEBUG
		global_log->debug() << "      * N = " << this->_localThermostatN[0]
			<< "   rotDOF = " << this->_localRotationalDOF[0] << "   mvv = "
			<< _local2KETrans[0] << " Iww = " << _local2KERot[0] << endl;
#endif
	}
}

void Domain::writeCheckpoint( string filename, 
		ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp )
{
	domainDecomp->assertDisjunctivity(particleContainer);
	if(!this->_localRank)
	{
		ofstream checkpointfilestream(filename.c_str());
		checkpointfilestream << "mardyn trunk " << VERSION;
		checkpointfilestream << "\n";
		checkpointfilestream << " currentTime\t"  << this->_currentTime << "\n";
		checkpointfilestream << " Length\t" << setprecision(9) << _globalLength[0] << " " << _globalLength[1] << " " << _globalLength[2] << "\n";
		if(this->_componentwiseThermostat)
		{
			for( map<int, int>::iterator thermit = this->_componentToThermostatIdMap.begin();
					thermit != this->_componentToThermostatIdMap.end();
					thermit++ )
			{
				if(0 >= thermit->second) continue;
				checkpointfilestream << " CT\t" << 1+thermit->first
					<< "\t" << thermit->second << "\n";
			}
			for( map<int, double>::iterator Tit = this->_universalTargetTemperature.begin();
					Tit != this->_universalTargetTemperature.end();
					Tit++ )
			{
				if((0 >= Tit->first) || (0 >= Tit->second)) continue;
				checkpointfilestream << " ThT " << Tit->first << "\t" << Tit->second << "\n";
			}
		}
		else
		{
			checkpointfilestream << " Temperature\t" << _universalTargetTemperature[0] << endl;
		}
#ifndef NDEBUG
		checkpointfilestream << "# rho\t" << this->_globalRho << "\n";
		checkpointfilestream << "# rc\t" << particleContainer->getCutoff() << "\n";
		checkpointfilestream << "# rcT\t" << particleContainer->getTersoffCutoff() << "\n";
#endif
		if(this->_globalUSteps > 1)
		{
			checkpointfilestream << setprecision(13);
			checkpointfilestream << " I\t" << this->_globalUSteps << " "
				<< this->_globalSigmaU << " " << this->_globalSigmaUU << "\n";
			checkpointfilestream << setprecision(8);
		}
		checkpointfilestream << " NumberOfComponents\t" << _components.size() << endl;
		for(vector<Component>::const_iterator pos=_components.begin();pos!=_components.end();++pos){
			pos->write(checkpointfilestream);
		}
		unsigned int numperline=_components.size();
		unsigned int iout=0;
		for(vector<double>::const_iterator pos=_mixcoeff.begin();pos!=_mixcoeff.end();++pos){
			checkpointfilestream << *pos;
			iout++;
			// 2 parameters (xi and eta)
			if(iout/2>=numperline) {
				checkpointfilestream << endl;
				iout=0;
				--numperline;
			}
			else if(!(iout%2)) {
				checkpointfilestream << "\t";
			}
			else {
				checkpointfilestream << " ";
			}
		}
		checkpointfilestream << _epsilonRF << endl;
		for( map<unsigned, unsigned>::const_iterator uCSIDit = this->_universalComponentSetID.begin();
				uCSIDit != this->_universalComponentSetID.end();
				uCSIDit++ )
		{ 
			if(uCSIDit->first > 100) continue;
			checkpointfilestream << " S\t" << 1+uCSIDit->first << "\t" << uCSIDit->second << "\n";
		}
		for( map<unsigned, double>::const_iterator gTit = _universalTau.begin();
				gTit != this->_universalTau.end();
				gTit++ )
		{
			unsigned cosetid = gTit->first;
			checkpointfilestream << " A\t" << cosetid << "\t"
				<< this->_globalTargetVelocity[0][cosetid] << " "
				<< this->_globalTargetVelocity[1][cosetid] << " "
				<< this->_globalTargetVelocity[2][cosetid] << "\t"
				<< gTit->second << "\t"
				<< this->_universalAdditionalAcceleration[0][cosetid] << " "
				<< this->_universalAdditionalAcceleration[1][cosetid] << " "
				<< this->_universalAdditionalAcceleration[2][cosetid] << "\n";
		}
		for( map<int, bool>::iterator uutit = this->_universalUndirectedThermostat.begin();
				uutit != this->_universalUndirectedThermostat.end();
				uutit++ )
		{
			if(0 > uutit->first) continue;
			if(uutit->second) checkpointfilestream << " U\t" << uutit->first << "\n";
		}
		checkpointfilestream << " NumberOfMolecules\t" << _globalNumMolecules << endl;

		checkpointfilestream << " MoleculeFormat\t" << "ICRVQD" << endl;
		checkpointfilestream.close();
	}

	domainDecomp->writeMoleculesToFile(filename, particleContainer); 
}

void Domain::initParameterStreams(double cutoffRadius, double cutoffRadiusLJ){
	_comp2params.initialize(_components, _mixcoeff, _epsilonRF, cutoffRadius, cutoffRadiusLJ); 
}

void Domain::initFarFieldCorr(double cutoffRadius, double cutoffRadiusLJ) {
	double UpotCorrLJ=0.;
	double VirialCorrLJ=0.;
	double MySelbstTerm=0.;
	unsigned int numcomp=_components.size();
	unsigned long nummolecules=0;
	for(unsigned int i=0;i<numcomp;++i) {
		Component& ci=_components[i];
		nummolecules+=ci.numMolecules();
		unsigned int numljcentersi=ci.numLJcenters();
		unsigned int numchargesi = ci.numCharges();
		unsigned int numdipolesi=ci.numDipoles();
		unsigned int numtersoffi = ci.numTersoff();

		// effective dipoles computed from point charge distributions
		double chargeBalance[3];
		for(unsigned d = 0; d < 3; d++) chargeBalance[d] = 0;
		for(unsigned int si = 0; si < numchargesi; si++)
		{
			double tq = ci.charge(si).q();
			for(unsigned d = 0; d < 3; d++) chargeBalance[d] += tq * ci.charge(si).r()[d];
		}
		// point dipoles
		for(unsigned int si=0;si<numdipolesi;++si)
		{
			double tmy = ci.dipole(si).absMy();
			double evect = 0;
			for(unsigned d = 0; d < 3; d++) evect += ci.dipole(si).e()[d] * ci.dipole(si).e()[d];
			double norm = 1.0 / sqrt(evect);
			for(unsigned d = 0; d < 3; d++) chargeBalance[d] += tmy * ci.dipole(si).e()[d] * norm;
		}
		double my2 = 0.0;
		for(unsigned d = 0; d < 3; d++) my2 += chargeBalance[d] * chargeBalance[d];
		MySelbstTerm += my2 * ci.numMolecules();

		for(unsigned int j=0;j<numcomp;++j) {
			Component& cj=_components[j];
			unsigned numtersoffj = cj.numTersoff();
			// no LJ interaction between Tersoff components
			if(numtersoffi && numtersoffj) continue;
			unsigned int numljcentersj=cj.numLJcenters();
			ParaStrm& params=_comp2params(i,j);
			params.reset_read();
			// LJ centers
			for(unsigned int si=0;si<numljcentersi;++si) {
				double xi=ci.ljcenter(si).rx();
				double yi=ci.ljcenter(si).ry();
				double zi=ci.ljcenter(si).rz();
				double tau1=sqrt(xi*xi+yi*yi+zi*zi);
				for(unsigned int sj=0;sj<numljcentersj;++sj) {
					double xj=cj.ljcenter(sj).rx();
					double yj=cj.ljcenter(sj).ry();
					double zj=cj.ljcenter(sj).rz();
					double tau2=sqrt(xj*xj+yj*yj+zj*zj);
					if(tau1+tau2>=cutoffRadiusLJ){
						global_log->error() << "Error calculating cutoff corrections, rc too small" << endl;
						exit(1);
					}
					double eps24;
					params >> eps24;
					double sig2;
					params >> sig2;
					double uLJshift6;
					params >> uLJshift6;  // 0 unless TRUNCATED_SHIFTED

					if(uLJshift6 == 0.0)
					{
						double fac=double(ci.numMolecules())*double(cj.numMolecules())*eps24;
						if(tau1==0. && tau2==0.)
						{
							UpotCorrLJ+=fac*(TICCu(-6,cutoffRadiusLJ,sig2)-TICCu(-3,cutoffRadiusLJ,sig2));
							VirialCorrLJ+=fac*(TICCv(-6,cutoffRadiusLJ,sig2)-TICCv(-3,cutoffRadiusLJ,sig2));
						}
						else if(tau1!=0. && tau2!=0.)
						{
							UpotCorrLJ += fac*( TISSu(-6,cutoffRadiusLJ,sig2,tau1,tau2)
									- TISSu(-3,cutoffRadiusLJ,sig2,tau1,tau2) );
							VirialCorrLJ += fac*( TISSv(-6,cutoffRadiusLJ,sig2,tau1,tau2)
									- TISSv(-3,cutoffRadiusLJ,sig2,tau1,tau2) );
						}
						else {
							if(tau2==0.) 
								tau2=tau1;
							UpotCorrLJ+=fac*(TICSu(-6,cutoffRadiusLJ,sig2,tau2)-TICSu(-3,cutoffRadiusLJ,sig2,tau2));
							VirialCorrLJ+=fac*(TICSv(-6,cutoffRadiusLJ,sig2,tau2)-TICSv(-3,cutoffRadiusLJ,sig2,tau2));
						}
					}
				}
			}
		}
	}

	double fac=M_PI*_globalRho/(3.*_globalNumMolecules);
	UpotCorrLJ*=fac;
	VirialCorrLJ*=-fac;

	double epsRFInvrc3=2.*(_epsilonRF-1.)/((cutoffRadius*cutoffRadius*cutoffRadius)*(2.*_epsilonRF+1.));
	MySelbstTerm*=-0.5*epsRFInvrc3;

	_UpotCorr=UpotCorrLJ+MySelbstTerm;
	_VirialCorr=VirialCorrLJ+3.*MySelbstTerm;

	global_log->info() << "Far field terms: U_pot_correction  = " << _UpotCorr << " virial_correction = " << _VirialCorr << endl;
}

void Domain::specifyComponentSet(unsigned cosetid, double v[3], double tau, double ainit[3], double timestep)
{
	this->_localN[cosetid] = 0;
	for(unsigned d = 0; d < 3; d++)
	{
		this->_universalAdditionalAcceleration[d][cosetid] = ainit[d];
		this->_localVelocitySum[d][cosetid] = 0.0;
	}
	this->_universalTau[cosetid] = tau;
	if(this->_localRank == 0)
	{
		this->_globalN[cosetid] = 0;
		for(unsigned d = 0; d < 3; d++)
		{
			this->_globalTargetVelocity[d][cosetid] = v[d];
			this->_globalVelocitySum[d][cosetid] = 0.0;
		}
		this->_globalVelocityQueuelength[cosetid] = (unsigned)ceil(
				sqrt(this->_universalTau[cosetid] / (timestep*this->_universalConstantAccelerationTimesteps))
				);
		global_log->info() << "coset " << cosetid << " will receive "
			<< _globalVelocityQueuelength[cosetid] << " velocity queue entries." << endl;
	}
}

/*
 * diese Version beschleunigt in alle Raumrichtungen
 */
void Domain::determineAdditionalAcceleration
(
 DomainDecompBase* domainDecomp,
 ParticleContainer* molCont, double dtConstantAcc )
{
	for( map<unsigned, double>::iterator uAAit = _universalAdditionalAcceleration[0].begin();
			uAAit != _universalAdditionalAcceleration[0].end();
			uAAit++ )
	{
		this->_localN[uAAit->first] = 0;
		for(unsigned d = 0; d < 3; d++)
			this->_localVelocitySum[d][uAAit->first] = 0.0;
	}
	for(Molecule* thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
	{
		unsigned cid = thismol->componentid();
		map<unsigned, unsigned>::iterator uCSIDit = this->_universalComponentSetID.find(cid);
		if(uCSIDit == _universalComponentSetID.end()) continue;
		unsigned cosetid = uCSIDit->second;
		this->_localN[cosetid]++;
		for(unsigned d = 0; d < 3; d++)
			this->_localVelocitySum[d][cosetid] += thismol->v(d);
	}

	// domainDecomp->collectCosetVelocity(&_localN, _localVelocitySum, &_globalN, _globalVelocitySum);
	//
	domainDecomp->collCommInit( 4 * _localN.size() );
	for( map<unsigned, unsigned long>::iterator lNit = _localN.begin(); lNit != _localN.end(); lNit++ )
	{
		domainDecomp->collCommAppendUnsLong(lNit->second);
		for(int d = 0; d < 3; d++)
			domainDecomp->collCommAppendDouble( this->_localVelocitySum[d][lNit->first]);
	}
	domainDecomp->collCommAllreduceSum();
	for( map<unsigned, unsigned long>::iterator lNit = _localN.begin(); lNit != _localN.end(); lNit++ )
	{
		_globalN[lNit->first] = domainDecomp->collCommGetUnsLong();
		for(int d = 0; d < 3; d++)
			_globalVelocitySum[d][lNit->first] = domainDecomp->collCommGetDouble();
	}
	domainDecomp->collCommFinalize();

	map<unsigned, long double>::iterator gVSit;
	if(!this->_localRank)
	{
		for(gVSit = _globalVelocitySum[0].begin(); gVSit != _globalVelocitySum[0].end(); gVSit++)
		{
#ifndef NDEBUG
			global_log->debug() << "required entries in velocity queue: " << _globalVelocityQueuelength[gVSit->first] << endl;
			global_log->debug() << "entries in velocity queue: " << _globalPriorVelocitySums[0][gVSit->first].size() << endl;
#endif
			for(unsigned d = 0; d < 3; d++)
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
			for(unsigned d = 0; d < 3; d++)
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
	for( map<unsigned, unsigned long>::iterator lNit = _localN.begin();
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
	for( map<unsigned, unsigned long>::iterator lNit = _localN.begin();
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

double Domain::getDirectedVelocity(unsigned cosetid)
{
	double vv = 0.0;
	if(!this->_localRank)
	{
		for(unsigned d = 0; d < 3; d++)
		{
			double vd = this->_globalVelocitySum[d][cosetid] / this->_globalN[cosetid];
			vv += vd*vd;
		}
	}
	return sqrt(vv);
}
double Domain::getDirectedVelocity(unsigned cosetid, unsigned d)
{
	if(!this->_localRank) 
		return this->_globalVelocitySum[d][cosetid] / this->_globalN[cosetid];
	else return 0.0;
}

double Domain::getUniformAcceleration(unsigned cosetid)
{
	double aa = 0.0;
	for(unsigned d = 0; d < 3; d++)
		aa += _universalAdditionalAcceleration[d][cosetid] * _universalAdditionalAcceleration[d][cosetid];
	return sqrt(aa);
}
double Domain::getUniformAcceleration(unsigned cosetid, unsigned d)
{
	return this->_universalAdditionalAcceleration[d][cosetid];
}

double Domain::getMissingVelocity(unsigned cosetid, unsigned d)
{
	double vd = this->_globalVelocitySum[d][cosetid] / this->_globalN[cosetid];
	return this->_globalTargetVelocity[d][cosetid] - vd;
}

void Domain::setupProfile(unsigned xun, unsigned yun, unsigned zun)
{
	this->_universalNProfileUnits[0] = xun;
	this->_universalNProfileUnits[1] = yun;
	this->_universalNProfileUnits[2] = zun;
	for(unsigned d = 0; d < 3; d++)
	{
		_universalInvProfileUnit[d] = _universalNProfileUnits[d] / _globalLength[d];
	}
	this->resetProfile();
}

void Domain::considerComponentInProfile(int cid)
{
	this->_universalProfiledComponents[cid] = true;
}

void Domain::recordProfile(ParticleContainer* molCont)
{
	int cid;
	unsigned xun, yun, zun, unID;
	double mv2, Iw2;
	for(Molecule* thismol = molCont->begin(); thismol != molCont->end(); thismol = molCont->next())
	{
		cid = thismol->componentid();
		if(this->_universalProfiledComponents[cid])
		{
			xun = (unsigned)floor(thismol->r(0) * this->_universalInvProfileUnit[0]);
			yun = (unsigned)floor(thismol->r(1) * this->_universalInvProfileUnit[1]);
			zun = (unsigned)floor(thismol->r(2) * this->_universalInvProfileUnit[2]);
			unID = xun * this->_universalNProfileUnits[1] * this->_universalNProfileUnits[2]
				+ yun * this->_universalNProfileUnits[2] + zun;
			this->_localNProfile[unID] += 1.0;
			for(int d=0; d<3; d++) this->_localvProfile[d][unID] += thismol->v(d);
			this->_localDOFProfile[unID] += 3.0 + (long double)(this->_components[cid].rot_dof());

			// record _twice_ the total (ordered + unordered) kinetic energy
			mv2 = 0.0;
			Iw2 = 0.0;
			thismol->calculate_mv2_Iw2(mv2, Iw2);
			this->_localKineticProfile[unID] += mv2+Iw2;
		}
	}
	this->_globalAccumulatedDatasets++;
}

void Domain::collectProfile(DomainDecompBase* dode)
{
	unsigned unIDs = this->_universalNProfileUnits[0] * this->_universalNProfileUnits[1]
		* this->_universalNProfileUnits[2];
	dode->collCommInit(6*unIDs);
	for(unsigned unID = 0; unID < unIDs; unID++)
	{
		dode->collCommAppendLongDouble(this->_localNProfile[unID]);
		for(int d=0; d<3; d++)
			dode->collCommAppendLongDouble(_localvProfile[d][unID]);
		dode->collCommAppendLongDouble(this->_localDOFProfile[unID]);
		dode->collCommAppendLongDouble(_localKineticProfile[unID]);
	}
	dode->collCommAllreduceSum();
	for(unsigned unID = 0; unID < unIDs; unID++)
	{
		_globalNProfile[unID] = (double)dode->collCommGetLongDouble();
		for(int d=0; d<3; d++)
			this->_globalvProfile[d][unID]
				= (double)dode->collCommGetLongDouble();
		this->_globalDOFProfile[unID]
			= (double)dode->collCommGetLongDouble();
		this->_globalKineticProfile[unID]
			= (double)dode->collCommGetLongDouble();
	}
	dode->collCommFinalize();
}

void Domain::outputProfile(const char* prefix)
{
	if(this->_localRank) return;

	string vzpryname(prefix);
	string Tpryname(prefix);
	string rhpryname(prefix);
	rhpryname += ".rhpry";
	vzpryname += ".vzpry";
	Tpryname += ".Tpry";
	ofstream rhpry(rhpryname.c_str());
	ofstream vzpry(vzpryname.c_str());
	ofstream Tpry(Tpryname.c_str());
	if (!(vzpry && Tpry && rhpry))
	{
		return;
	}
	rhpry.precision(4);
	rhpry << "# y\trho\ttotal DOF\n# \n";
	vzpry.precision(4);
	vzpry << "# y\tvz\tv\n# \n";
	Tpry.precision(5);
	Tpry << "# y\t2Ekin/#DOF\n# \n";

	double layerVolume = this->_globalLength[0] * this->_globalLength[1] * this->_globalLength[2]
		/ this->_universalNProfileUnits[1];
	for(unsigned y = 0; y < this->_universalNProfileUnits[1]; y++)
	{
		double yval = (y + 0.5) / this->_universalInvProfileUnit[1];

		long double Ny = 0.0;
		long double DOFy = 0.0;
		long double twoEkiny = 0.0;
		long double velocitysumy[3];
		for(unsigned d = 0; d < 3; d++) velocitysumy[d] = 0.0;
		for(unsigned x = 0; x < this->_universalNProfileUnits[0]; x++)
		{
			for(unsigned z = 0; z < this->_universalNProfileUnits[2]; z++)
			{
				unsigned unID = x * this->_universalNProfileUnits[1] * this->_universalNProfileUnits[2]
					+ y * this->_universalNProfileUnits[2] + z;
				Ny += this->_globalNProfile[unID];
				DOFy += this->_globalDOFProfile[unID];
				twoEkiny += this->_globalKineticProfile[unID];
				for(unsigned d = 0; d < 3; d++) velocitysumy[d] += this->_globalvProfile[d][unID];
			}
		}

		if(Ny >= 64.0)
		{
			double vvdir = 0.0;
			for(unsigned d = 0; d < 3; d++)
			{
				double vd = velocitysumy[d] / Ny;
				vvdir += vd*vd;
			}
			rhpry << yval << "\t" << (Ny / (layerVolume * this->_globalAccumulatedDatasets))
				<< "\t" << DOFy << "\n";
			vzpry << yval << "\t" << (velocitysumy[2] / Ny) << "\t" << sqrt(vvdir) << "\n";
			Tpry << yval << "\t" << (twoEkiny / DOFy) << "\n";
		}
		else
		{
			rhpry << yval << "\t0.000\t" << DOFy << "\n";
		}
	}

	rhpry.close();
	vzpry.close();
	Tpry.close();
}

void Domain::resetProfile()
{
	unsigned unIDs = this->_universalNProfileUnits[0] * this->_universalNProfileUnits[1]
		* this->_universalNProfileUnits[2];
	for(unsigned unID = 0; unID < unIDs; unID++)
	{
		this->_localNProfile[unID] = 0.0;
		this->_globalNProfile[unID] = 0.0;
		for(int d=0; d<3; d++)
		{
			this->_localvProfile[d][unID] = 0.0;
			this->_globalvProfile[d][unID] = 0.0;
		}
		this->_localDOFProfile[unID] = 0.0;
		this->_globalDOFProfile[unID] = 0.0;
		this->_localKineticProfile[unID] = 0.0;
		this->_globalKineticProfile[unID] = 0.0;
	}
	this->_globalAccumulatedDatasets = 0;
}

void Domain::setupRDF(double interval, unsigned bins)
{
	unsigned numc = this->_components.size();

	this->_doCollectRDF = true;
	this->_universalInterval = interval;
	this->_universalBins = bins;
	this->_universalRDFTimesteps = -1;
	this->_universalAccumulatedTimesteps = 0;
	this->ddmax = interval*interval*bins*bins;

	this->_localCtr = new unsigned long[numc];
	this->_globalCtr = new unsigned long[numc];
	this->_globalAccumulatedCtr = new unsigned long[numc];
	this->_localDistribution = new unsigned long**[numc];
	this->_globalDistribution = new unsigned long**[numc];
	this->_globalAccumulatedDistribution = new unsigned long**[numc];
	for(unsigned i = 0; i < numc; i++)
	{
		this->_localCtr[i] = 0;
		this->_globalCtr[i] = 0;
		this->_globalAccumulatedCtr[i] = 0;

		this->_localDistribution[i] = new unsigned long*[numc-i];
		this->_globalDistribution[i] = new unsigned long*[numc-i];
		this->_globalAccumulatedDistribution[i] = new unsigned long*[numc-i];

		for(unsigned k=0; i+k < numc; k++)
		{
			this->_localDistribution[i][k] = new unsigned long[bins];
			this->_globalDistribution[i][k] = new unsigned long[bins];
			this->_globalAccumulatedDistribution[i][k] = new unsigned long[bins];

			for(unsigned l=0; l < bins; l++)
			{
				this->_localDistribution[i][k][l] = 0;
				this->_globalDistribution[i][k][l] = 0;
				this->_globalAccumulatedDistribution[i][k][l] = 0;
			}
		}
	}
}

void Domain::resetRDF()
{
	this->_universalRDFTimesteps = 0;
	for(unsigned i=0; i < this->_components.size(); i++)
	{
		this->_localCtr[i] = 0;
		this->_globalCtr[i] = 0;
		for(unsigned k=0; i+k < this->_components.size(); k++)
		{
			for(unsigned l=0; l < this->_universalBins; l++)
			{
				this->_localDistribution[i][k][l] = 0;
				this->_globalDistribution[i][k][l] = 0;
			}
		}
	}
}

void Domain::accumulateRDF()
{
	if(0 >= this->_universalRDFTimesteps) return;
	this->_universalAccumulatedTimesteps += this->_universalRDFTimesteps;
	if(!this->_localRank)
	{
		for(unsigned i=0; i < this->_components.size(); i++)
		{
			this->_globalAccumulatedCtr[i] += this->_globalCtr[i];
			for(unsigned k=0; i+k < this->_components.size(); k++)
				for(unsigned l=0; l < this->_universalBins; l++)
					this->_globalAccumulatedDistribution[i][k][l]
						+= this->_globalDistribution[i][k][l];
		}
	}
}

void Domain::collectRDF(DomainDecompBase* dode)
{
	unsigned nb = this->_universalBins;
	unsigned nc = this->_components.size();
	dode->collCommInit(nc + nb*nc*(nc+1)/2);
	for(unsigned i=0; i < nc; i++)
	{
		dode->collCommAppendUnsLong(this->_localCtr[i]);
		for(unsigned k=0; i+k < nc; k++)
		{
			for(unsigned l=0; l < nb; l++)
			{
				dode->collCommAppendUnsLong(_localDistribution[i][k][l]);
			}
		}
	}
	dode->collCommAllreduceSum();
	for(unsigned i=0; i < nc; i++)
	{
		this->_globalCtr[i] = dode->collCommGetUnsLong();
		for(unsigned k=0; i+k < nc; k++)
		{
			for(unsigned l=0; l < nb; l++)
			{
				_globalDistribution[i][k][l] = dode->collCommGetUnsLong();
			}
		}
	}
	dode->collCommFinalize();
}

void Domain::outputRDF(const char* prefix, unsigned i, unsigned j)
{
	/* Output only from process with rank 0 */
	if( 0 == _localRank ) 
		return;

	/* FIXME: MINIMAL >= xyz ?! */
	if(RDF_MINIMAL_OUTPUT_STEPS >= _universalRDFTimesteps) 
		return;

	string rdfname(prefix);
	rdfname += ".rdf";
	ofstream rdfout(rdfname.c_str());
	if (!rdfout) return;

	double V = _globalLength[0] * _globalLength[1] * _globalLength[2];
	double N_i = _globalCtr[i] / _universalRDFTimesteps;
	double N_Ai = _globalAccumulatedCtr[i] / _universalAccumulatedTimesteps;
	double N_j = _globalCtr[j] / _universalRDFTimesteps;
	double N_Aj = _globalAccumulatedCtr[j] / _universalAccumulatedTimesteps;
	double rho_i = N_i / V;
	double rho_Ai = N_Ai / V;
	double rho_j = N_j / V;
	double rho_Aj = N_Aj / V;

	rdfout.precision(5);
	rdfout << "# r\tcurr.\taccu.\t\tdV\tNpair(curr.)\tNpair(accu.)\t\tnorm(curr.)\tnorm(accu.)\n";
	rdfout << "# \n# ctr_i: " << _globalCtr[i] << "\n# ctr_j: " << _globalCtr[j]
		<< "\n# V: " << V << "\n# _universalRDFTimesteps: " << _universalRDFTimesteps
		<< "\n# _universalAccumulatedTimesteps: " << _universalAccumulatedTimesteps
		<< "\n# rho_i: " << rho_i << " (acc. " << rho_Ai << ")"
		<< "\n# rho_j: " << rho_j << " (acc. " << rho_Aj << ")"
		<< "\n# \n";

	for(unsigned l=0; l < this->_universalBins; l++)
	{
		double rmin = l * _universalInterval;
		double rmid = (l+0.5) * _universalInterval;
		double rmax = (l+1.0) * _universalInterval;
		double r3min = rmin*rmin*rmin;
		double r3max = rmax*rmax*rmax;
		/* FIXME: uncommented constant */
		double dV = 4.1887902 * (r3max - r3min);

		unsigned long N_pair = _globalDistribution[i][j-i][l] / _universalRDFTimesteps;
		unsigned long N_Apair = _globalAccumulatedDistribution[i][j-i][l]
			/ _universalAccumulatedTimesteps;
		double N_pair_norm;
		double N_Apair_norm;
		if(i == j)
		{
			N_pair_norm = 0.5*N_i*rho_i*dV;
			N_Apair_norm = 0.5*N_Ai*rho_Ai*dV;
		}
		else
		{
			N_pair_norm = N_i*rho_j*dV;
			N_Apair_norm = N_Ai*rho_Aj*dV;
		}

		rdfout << rmid << "\t" << N_pair/N_pair_norm
			<< "\t" << N_Apair/N_Apair_norm
			<< "\t\t" << dV << "\t" << N_pair << "\t" << N_Apair
			<< "\t\t" << N_pair_norm << "\t" << N_Apair_norm << "\n";
	}
	rdfout.close();
}

void Domain::Nadd(unsigned cid, int N, int localN)
{
	this->_components[cid].Nadd(N);
	this->_globalNumMolecules += N;
	this->_localRotationalDOF[0] += localN * this->_components[cid].rot_dof();
	this->_universalRotationalDOF[0] += N * this->_components[cid].rot_dof();
	if( (this->_componentwiseThermostat)
			&& (this->_componentToThermostatIdMap[cid] > 0) )
	{
		int thid = this->_componentToThermostatIdMap[cid];
		this->_localThermostatN[thid] += localN;
		this->_universalThermostatN[thid] += N;
		this->_localRotationalDOF[thid] += localN * this->_components[cid].rot_dof();
		this->_universalRotationalDOF[thid] += N * this->_components[cid].rot_dof();
	}
	this->_localThermostatN[0] += localN;
	this->_universalThermostatN[0] += N;
	this->_localRotationalDOF[0] += localN * this->_components[cid].rot_dof();
	this->_universalRotationalDOF[0] += N * this->_components[cid].rot_dof();
}

void Domain::evaluateRho(
		unsigned long localN, DomainDecompBase* domainDecomp
		) {
	domainDecomp->collCommInit(1);
	domainDecomp->collCommAppendUnsLong(localN);
	domainDecomp->collCommAllreduceSum();
	this->_globalNumMolecules = domainDecomp->collCommGetUnsLong();
	domainDecomp->collCommFinalize();

	this->_globalRho = this->_globalNumMolecules /
		(this->_globalLength[0] * this->_globalLength[1] * this->_globalLength[2]);
}

void Domain::setTargetTemperature(int thermostat, double targetT)
{
	if(thermostat < 0)
	{
		global_log->warning() << "Warning: thermostat \'" << thermostat << "\' (T = "
			<< targetT << ") will be ignored." << endl;
		return;
	}

	this->_universalTargetTemperature[thermostat] = targetT;
	if(!(this->_universalUndirectedThermostat[thermostat] == true))
		this->_universalUndirectedThermostat[thermostat] = false;

	/* FIXME: Substantial change in program behaviour without any info to the user! */
	if(thermostat == 0) this->disableComponentwiseThermostat();
	if(thermostat >= 1)
	{
		if(!this->_componentwiseThermostat)
		{
			/* FIXME: Substantial change in program behaviour without any info to the user! */
			this->_componentwiseThermostat = true;
			this->_universalTargetTemperature.erase(0);
			this->_universalUndirectedThermostat.erase(0);
			for(int d=0; d < 3; d++) this->_universalThermostatDirectedVelocity[d].erase(0);
			for( vector<Component>::iterator tc = this->_components.begin();
					tc != this->_components.end();
					tc ++ )
				if(!(this->_componentToThermostatIdMap[ tc->ID() ] > 0))
					this->_componentToThermostatIdMap[ tc->ID() ] = -1;
		}
	}
}

void Domain::enableComponentwiseThermostat()
{
	if(this->_componentwiseThermostat) return;

	this->_componentwiseThermostat = true;
	this->_universalTargetTemperature.erase(0);
	for( vector<Component>::iterator tc = this->_components.begin();
			tc != this->_components.end();
			tc ++ )
		if(!(this->_componentToThermostatIdMap[ tc->ID() ] > 0))
			this->_componentToThermostatIdMap[ tc->ID() ] = -1;
}

void Domain::enableUndirectedThermostat(int tst)
{
	this->_universalUndirectedThermostat[tst] = true;
	for(int d=0; d < 3; d++)
	{
		this->_universalThermostatDirectedVelocity[d][tst] = 0.0;
		this->_localThermostatDirectedVelocity[d][tst] = 0.0;
	}
}

void Domain::setGlobalTemperature(double temp)
{
	this->disableComponentwiseThermostat();
	this->_universalTargetTemperature[0] = temp;
}

vector<double> & Domain::getmixcoeff() { return _mixcoeff; }

double Domain::getepsilonRF() const { return _epsilonRF; }

void Domain::setepsilonRF(double erf) { _epsilonRF = erf; }

unsigned long Domain::getglobalNumMolecules() const { return _globalNumMolecules; }

void Domain::setglobalNumMolecules(unsigned long glnummol) { _globalNumMolecules = glnummol; }

int Domain::getlocalRank(){ return _localRank;}

unsigned long Domain::getinpversion(){ return _inpversion;}

void Domain::setinpversion(unsigned long inpv){ _inpversion = inpv;}

double Domain::getglobalRho(){ return _globalRho;}

void Domain::setglobalRho(double grho){ _globalRho = grho;}

unsigned long Domain::getglobalRotDOF()
{
	return this->_universalRotationalDOF[0]; 
}

void Domain::setglobalRotDOF(unsigned long grotdof)
{
	this->_universalRotationalDOF[0] = grotdof;
}

void Domain::setGlobalLength(int index, double length) {
	_globalLength[index] = length;
}

void Domain::record_cv()
{
	if(_localRank != 0) return;

	this->_globalUSteps ++;
	this->_globalSigmaU += this->_globalUpot;
	this->_globalSigmaUU += this->_globalUpot*_globalUpot;
}

double Domain::cv()
{
	if((_localRank != 0) || (_globalUSteps == 0)) return 0.0;

	double id = 1.5 + 0.5*_universalRotationalDOF[0]/_globalNumMolecules;
	double conf = (_globalSigmaUU - _globalSigmaU*_globalSigmaU/_globalUSteps)
		/ (_globalUSteps * _globalNumMolecules * _globalTemperatureMap[0] * _globalTemperatureMap[0]);

	return id + conf;
}

