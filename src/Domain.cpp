
#include <iostream>
#include <string>
#include <cmath>

#include "Domain.h"

#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "ensemble/GrandCanonical.h"
#include "ensemble/PressureGradient.h"
//#include "CutoffCorrections.h"
#include "Simulation.h"
#include "ensemble/EnsembleBase.h"

#include "utils/FileUtils.h"
#include "utils/Logger.h"
using Log::global_log;

using namespace std;


Domain::Domain(int rank, PressureGradient* pg){
	_localRank = rank;
	_localUpot = 0;
	_localVirial = 0;
	_globalUpot = 0;
	_globalVirial = 0;
	_globalRho = 0;
	this->_Gamma = map<unsigned, double>();

	this->_universalPG = pg;

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

	this->_universalNVE = false;
	this->_globalUSteps = 0;
	this->_globalSigmaU = 0.0;
	this->_globalSigmaUU = 0.0;
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

    // explosion heuristics, NOTE: turn off when using slab thermostat
    _bDoExplosionHeuristics = true;
}

void Domain::readXML(XMLfileUnits& xmlconfig) {
	string originalpath = xmlconfig.getcurrentnodepath();

	/* volume */
	if ( xmlconfig.changecurrentnode( "volume" )) {
		std::string type;
		xmlconfig.getNodeValue( "@type", type );
		global_log->info() << "Volume type: " << type << endl;
		if( type == "box" ) {
			xmlconfig.getNodeValueReduced( "lx", _globalLength[0] );
			xmlconfig.getNodeValueReduced( "ly", _globalLength[1] );
			xmlconfig.getNodeValueReduced( "lz", _globalLength[2] );
			global_log->info() << "Box size: " << _globalLength[0] << ", "
				<< _globalLength[1] << ", "
				<< _globalLength[2] << endl;
		}
		else {
			global_log->error() << "Unsupported volume type " << type << endl;
		}
	}
	xmlconfig.changecurrentnode(originalpath);

	/* temperature */
	double temperature = 0.;
	xmlconfig.getNodeValueReduced("temperature", temperature);
	setGlobalTemperature(temperature);
	global_log->info() << "Temperature: " << temperature << endl;
	xmlconfig.changecurrentnode(originalpath);

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

void Domain::setLocalSumIw2(double sumIw2, int thermostat)
{
	_local2KERot[thermostat] = sumIw2;
}

double Domain::getGlobalPressure()
{
	double globalTemperature = _globalTemperatureMap[0];
	return globalTemperature * _globalRho + _globalRho * getAverageGlobalVirial()/3.;
}

double Domain::getAverageGlobalVirial() const { return _globalVirial/_globalNumMolecules; }

double Domain::getAverageGlobalUpot() const { return getGlobalUpot()/_globalNumMolecules; }
double Domain::getGlobalUpot() const { return _globalUpot; }

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

	/* FIXME stuff for the ensemble class */
	domainDecomp->collCommInit(2);
	domainDecomp->collCommAppendDouble(Upot);
	domainDecomp->collCommAppendDouble(Virial);
	domainDecomp->collCommAllreduceSum();
	Upot = domainDecomp->collCommGetDouble();
	Virial = domainDecomp->collCommGetDouble();
	domainDecomp->collCommFinalize();

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
		global_log->debug() << "* applying a component-wise thermostat" << endl;
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
		global_log->debug() << "[ thermostat ID " << thermit->first << "]\tN = " << numMolecules << "\trotDOF = " << rotDOF
			<< "\tmv2 = " <<  summv2 << "\tIw2 = " << sumIw2 << endl;

		this->_universalThermostatN[thermit->first] = numMolecules;
		this->_universalRotationalDOF[thermit->first] = rotDOF;
		mardyn_assert((summv2 > 0.0) || (numMolecules == 0));

		/* calculate the temperature of the entire system */
		if(numMolecules > 0)
			_globalTemperatureMap[thermit->first] =
				(summv2 + sumIw2) / (double)(3*numMolecules + rotDOF);
		else
			_globalTemperatureMap[thermit->first] = _universalTargetTemperature[thermit->first];

		double Ti = Tfactor * _universalTargetTemperature[thermit->first];
		if((Ti > 0.0) && (numMolecules > 0) && !_universalNVE)
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
				&& (0 >= _universalSelectiveThermostatError)  && _bDoExplosionHeuristics == true)
		{
			global_log->warning() << "Explosion!" << endl;
			global_log->debug() << "Selective thermostat will be applied to set " << thermit->first
				<< " (beta_trans = " << this->_universalBTrans[thermit->first]
				<< ", beta_rot = " << this->_universalBRot[thermit->first] << "!)" << endl;
			int rot_dof;
			double Utrans, Urot;
			double limit_energy =  KINLIMIT_PER_T * Ti;
			double limit_rot_energy;
			double vcorr, Dcorr;
			ParticleIterator tM;
			for( tM = particleContainer->iteratorBegin();
					tM != particleContainer->iteratorEnd();
					++tM)
			{
				Utrans = tM->U_trans();
				if(Utrans > limit_energy)
				{
					vcorr = sqrt(limit_energy / Utrans);
					global_log->debug() << ": v(m" << tM->id() << ") *= " << vcorr << endl;
					tM->scale_v(vcorr);
					tM->scale_F(vcorr);
				}

				rot_dof = tM->component()->getRotationalDegreesOfFreedom();
				if(rot_dof > 0)
				{
					limit_rot_energy = 3.0*rot_dof * Ti;
					Urot = tM->U_rot();
					if(Urot > limit_rot_energy)
					{
						Dcorr = sqrt(limit_rot_energy / Urot);
						global_log->debug() << "D(m" << tM->id() << ") *= " << Dcorr << endl;
						tM->scale_D(Dcorr);
						tM->scale_M(Dcorr);
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
	ParticleIterator tM;
	if(this->_componentwiseThermostat)
	{
		for( map<int, bool>::iterator thit = _universalUndirectedThermostat.begin();
				thit != _universalUndirectedThermostat.end();
				thit ++ )
		{
			if(thit->second)
				for(int d=0; d < 3; d++) _localThermostatDirectedVelocity[d][thit->first] = 0.0;
		}
		for(tM = partCont->iteratorBegin(); tM != partCont->iteratorEnd(); ++tM)
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
		for(tM = partCont->iteratorBegin(); tM != partCont->iteratorEnd(); ++tM)
		{
			for(int d=0; d < 3; d++)
				_localThermostatDirectedVelocity[d][0] += tM->v(d);
		}
	}
}

void Domain::calculateVelocitySums(ParticleContainer* partCont)
{
	ParticleIterator tM;
	if(this->_componentwiseThermostat)
	{
		for(tM = partCont->iteratorBegin(); tM != partCont->iteratorEnd(); ++tM)
		{
			int cid = tM->componentid();
			int thermostat = this->_componentToThermostatIdMap[cid];
			this->_localThermostatN[thermostat]++;
			this->_localRotationalDOF[thermostat] += tM->component()->getRotationalDegreesOfFreedom();
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
		for(tM = partCont->iteratorBegin(); tM != partCont->iteratorEnd(); ++tM)
		{
			this->_localThermostatN[0]++;
			this->_localRotationalDOF[0] += tM->component()->getRotationalDegreesOfFreedom();
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
		global_log->debug() << "      * N = " << this->_localThermostatN[0]
			<< "rotDOF = " << this->_localRotationalDOF[0] << "   mv2 = "
			<< _local2KETrans[0] << " Iw2 = " << _local2KERot[0] << endl;
	}
}

void Domain::writeCheckpointHeader(string filename,
		ParticleContainer* particleContainer, const DomainDecompBase* domainDecomp, double currentTime) {
		domainDecomp->assertDisjunctivity(particleContainer);
		/* Rank 0 writes file header */
		if(0 == this->_localRank) {
			ofstream checkpointfilestream(filename.c_str());
			checkpointfilestream << "mardyn trunk " << CHECKPOINT_FILE_VERSION;
			checkpointfilestream << "\n";  // store default format flags
			ios::fmtflags f( checkpointfilestream.flags() );
			checkpointfilestream << "currentTime\t"  << FORMAT_SCI_MAX_DIGITS << currentTime << "\n"; //edited by Michaela Heier
			checkpointfilestream.flags(f);  // restore default format flags
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
			//checkpointfilestream << "# rc\t" << global_simulation->getcutoffRadius() << "\n";
	        checkpointfilestream << "# \n# Please address your questions and suggestions to\n# the ls1 mardyn contact point: <contact@ls1-mardyn.de>.\n# \n";
	#endif
			/* by Stefan Becker: the output line "I ..." causes an error: the restart run does not start!!!
			if(this->_globalUSteps > 1)
			if(this->_globalUSteps > 1)
			{
				checkpointfilestream << setprecision(13);
				checkpointfilestream << " I\t" << this->_globalUSteps << " "
					<< this->_globalSigmaU << " " << this->_globalSigmaUU << "\n";
				checkpointfilestream << setprecision(8);
			}
			*/
			vector<Component>* components = _simulation.getEnsemble()->getComponents();
			checkpointfilestream << " NumberOfComponents\t" << components->size() << endl;
			for(vector<Component>::const_iterator pos=components->begin();pos!=components->end();++pos){
				pos->write(checkpointfilestream);
			}
			unsigned int numperline=_simulation.getEnsemble()->getComponents()->size();
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
			map<unsigned, unsigned> componentSets = this->_universalPG->getComponentSets();
			for( map<unsigned, unsigned>::const_iterator uCSIDit = componentSets.begin();
					uCSIDit != componentSets.end();
					uCSIDit++ )
			{
				if(uCSIDit->first > 100) continue;
				checkpointfilestream << " S\t" << 1+uCSIDit->first << "\t" << uCSIDit->second << "\n";
			}
			map<unsigned, double> tau = this->_universalPG->getTau();
			for( map<unsigned, double>::const_iterator gTit = tau.begin();
					gTit != tau.end();
					gTit++ )
			{
				unsigned cosetid = gTit->first;
				double* ttargetv = this->_universalPG->getTargetVelocity(cosetid);
				double* tacc = this->_universalPG->getAdditionalAcceleration(cosetid);
				checkpointfilestream << " A\t" << cosetid << "\t"
					<< ttargetv[0] << " " << ttargetv[1] << " " << ttargetv[2] << "\t"
					<< gTit->second << "\t"
					<< tacc[0] << " " << tacc[1] << " " << tacc[2] << "\n";
				delete ttargetv;
				delete tacc;
			}
			for( map<int, bool>::iterator uutit = this->_universalUndirectedThermostat.begin();
					uutit != this->_universalUndirectedThermostat.end();
					uutit++ )
			{
				if(0 > uutit->first) continue;
				if(uutit->second) checkpointfilestream << " U\t" << uutit->first << "\n";
			}
			checkpointfilestream << " NumberOfMolecules\t" << _globalNumMolecules << endl;

			checkpointfilestream << " MoleculeFormat\t" << Molecule::getWriteFormat() << endl;
			checkpointfilestream.close();
		}

}

void Domain::writeCheckpoint(string filename,
		ParticleContainer* particleContainer, const DomainDecompBase* domainDecomp, double currentTime,
		bool binary) {

#ifdef MARDYN_WR
	global_log->warning() << "The checkpoints are not adapted for WR-mode. Velocity will be one half-timestep ahead!" << std::endl;
	global_log->warning() << "See Domain::writeCheckpoint() for a suggested workaround." << std::endl;
	//TODO: desired correctness (compatibility to normal mode) should be achievable by:
	// 1. integrating positions by half a timestep forward (+ delta T / 2)
	// 2. writing the checkpoint (with currentTime + delta T ? )
	// 3. integrating positions by half a timestep backward (- delta T / 2)
#endif

	if (binary == true) {
		this->writeCheckpointHeader((filename + ".header.xdr"), particleContainer, domainDecomp, currentTime);
	} else {
		this->writeCheckpointHeader(filename, particleContainer, domainDecomp, currentTime);
	}

	if (binary == true) {
		domainDecomp->writeMoleculesToFile(filename, particleContainer, true);
	} else {
		domainDecomp->writeMoleculesToFile(filename, particleContainer, false);
	}
}


void Domain::initParameterStreams(double cutoffRadius, double cutoffRadiusLJ){
	_comp2params.initialize(*(_simulation.getEnsemble()->getComponents()), _mixcoeff, _epsilonRF, cutoffRadius, cutoffRadiusLJ);
}

/*void Domain::initFarFieldCorr(double cutoffRadius, double cutoffRadiusLJ) {
	double UpotCorrLJ=0.;
	double VirialCorrLJ=0.;
	double MySelbstTerm=0.;
	vector<Component>* components = _simulation.getEnsemble()->getComponents();
	unsigned int numcomp=components->size();
	for(unsigned int i=0;i<numcomp;++i) {
		Component& ci=(*components)[i];
		unsigned int numljcentersi=ci.numLJcenters();
		unsigned int numchargesi = ci.numCharges();
		unsigned int numdipolesi=ci.numDipoles();

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
		MySelbstTerm += my2 * ci.getNumMolecules();

		for(unsigned int j=0;j<numcomp;++j) {
			Component& cj=(*components)[j];
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
						Simulation::exit(1);
					}
					double eps24;
					params >> eps24;
					double sig2;
					params >> sig2;
					double uLJshift6;
					params >> uLJshift6;  // 0 unless TRUNCATED_SHIFTED

					if(uLJshift6 == 0.0)
					{
						double fac=double(ci.getNumMolecules())*double(cj.getNumMolecules())*eps24;
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
}*/

void Domain::setupProfile(unsigned xun, unsigned yun, unsigned zun)
{
	this->_universalNProfileUnits[0] = xun;
	this->_universalNProfileUnits[1] = yun;
	this->_universalNProfileUnits[2] = zun;
	for(unsigned d = 0; d < 3; d++)
	{
		_universalInvProfileUnit[d] = (double)_universalNProfileUnits[d] / _globalLength[d];
	}
	this->resetProfile(true);
}

// author: Stefan Becker. Method called by Simulation::output() in order to decide wheter or not a cylindrical profile is to be written out,
//i.e. wheter the method outputCylProfile() (isCylindrical==true) or the method outputProfile() (isCylindrical==false) is called.
bool Domain::isCylindrical(){
	return this->_universalCylindricalGeometry;
}

void Domain::considerComponentInProfile(int cid)
{
	this->_universalProfiledComponents[cid] = true;
}

void Domain::recordProfile(ParticleContainer* molCont, bool virialProfile)
{
	int cid;
	unsigned xun, yun, zun;
	long int unID;
	double mv2, Iw2;
	unID = 0;
	unsigned lNin = 0;
	unsigned lNout = 0;
	for(ParticleIterator thismol = molCont->iteratorBegin(); thismol != molCont->iteratorEnd(); ++thismol)
	{
		cid = thismol->componentid();
		if(this->_universalProfiledComponents[cid])
		{

// by Stefan Becker: enquiry if(_universalCylindricalGeometry) ... implemented
// possible???: if no cylindrical profile is recorded => permanent calculation (before the "if..." ) of "distFor_unID" slows down the code
// if the profile is however recorded in cylindrical coordinates, the current implementation is fast (?)   => solution?

			if (this->_universalCylindricalGeometry) {
				double distFor_unID = pow((thismol->r(0) - this->_universalCentre[0]), 2.0)
						+ pow((thismol->r(2) - this->_universalCentre[2]), 2.0);
				if (distFor_unID <= this->_universalR2max) {
					unID = this->unID(thismol->r(0), thismol->r(1), thismol->r(2));
					mardyn_assert(unID >= 0);
					mardyn_assert(
							unID
									< (this->_universalNProfileUnits[0] * this->_universalNProfileUnits[1]
											* this->_universalNProfileUnits[2]));
				} else {
					lNout++;
					continue;
				}
			} else {
				xun = (unsigned) floor(thismol->r(0) * this->_universalInvProfileUnit[0]);
				yun = (unsigned) floor(thismol->r(1) * this->_universalInvProfileUnit[1]);
				zun = (unsigned) floor(thismol->r(2) * this->_universalInvProfileUnit[2]);
				unID = xun * this->_universalNProfileUnits[1] * this->_universalNProfileUnits[2]
						+ yun * this->_universalNProfileUnits[2] + zun;
			}

// @TODO: (by Stefan Becker)  differentiation of _localNProfile by the component number cid => _localNProfile[cid][unID]!!!
			lNin++;
			this->_localNProfile[unID] += 1.0;
			for (int d = 0; d < 3; d++) {
				this->_localvProfile[d][unID] += thismol->v(d);
			}
			this->_localDOFProfile[unID] += 3.0 + (long double) (thismol->component()->getRotationalDegreesOfFreedom());

			// record _twice_ the total (ordered + unordered) kinetic energy
			mv2 = 0.0;
			Iw2 = 0.0;
			thismol->calculate_mv2_Iw2(mv2, Iw2);
			this->_localKineticProfile[unID] += mv2 + Iw2;

			if (virialProfile) {
				// this->_localPDProfile[unID] += thismol->Vi(1)-0.5*(thismol->Vi(0)+thismol->Vi(2)); // unnecessary redundancy
				this->_localPXProfile[unID] += thismol->Vi(0);
				this->_localPYProfile[unID] += thismol->Vi(1);
				this->_localPZProfile[unID] += thismol->Vi(2);
			}
		}
	}
	this->_globalAccumulatedDatasets++;
#ifndef NDEBUG
	//cout << "Rank " << this->_localRank << " counted " << lNin << " molecules inside and " << lNout << " outside.\n";
	// cout << "Universal centre situated at (" << this->_universalCentre[0] << " / " << this->_universalCentre[1] << " / " << this->_universalCentre[2] << ").\n";
#endif
}

void Domain::collectProfile(DomainDecompBase* dode, bool virialProfile)
{
	unsigned unIDs = this->_universalNProfileUnits[0] * this->_universalNProfileUnits[1]
		* this->_universalNProfileUnits[2];
	if(virialProfile) dode->collCommInit(13*unIDs);
        else dode->collCommInit(10*unIDs);
	for(unsigned unID = 0; unID < unIDs; unID++)
	{
		dode->collCommAppendLongDouble(this->_localNProfile[unID]);
		for(int d=0; d<3; d++){
			dode->collCommAppendLongDouble(_localvProfile[d][unID]);

		}
		dode->collCommAppendLongDouble(this->_localDOFProfile[unID]);
		dode->collCommAppendLongDouble(_localKineticProfile[unID]);

                if(virialProfile)
                {
		   // dode->collCommAppendLongDouble(_localPDProfile[unID]); // unnecessary redundancy
		   dode->collCommAppendLongDouble(_localPXProfile[unID]);
		   dode->collCommAppendLongDouble(_localPYProfile[unID]);
		   dode->collCommAppendLongDouble(_localPZProfile[unID]);
                }

                dode->collCommAppendLongDouble(this->_localWidomProfile[unID]);
                dode->collCommAppendLongDouble(this->_localWidomInstances[unID]);
                dode->collCommAppendLongDouble(this->_localWidomProfileTloc[unID]);
                dode->collCommAppendLongDouble(this->_localWidomInstancesTloc[unID]);

	}
	dode->collCommAllreduceSum();
	for(unsigned unID = 0; unID < unIDs; unID++)
	{
		_universalNProfile[unID] = (double)dode->collCommGetLongDouble();
		for(int d=0; d<3; d++){
			this->_universalvProfile[d][unID]
				= (double)dode->collCommGetLongDouble();

		}
		this->_universalDOFProfile[unID]
			= (double)dode->collCommGetLongDouble();
		this->_universalKineticProfile[unID]
			= (double)dode->collCommGetLongDouble();

                if(virialProfile)
                {
		   // this->_universalPDProfile[unID]
			   // = (double)dode->collCommGetLongDouble(); // unnecessary redundancy
		   this->_universalPXProfile[unID]
			   = (double)dode->collCommGetLongDouble();
		   this->_universalPYProfile[unID]
			   = (double)dode->collCommGetLongDouble();
		   this->_universalPZProfile[unID]
			   = (double)dode->collCommGetLongDouble();
                }

		this->_globalWidomProfile[unID]
			= (double)dode->collCommGetLongDouble();
		this->_globalWidomInstances[unID]
			= (double)dode->collCommGetLongDouble();
		this->_globalWidomProfileTloc[unID]
			= (double)dode->collCommGetLongDouble();
		this->_globalWidomInstancesTloc[unID]
			= (double)dode->collCommGetLongDouble();

		/*
		 * construct the temperature profile
		 */
		double Tun = this->getGlobalCurrentTemperature();
		if((!this->_localRank) && (_universalDOFProfile[unID] >= 1000))
		{
		  double twoEkin = this->_universalKineticProfile[unID];

		  double vvNN = 0.0;
		  for(unsigned d = 0; d < 3; d++)
		     vvNN += _universalvProfile[d][unID] * _universalvProfile[d][unID];
		  double twoEkindir = _universalProfiledComponentMass * vvNN / _universalNProfile[unID];

		  Tun = (twoEkin - twoEkindir) / _universalDOFProfile[unID];
		}
		// domainDecomp->doBroadcast(&Tun); // no longer needed, since MPI_Reduce (branch) was replaced with MPI_Allreduce (trunk)
		this->_universalTProfile[unID] = Tun;
	}
	dode->collCommFinalize();
}

void Domain::outputProfile(const char* prefix, bool virialProfile)
{
	if(this->_localRank) return;

	string vzpryname(prefix);
	string Tpryname(prefix);
	string rhpryname(prefix);
        string upryname(prefix);
	string Vipryname(prefix);
	rhpryname += ".rhpry";
	vzpryname += ".vzpry";
	Tpryname += ".Tpry";
        upryname += ".upr";
	Vipryname += ".Vpry";
	ofstream rhpry(rhpryname.c_str());
	ofstream vzpry(vzpryname.c_str());
	ofstream Tpry(Tpryname.c_str());
	ofstream upry(upryname.c_str());
        ofstream* Vipry = NULL;
	if(virialProfile)
    {
       Vipry = new ofstream(Vipryname.c_str());
       Vipry->precision(5);
       *Vipry << "# y\tvn-vt\tpx\tpy\tpz\n# \n";
    }
	if (!(vzpry && Tpry && rhpry && upry))
	{
		return;
	}
	rhpry.precision(6);
	rhpry << "# y\trho\ttotal DOF\n# \n";
	vzpry.precision(5);
	vzpry << "# y\tvz\tv\n# \n";
	Tpry.precision(6);
        Tpry << "# coord\t2Ekin/#DOF\t2Ekin/3N (dir.)\tT\n# \n";
	upry.precision(5);
        upry << "# coord\t\tmu_conf(loc)  mu_conf(glob) \t\t mu_rho(loc)  mu_rho(glob) \t "
             << "mu_at(loc)  mu_at(glob) \t\t mu_T(loc)  mu_T(glob) \t mu_id(loc)  "
             << "mu_id(glob) \t\t mu_res(loc)  mu_res(glob) \t mu(loc)  mu(glob) \t\t #(loc)  "
             << "#(glob)\n";

	double layerVolume = this->_globalLength[0] * this->_globalLength[1] * this->_globalLength[2]
		/ this->_universalNProfileUnits[1];
	for(unsigned y = 0; y < this->_universalNProfileUnits[1]; y++)
	{
		double yval = (y + 0.5) / this->_universalInvProfileUnit[1];

		long double Ny = 0.0;
		long double DOFy = 0.0;
		long double twoEkiny = 0.0;
                long double twoEkindiry = 0.0;
		long double velocitysumy[3];
                long double widomSigExpy = 0.0;
                long double widomInstancesy = 0.0;
                long double widomSigExpyTloc = 0.0;
                long double widomInstancesyTloc = 0.0;
		long double Pd = 0.0;
		long double Px = 0.0;
		long double Py = 0.0;
		long double Pz = 0.0;
		for(unsigned d = 0; d < 3; d++) velocitysumy[d] = 0.0;
		for(unsigned x = 0; x < this->_universalNProfileUnits[0]; x++)
		{
			for(unsigned z = 0; z < this->_universalNProfileUnits[2]; z++)
			{
				unsigned unID = x * this->_universalNProfileUnits[1] * this->_universalNProfileUnits[2]
					+ y * this->_universalNProfileUnits[2] + z;
				Ny += this->_universalNProfile[unID];
				DOFy += this->_universalDOFProfile[unID];
				twoEkiny += this->_universalKineticProfile[unID];
				for(unsigned d = 0; d < 3; d++) velocitysumy[d] += this->_universalvProfile[d][unID];
                                widomSigExpy += this->_globalWidomProfile[unID];
                                widomInstancesy += this->_globalWidomInstances[unID];
                                widomSigExpyTloc += this->_globalWidomProfileTloc[unID];
                                widomInstancesyTloc += this->_globalWidomInstancesTloc[unID];

                                if(virialProfile)
                                {
				   // Pd += this->_universalPDProfile[unID]; // unnecessary redundancy
				   Px += this->_universalPXProfile[unID];
				   Py += this->_universalPYProfile[unID];
				   Pz += this->_universalPZProfile[unID];
                                }
			}
		}
		Pd = Py - 0.5*(Px + Pz);
                double rho_loc = Ny / (layerVolume * this->_globalAccumulatedDatasets);

		if(Ny >= 64.0)
		{
		   double vvdir = 0.0;
		   for(unsigned d = 0; d < 3; d++)
		   {
		      double vd = velocitysumy[d] / Ny;
		      vvdir += vd*vd;
		   }
                   twoEkindiry = Ny * _universalProfiledComponentMass * vvdir;

		   rhpry << yval << "\t" << rho_loc << "\t" << DOFy << "\n";
		   vzpry << yval << "\t" << (velocitysumy[2] / Ny) << "\t" << sqrt(vvdir) << "\n";
                   Tpry << yval << "\t" << (twoEkiny / DOFy) << "\t"
                        << (twoEkindiry / (3.0*Ny)) << "\t" << ((twoEkiny - twoEkindiry) / DOFy) << "\n";

		   if(virialProfile) *Vipry << yval << "\t" << Pd / (layerVolume * this->_globalAccumulatedDatasets) << "\t" << (_globalTemperatureMap[0]*Ny + Px) / (layerVolume * this->_globalAccumulatedDatasets) << "\t" << (_globalTemperatureMap[0]*Ny + Py) / (layerVolume * this->_globalAccumulatedDatasets) << "\t" << (_globalTemperatureMap[0]*Ny + Pz) / (layerVolume * this->_globalAccumulatedDatasets) << "\n";
                   if(widomInstancesy >= 100.0)
                   {
                      double mu_res_glob = -log(widomSigExpy / widomInstancesy);
                      double mu_conf_glob = mu_res_glob * _globalTemperatureMap[0];
                      double mu_id_glob = _globalTemperatureMap[0] * log(_globalDecisiveDensity);

                      double mu_T_glob = 0.0;
                      double mu_rho_loc = 0.0;
                      double mu_T_loc = 0.0;
                      double mu_res_loc = 0.0;
                      double mu_conf_loc = 0.0;
                      if(Ny >= 10.0)
                      {
                         mu_T_glob = 3.0*_globalTemperatureMap[0] * log(_universalLambda);
                         mu_res_loc = -log(widomSigExpyTloc / widomInstancesyTloc);

                         if(widomInstancesyTloc >= 100.0)
                         {
                            double Tloc = (twoEkiny - twoEkindiry) / DOFy;
                            // mu_id_loc = Tloc * log(rho_loc * _universalLambda*_universalLambda*_universalLambda);
                            mu_rho_loc = Tloc * log(rho_loc);
                            mu_T_loc = 3.0*Tloc * (log(_universalLambda) + 0.5*log(_globalTemperatureMap[0] / Tloc));
                            mu_conf_loc = Tloc * mu_res_loc;
                         }
                      }

                      upry << yval << " \t\t " << mu_conf_loc << "  " << mu_conf_glob << " \t\t "
                           << mu_rho_loc << "  " << mu_id_glob - mu_T_glob << " \t "
                           << mu_conf_loc + mu_rho_loc << "  " << mu_conf_glob + mu_id_glob - mu_T_glob << " \t\t "
                           << mu_T_loc << "  " << mu_T_glob << " \t "
                           << mu_rho_loc + mu_T_loc << "  " << mu_id_glob << " \t\t "
                           << mu_res_loc << "  " << mu_res_glob << " \t "
                           << mu_rho_loc + mu_T_loc + mu_conf_loc << "  " << mu_id_glob + mu_conf_glob << " \t\t "
                           << widomInstancesy << "  " << widomInstancesyTloc << "\n";
                   }
		}
		else
		{
			rhpry << yval << "\t0.000\t" << DOFy << "\n";
		}
	}

	rhpry.close();
	vzpry.close();
	Tpry.close();
	upry.close();
	if(virialProfile)
        {
           Vipry->close();
           delete Vipry;
        }
}


void Domain::outputKartesian2DProfile(const char* prefix, bool virialProfile)
{
	if(this->_localRank) return;

	string vzpryname(prefix);
	string Tprname(prefix);
	string Tpryname(prefix);
	string rhprname(prefix);
	string rhpryname(prefix);
    string upryname(prefix);
	string Viprname(prefix);
	string Vipryname(prefix);
	rhprname += "_Kart2D.rhpr";
	rhpryname += ".rhpry";
	vzpryname += ".vzpry";
	Tprname += "_Kart2D.Tpr";
	Tpryname += ".Tpry";
    upryname += ".upr";
	Viprname += "_Kart2D.Vpr";
	Vipryname += ".Vpry";
	ofstream rhpry(rhpryname.c_str());
	ofstream rhpr(rhprname.c_str());
	ofstream vzpry(vzpryname.c_str());
	ofstream Tpry(Tpryname.c_str());
	ofstream Tpr(Tprname.c_str());
	ofstream upry(upryname.c_str());
    ofstream Vipry(Vipryname.c_str());
	ofstream Vipr(Viprname.c_str());
	if(virialProfile)
    {
       //Vipry = new ofstream(Vipryname.c_str());
       Vipry.precision(5);
	   //Vipr = new ofstream(Viprname.c_str());
       Vipr.precision(5);
       Vipry << "# y\tvn-vt\tpx\tpy\tpz\n# \n";
    }
	if (!(vzpry && Tpry && rhpry && upry && rhpr && Tpr))
	{
		return;
	}
	rhpry.precision(6);
	rhpr.precision(6);
	rhpry << "# y\trho\ttotal DOF\n# \n";
	vzpry.precision(5);
	vzpry << "# y\tvz\tv\n# \n";
	Tpry.precision(6);
	Tpr.precision(6);
    Tpry << "# coord\t2Ekin/#DOF\t2Ekin/3N (dir.)\tT\n# \n";
	upry.precision(5);
        upry << "# coord\t\tmu_conf(loc)  mu_conf(glob) \t\t mu_rho(loc)  mu_rho(glob) \t "
             << "mu_at(loc)  mu_at(glob) \t\t mu_T(loc)  mu_T(glob) \t mu_id(loc)  "
             << "mu_id(glob) \t\t mu_res(loc)  mu_res(glob) \t mu(loc)  mu(glob) \t\t #(loc)  "
             << "#(glob)\n";

	double layerVolume = this->_globalLength[0] * this->_globalLength[1] * this->_globalLength[2]
		/ this->_universalNProfileUnits[1];
	double segmentVolume = this->_globalLength[0] * this->_globalLength[1] * this->_globalLength[2]
		/ (this->_universalNProfileUnits[1]*this->_universalNProfileUnits[2]);

	rhpr << "//Segment volume: " << segmentVolume << "\n//Accumulated data sets: " << _globalAccumulatedDatasets << "\n//Local profile of the number density. Output file generated by the \"outputCylProfile\" method, located in Domain.cpp. \n";
	rhpr << "//local density profile: Each matrix corresponds to a single value of \"phi_i\", measured in [rad]\n";
	rhpr << "//one single matrix of the local number density rho'(phi_i;r_i',h_i') \n//       | r_i'\n// ---------------------\n//   h_i'| rho'(r_i',h_i')\n//       | \n";
	rhpr << "//  T' \t sigma_ii' \t eps_ii' \t yOffset \t DELTA_phi \t DELTA_r2' \t DELTA_h' \t quantities in atomic units are denoted by an apostrophe '\n";
	rhpr << this->_universalTargetTemperature[_simulation.getEnsemble()->getComponents()->size()]<<"\t"<<_simulation.getEnsemble()->getComponent(0)->getSigma(0)<<"\t"<<_simulation.getEnsemble()->getComponent(0)->getEps(0)<<"\t"; //changed by Michaela Heier
	rhpr << _yOff << "\t" << 1/this->_universalInvProfileUnit[0] << "\t" << 1/this->_universalInvProfileUnit[1] << "\t" << 1/this->_universalInvProfileUnit[2]<< "\n";
	rhpr << "0 \t";

	Tpr << "//Local temperature profile generated by the \"outputCylProfile\" method.\n";
	Tpr << "//Temperature expressed by 2Ekin/#DOF\n";
	Tpr << "//one single matrix of the local temperature T(r_i,h_i) \n//      | r_i\n//---------------------\n//  h_i| T(r_i,h_i)\n//      | \n";
	Tpr << "// T \t sigma_ii \t eps_ii \t yOffset \t DELTA_phi \t DELTA_r2 \t DELTA_h \n";
	Tpr << this->_universalTargetTemperature[_simulation.getEnsemble()->getComponents()->size()] <<"\t"<<_simulation.getEnsemble()->getComponent(0)->getSigma(0)<<"\t"<<_simulation.getEnsemble()->getComponent(0)->getEps(0)<<"\t"; //changed by Michaela Heier
	Tpr << _yOff << "\t" << 1/this->_universalInvProfileUnit[0] << "\t" << 1/this->_universalInvProfileUnit[1] << "\t" << 1/this->_universalInvProfileUnit[2]<< "\n";
	Tpr <<"0 \t";
	
	if(virialProfile)
    {
	Vipr << "//Local profile of the pressure. Output file generated by the \"outputCylProfile\" method, located in Domain.cpp. \n";
	Vipr << "//local pressure profile: Each matrix corresponds to a single value of \"phi_i\", measured in [rad]\n";
	Vipr << "//one single matrix of the local number pressure p'(phi_i;r_i',h_i') \n//      | r_i'\n//---------------------\n//  h_i'| p'(r_i',h_i')\n//      | \n";
	Vipr << "// T' \t sigma_ii' \t eps_ii' \t yOffset \t DELTA_phi \t DELTA_r2' \t DELTA_h' \t quantities in atomic units are denoted by an apostrophe '\n";
	Vipr << this->_universalTargetTemperature[_simulation.getEnsemble()->getComponents()->size()]<<"\t"<<_simulation.getEnsemble()->getComponent(0)->getSigma(0)<<"\t"<<_simulation.getEnsemble()->getComponent(0)->getEps(0)<<"\t"; 
	Vipr << _yOff << "\t" << 1/this->_universalInvProfileUnit[0] << "\t" << 1/this->_universalInvProfileUnit[1] << "\t" << 1/this->_universalInvProfileUnit[2]<< "\n";
	Vipr << "0 \t";
	}
//for(unsigned x = 0; x < this->_universalNProfileUnits[0]; x++)
	//{
		//rhpry <<"> "<< (z+0.5)/this->_universalInvProfileUnit[2] <<"\n 0 \t";
	   	for(unsigned z = 0; z< this->_universalNProfileUnits[2]; z++){
	   		rhpr << (z+0.5)/this->_universalInvProfileUnit[2] <<"  \t"; // Eintragen der radialen Koordinaten r_i in Header
			Tpr << (z+0.5)/this->_universalInvProfileUnit[2] <<"\t";; // Eintragen der radialen Koordinaten r_i
			if(virialProfile) Vipr << (z+0.5)/this->_universalInvProfileUnit[2] <<"  \t"; // Eintragen der radialen Koordinaten r_i in Header 

	   	}
	   	rhpr << "\n";
		Tpr << "\n";
		if(virialProfile) Vipr << "\n";

	   	for(unsigned y = 0; y < this->_universalNProfileUnits[1]; y++)
	   		{
			//long double Ny = 0.0;
			//long double DOFy = 0.0;
			//long double twoEkiny = 0.0;
			long double twoEkindiry = 0.0;
			//long double velocitysumy[3];
			
			double hval = (y + 0.5) / this->_universalInvProfileUnit[1];
	   	    Tpr << hval<< "  \t";
			rhpr << hval<< "  \t";
			if(virialProfile) Vipr << hval<< "  \t";

	   		for(unsigned z = 0; z < this->_universalNProfileUnits[2]; z++)
	   	    {
				for(unsigned x = 0; x < this->_universalNProfileUnits[0];x++)
				{
	   	        	unsigned unID = x * this->_universalNProfileUnits[0] * this->_universalNProfileUnits[2] + y * this->_universalNProfileUnits[1] + z;
																	
	   				double rho_loc = this->_universalNProfile[unID] / (segmentVolume * this->_globalAccumulatedDatasets);
	   				rhpr << rho_loc << "\t";
					//Ny += this->_universalNProfile[unID];
					//DOFy += this->_universalDOFProfile[unID];
					//twoEkiny += this->_universalKineticProfile[unID];
					//for(unsigned d = 0; d < 3; d++) velocitysumy[d] += this->_universalvProfile[d][unID];
					if(virialProfile)
					{
						double virial =(_globalTemperatureMap[0]*this->_universalNProfile[unID]+this->_universalPYProfile[unID]) / (segmentVolume * this->_globalAccumulatedDatasets);
	   	        		Vipr << virial << "\t";
					}
				
				if(this->_universalDOFProfile[unID] == 0.0){
					Tpr << 0 << "\t";
				}
				else{
					//if(this->_universalNProfile[unID] >= 64.0)
					//{
						double vvdir = 0.0;
						for(unsigned d = 0; d < 3; d++)
						{
							double vd = this->_universalvProfile[d][unID] / (this->_universalNProfile[unID]);
							vvdir += vd*vd;
						}

						twoEkindiry = this->_universalNProfile[unID] * _universalProfiledComponentMass * vvdir;
						Tpr << ((this->_universalKineticProfile[unID] - twoEkindiry) / (this->_universalDOFProfile[unID])) << "\t";

					//}
						
				}
				}
			}
	   	    rhpr << "\n";
			Tpr << "\n";
			if(virialProfile) Vipr << "\n";
	   	}
//}
		       
	for(unsigned y = 0; y < this->_universalNProfileUnits[1]; y++)
	{
		double yval = (y + 0.5) / this->_universalInvProfileUnit[1];

		long double Ny = 0.0;
		long double DOFy = 0.0;
		long double twoEkiny = 0.0;
                long double twoEkindiry = 0.0;
		long double velocitysumy[3];
                long double widomSigExpy = 0.0;
                long double widomInstancesy = 0.0;
                long double widomSigExpyTloc = 0.0;
                long double widomInstancesyTloc = 0.0;
		long double Pd = 0.0;
		long double Px = 0.0;
		long double Py = 0.0;
		long double Pz = 0.0;
		for(unsigned d = 0; d < 3; d++) velocitysumy[d] = 0.0;
		for(unsigned x = 0; x < this->_universalNProfileUnits[0]; x++)
		{
			for(unsigned z = 0; z < this->_universalNProfileUnits[2]; z++)
			{
				unsigned unID = x * this->_universalNProfileUnits[1] * this->_universalNProfileUnits[2]
					+ y * this->_universalNProfileUnits[2] + z;
				Ny += this->_universalNProfile[unID];
				DOFy += this->_universalDOFProfile[unID];
				twoEkiny += this->_universalKineticProfile[unID];
				for(unsigned d = 0; d < 3; d++) velocitysumy[d] += this->_universalvProfile[d][unID];
                                widomSigExpy += this->_globalWidomProfile[unID];
                                widomInstancesy += this->_globalWidomInstances[unID];
                                widomSigExpyTloc += this->_globalWidomProfileTloc[unID];
                                widomInstancesyTloc += this->_globalWidomInstancesTloc[unID];
                                
                                if(virialProfile)
                                {
				  // Pd += this->_universalPDProfile[unID]; // unnecessary redundancy
				   Px += this->_universalPXProfile[unID];
				   Py += this->_universalPYProfile[unID];
				   Pz += this->_universalPZProfile[unID];
                               }
			}
		}

		Pd = Py - 0.5*(Px + Pz);
        double rho_loc = Ny / (layerVolume * this->_globalAccumulatedDatasets);

		if(Ny >= 64.0)
		{
		   double vvdir = 0.0;
		   for(unsigned d = 0; d < 3; d++)
		   {
		      double vd = velocitysumy[d] / Ny;
		      vvdir += vd*vd;
		   }
                   twoEkindiry = Ny * _universalProfiledComponentMass * vvdir;

		   rhpry << yval << "\t" << rho_loc << "\t" << DOFy << "\n";
		   vzpry << yval << "\t" << (velocitysumy[2] / Ny) << "\t" << sqrt(vvdir) << "\n";
                   Tpry << yval << "\t" << (twoEkiny / DOFy) << "\t"
                        << (twoEkindiry / (3.0*Ny)) << "\t" << ((twoEkiny - twoEkindiry) / DOFy) << "\n";

		   if(virialProfile) Vipry << yval << "\t" << Pd / (layerVolume * this->_globalAccumulatedDatasets) << "\t" << (_globalTemperatureMap[0]*Ny + Px) / (layerVolume * this->_globalAccumulatedDatasets) << "\t" << (_globalTemperatureMap[0]*Ny + Py) / (layerVolume * this->_globalAccumulatedDatasets) << "\t" << (_globalTemperatureMap[0]*Ny + Pz) / (layerVolume * this->_globalAccumulatedDatasets) << "\n";
                   if(widomInstancesy >= 100.0)
                   {
                      double mu_res_glob = -log(widomSigExpy / widomInstancesy);
                      double mu_conf_glob = mu_res_glob * _globalTemperatureMap[0];
                      double mu_id_glob = _globalTemperatureMap[0] * log(_globalDecisiveDensity);

                      double mu_T_glob = 0.0;
                      double mu_rho_loc = 0.0;
                      double mu_T_loc = 0.0;
                      double mu_res_loc = 0.0;
                      double mu_conf_loc = 0.0;
                      if(Ny >= 10.0)
                      {
                         mu_T_glob = 3.0*_globalTemperatureMap[0] * log(_universalLambda);
                         mu_res_loc = -log(widomSigExpyTloc / widomInstancesyTloc);              

                         if(widomInstancesyTloc >= 100.0)
                         {
                            double Tloc = (twoEkiny - twoEkindiry) / DOFy;
                            // mu_id_loc = Tloc * log(rho_loc * _universalLambda*_universalLambda*_universalLambda);
                            mu_rho_loc = Tloc * log(rho_loc);
                            mu_T_loc = 3.0*Tloc * (log(_universalLambda) + 0.5*log(_globalTemperatureMap[0] / Tloc));
                            mu_conf_loc = Tloc * mu_res_loc;
                         }
                      }
               
                      upry << yval << " \t\t " << mu_conf_loc << "  " << mu_conf_glob << " \t\t "
                           << mu_rho_loc << "  " << mu_id_glob - mu_T_glob << " \t "
                           << mu_conf_loc + mu_rho_loc << "  " << mu_conf_glob + mu_id_glob - mu_T_glob << " \t\t "
                           << mu_T_loc << "  " << mu_T_glob << " \t "
                           << mu_rho_loc + mu_T_loc << "  " << mu_id_glob << " \t\t "
                           << mu_res_loc << "  " << mu_res_glob << " \t "
                           << mu_rho_loc + mu_T_loc + mu_conf_loc << "  " << mu_id_glob + mu_conf_glob << " \t\t "
                           << widomInstancesy << "  " << widomInstancesyTloc << "\n";
                   }
		}
		else
		{
			rhpry << yval << "\t0.000\t" << DOFy << "\n";
		}
	}

	rhpry.close();
	rhpr.close();
	vzpry.close();
	Tpry.close();
	Tpr.close();
	upry.close();
	if(virialProfile)
        {
           Vipry.close();
		   Vipr.close();
           //delete Vipry;
        }
}

void Domain::outputDropMove(double _universalRealignmentMotionX,double _universalRealignmentMotionZ,string prefix, unsigned long timestep,string diff_number)
{
	prefix = prefix.substr(0, prefix.length()-4);	
        string outDMname(prefix);
        outDMname += diff_number + ".dat";
        ofstream DM(outDMname.c_str(),ios::out|ios::app);

	DM.precision(6);

	_universalRealignmentMotionX = - _universalRealignmentMotionX;
	_universalRealignmentMotionZ = - _universalRealignmentMotionZ;

	DM << timestep << "\t" << _universalRealignmentMotionX << "\t" << _universalRealignmentMotionZ << "\n";
}

void Domain::resetProfile(bool virialProfile)
{
	unsigned unIDs = this->_universalNProfileUnits[0] * this->_universalNProfileUnits[1]
		* this->_universalNProfileUnits[2];
	for(unsigned unID = 0; unID < unIDs; unID++)
	{
		this->_localNProfile[unID] = 0.0;
		this->_universalNProfile[unID] = 0.0;
		for(int d=0; d<3; d++)
		{
			this->_localvProfile[d][unID] = 0.0;
			this->_universalvProfile[d][unID] = 0.0;
//			this->_localV2Profile[d][unID] = 0.0;
//			this->_universalV2Profile[d][unID] = 0.0;
		}
		this->_localDOFProfile[unID] = 0.0;
		this->_universalDOFProfile[unID] = 0.0;
		this->_localKineticProfile[unID] = 0.0;
		this->_universalKineticProfile[unID] = 0.0;

                if(virialProfile)
                {
		   // this->_localPDProfile[unID] = 0.0;
		   // this->_universalPDProfile[unID] = 0.0;
		   this->_localPXProfile[unID] = 0.0;
		   this->_universalPXProfile[unID] = 0.0;
		   this->_localPYProfile[unID] = 0.0;
		   this->_universalPYProfile[unID] = 0.0;
		   this->_localPZProfile[unID] = 0.0;
		   this->_universalPZProfile[unID] = 0.0;
                }

		this->_localWidomProfile[unID] = 0.0;
		this->_globalWidomProfile[unID] = 0.0;
		this->_localWidomInstances[unID] = 0.0;
		this->_globalWidomInstances[unID] = 0.0;
		this->_localWidomProfileTloc[unID] = 0.0;
		this->_globalWidomProfileTloc[unID] = 0.0;
		this->_localWidomInstancesTloc[unID] = 0.0;
		this->_globalWidomInstancesTloc[unID] = 0.0;
	}
	this->_globalAccumulatedDatasets = 0;
}


void Domain::Nadd(unsigned cid, int N, int localN)
{
	Ensemble* ensemble = _simulation.getEnsemble();
	Component* component = ensemble->getComponent(cid);
	component->incNumMolecules(N);
	unsigned int rotationDegreesOfFreeedom = component->getRotationalDegreesOfFreedom();

	this->_globalNumMolecules += N;
	if( (this->_componentwiseThermostat)
			&& (this->_componentToThermostatIdMap[cid] > 0) )
	{
		int thid = this->_componentToThermostatIdMap[cid];
		this->_localThermostatN[thid] += localN;
		this->_universalThermostatN[thid] += N;
		this->_localRotationalDOF[thid] += localN * rotationDegreesOfFreeedom;
		this->_universalRotationalDOF[thid] += N * rotationDegreesOfFreeedom;
	}
	this->_localThermostatN[0] += localN;
	this->_universalThermostatN[0] += N;
	this->_localRotationalDOF[0] += localN * rotationDegreesOfFreeedom;
	this->_universalRotationalDOF[0] += N * rotationDegreesOfFreeedom;
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

	/* FIXME: Substantial change in program behavior! */
	if(thermostat == 0) {
		global_log->warning() << "Disabling the component wise thermostat!" << endl;
		disableComponentwiseThermostat();
	}
	if(thermostat >= 1) {
		if( ! _componentwiseThermostat ) {
			/* FIXME: Substantial change in program behavior! */
			global_log->warning() << "Enabling the component wise thermostat!" << endl;
			_componentwiseThermostat = true;
			_universalTargetTemperature.erase(0);
			_universalUndirectedThermostat.erase(0);
			for(int d=0; d < 3; d++) this->_universalThermostatDirectedVelocity[d].erase(0);
			vector<Component>* components = _simulation.getEnsemble()->getComponents();
			for( vector<Component>::iterator tc = components->begin(); tc != components->end(); tc ++ ) {
				if(!(this->_componentToThermostatIdMap[ tc->ID() ] > 0)) {
					this->_componentToThermostatIdMap[ tc->ID() ] = -1;
				}
			}
		}
	}
}

void Domain::enableComponentwiseThermostat()
{
	if(this->_componentwiseThermostat) return;

	this->_componentwiseThermostat = true;
	this->_universalTargetTemperature.erase(0);
	vector<Component>* components = _simulation.getEnsemble()->getComponents();
	for( vector<Component>::iterator tc = components->begin(); tc != components->end(); tc ++ ) {
		if(!(this->_componentToThermostatIdMap[ tc->ID() ] > 0)) {
			this->_componentToThermostatIdMap[ tc->ID() ] = -1;
		}
	}
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
	double globalUpot = getGlobalUpot();
	this->_globalSigmaU += globalUpot;
	this->_globalSigmaUU += globalUpot * globalUpot;
}

double Domain::cv()
{
	if((_localRank != 0) || (_globalUSteps == 0)) return 0.0;

	double id = 1.5 + 0.5*_universalRotationalDOF[0]/_globalNumMolecules;
	double conf = (_globalSigmaUU - _globalSigmaU*_globalSigmaU/_globalUSteps)
		/ (_globalUSteps * _globalNumMolecules * _globalTemperatureMap[0] * _globalTemperatureMap[0]);

	return id + conf;
}

//! methods implemented by Stefan Becker <stefan.becker@mv.uni-kl.de>
// the following two methods are used by the MmspdWriter (writing the output file in a format used by MegaMol)
double Domain::getSigma(unsigned cid, unsigned nthSigma){
  return _simulation.getEnsemble()->getComponent(cid)->getSigma(nthSigma);
}
unsigned Domain::getNumberOfComponents(){
  return _simulation.getEnsemble()->getComponents()->size();
}

void Domain::submitDU(unsigned /*cid*/, double DU, double* r)
{
   unsigned xun, yun, zun;
   xun = (unsigned)floor(r[0] * this->_universalInvProfileUnit[0]);
   yun = (unsigned)floor(r[1] * this->_universalInvProfileUnit[1]);
   zun = (unsigned)floor(r[2] * this->_universalInvProfileUnit[2]);
   int unID = xun * this->_universalNProfileUnits[1] * this->_universalNProfileUnits[2]
                  + yun * this->_universalNProfileUnits[2] + zun;
   if(unID < 0) return;

   _localWidomProfile[unID] += exp(-DU / _globalTemperatureMap[0]);
   _localWidomInstances[unID] += 1.0;

   double Tloc = _universalTProfile[unID];
   if(Tloc != 0.0)
   {
      _localWidomProfileTloc[unID] += exp(-DU/Tloc);
      _localWidomInstancesTloc[unID] += 1.0;
   }
}

void Domain::resetGamma(){
	for (unsigned i=0; i<_simulation.getEnsemble()->getComponents()->size(); i++){
		_Gamma[i]=0;
	}
}

double Domain::getGamma(unsigned id){
	return (_Gamma[id]/(2*_globalLength[0]*_globalLength[2]));
}

void Domain::calculateGamma(ParticleContainer* _particleContainer, DomainDecompBase* _domainDecomposition){
	unsigned numComp = _simulation.getEnsemble()->getComponents()->size();
	double _localGamma[numComp];
	for (unsigned i=0; i<numComp; i++){
		_localGamma[i]=0;
	}
	for(ParticleIterator tempMol = _particleContainer->iteratorBegin(); tempMol != _particleContainer->iteratorEnd(); ++tempMol){
		unsigned cid=tempMol->componentid();
		_localGamma[cid]+=tempMol->Vi(1)-0.5*(tempMol->Vi(0)+tempMol->Vi(2));
	//	cout << _localGamma[cid] << "\t" << cid << endl;
	}
	_domainDecomposition->collCommInit(numComp);
	for (unsigned i=0; i<numComp; i++){
		_domainDecomposition->collCommAppendDouble(_localGamma[i]);
	}
	_domainDecomposition->collCommAllreduceSum();
	for (unsigned i=0; i<numComp; i++){
		_localGamma[i] = _domainDecomposition->collCommGetDouble();
	}
	_domainDecomposition->collCommFinalize();
	for (unsigned i=0; i<numComp; i++){
		_Gamma[i]+=_localGamma[i];
	}
}
void Domain::considerComponentForYShift(unsigned cidMin, unsigned cidMax){
    for (unsigned i = 0; i <= cidMax; i++) _componentForYShift[i] = false;
    for (unsigned i = 0; i <= cidMax-cidMin; i++) _componentForYShift[cidMin+i] = true;
#ifndef NDEBUG
    global_log->info() << "Components used for determination of shift in y-direction:\n";
    for(unsigned i = 0; i <= cidMax; i++)
      global_log->info() << "component id "<< i << "\t value: " << _componentForYShift[i] << "\n";
#endif
}


void Domain::sYOffset(double in_yOff){
  _yOff = in_yOff;
}

// author: Stefan Becker, method to determine the integer value of unID, extra method in order to not inflating the method "recordProfile()"
// method only matched for the case of a cylindrical profile (i.e. sessile drop)
long int Domain::unID(double qx, double qy, double qz){
		int xun,yun,zun;// (xun,yun,zun): bin number in a special direction, e.g. yun==5 corresponds to the 5th bin in the radial direction,
		long int unID;	// as usual
		double xc,yc,zc; // distance of a particle with respect to the origin of the cylindrical coordinate system

		unID = -1; // initialization, causes an error message, if unID is not calculated in this method but used in record profile

	    xc = qx - this->_universalCentre[0];
	    yc = qy - this->_universalCentre[1];
	    zc = qz - this->_universalCentre[2];

	    // transformation in polar coordinates
	    double R2 = xc*xc + zc*zc;
	    double phi = asin(zc/sqrt(R2)) + ((xc>=0.0) ? 0:M_PI);
	    if(phi<0.0) {phi = phi + 2.0*M_PI;}

	    xun = (int)floor(phi * this->_universalInvProfileUnit[0]);   // bin no. in phi-direction
	    yun = (int)floor(R2 *  this->_universalInvProfileUnit[1]);   // bin no. in R-direction
	    zun = (int)floor(yc *  this->_universalInvProfileUnit[2]);   // bin no. in H-direction

	    if((xun >= 0) && (yun >= 0) && (zun >= 0) &&
	          (xun < (int)_universalNProfileUnits[0]) && (yun < (int)_universalNProfileUnits[1]) && (zun < (int)_universalNProfileUnits[2]))
	       {
	          unID = xun * this->_universalNProfileUnits[1] * this->_universalNProfileUnits[2]
	               + yun * this->_universalNProfileUnits[2] + zun;
	       }
	       else
	       {
	          global_log->error() << "Severe error!! Invalid profile unit (" << xun << " / " << yun << " / " << zun << ").\n\n";
	          global_log->error() << "Coordinates (" << qx << " / " << qy << " / " << qz << ").\n";
		  global_log->error() << "unID = " << unID << "\n";
	          //Simulation::exit(707);
	       }
	       return unID;
}

// author: Stefan Becker, counterpart of method setupProfile(xun,yun,zun) => to be used for profile in cylindrical coordinates
// the radial component is located in the x,z-plane, the axial component is parallel to the y-direction.
// In the radial direction the spacing is linear in R^2 in order to obtain equal areas with increasing radius.
// @TODO: in the input file the token "profile" immedeatly causes a call of "setupProfile(xun,yun,zun)". This in turn sets up a cubic grid
// by default.
// Proposal: an enquiry that checks wheter a cubic, or cylindrical, sperical, etc. grid has to be set up. This first requires a check of
// the input file (what additional token is set that controls the kind of grid beeing used). => code more modular
void Domain::sesDrop(){
	this->_universalCylindricalGeometry = true;

	// R^2_max: (squared) maximum radius up to which the profile is recorded
	// (otherwise erroneous densities recorded at the boundary of the RECTANGULAR simulation box)
	double minXZ = this->_globalLength[0];
	if(this->_globalLength[2]<minXZ){
		minXZ = this->_globalLength[2];
	}
	this->_universalR2max = 0.24*minXZ*minXZ;

	// origin of the cylindrical coordinate system
	this->_universalCentre[0] = 0.5*this->_globalLength[0];
	this->_universalCentre[1] = 0;
	this->_universalCentre[2] = 0.5*this->_globalLength[2];

	_universalInvProfileUnit[0] = this->_universalNProfileUnits[0]/(2*M_PI);                   // delta_phi
	_universalInvProfileUnit[1] = this->_universalNProfileUnits[1]/(this->_universalR2max);  // delta_R^2
	_universalInvProfileUnit[2] = this->_universalNProfileUnits[2]/(this->_globalLength[1]); // delta_H
	global_log->info() << "\nInv Profile unit for sessile drop: (phi,R2,H) = (" << _universalInvProfileUnit[0] <<", " <<_universalInvProfileUnit[1]<<", "<<_universalInvProfileUnit[2]<<") \n";

	this->resetProfile(true);
}

/*
 *By Stefan Becker: 
 *realign tool, borrowed from Martin Horsch: it effects that the center of mass of all the particles (with respect to the (x,z)-plane)
 *corresponds with the box center (in the (x,z)plane). With respect to the y-direction, the shift is carried out so 
 *that the wall remains at the initial position.
 *Basically the current centre of mass of all the particles is determined and then all the particles are moved towards x==0.5*Lx and z==0*Lz
 *The method is frequently applied when the specified span of time (i.e. timesteps) is reached. 
 *Part of this tool in Domain.cpp: (i) determineShift() => determines the vector by which the particles are shifted
 *                                 (ii) realign() => carries out the actual realignment / shift
 * when this method is called, the halo should NOT be present
 */
void Domain::determineXZShift( DomainDecompBase* domainDecomp, ParticleContainer* molCont,
			     double fraction)
{
   double localBalance[3]; // localBalance[1] not considered, employed only for not confusing: index 0 => x-direction, index 1 => y-direction, index 2 => z-direction
   for(unsigned d = 0; d < 3; d ++) localBalance[d] = 0.0; // initialising the array by zeros
   double localMass = 0.0;
   int cid;

   for(ParticleIterator tm = molCont->iteratorBegin(); tm != molCont->iteratorEnd(); ++tm)
   {
      cid = tm->componentid();
      if(_universalProfiledComponents[cid])
      {
	double tmass = tm->gMass();
	localMass += tmass;
	for(unsigned d = 0; d < 3; d=d+2){ // counting d=d+2 causes the value d==1 (i.e. y-direction) to be skipped, this is performed by the method determineYShift()
	  localBalance[d] += tm->r(d) * tmass;
	}
      }
   }

   for(unsigned d = 0; d < 3; d++) _globalRealignmentBalance[d] = localBalance[d];
   _globalRealignmentMass[0] = localMass;
   // determining the global values by the use of collectiveCommunication
   domainDecomp->collCommInit(4);
   for(unsigned d = 0; d < 3; d++)  domainDecomp->collCommAppendDouble(_globalRealignmentBalance[d]);
   domainDecomp->collCommAppendDouble(_globalRealignmentMass[0]);
   domainDecomp->collCommAllreduceSum();
   for(unsigned d = 0; d < 3; d++)  _globalRealignmentBalance[d] = domainDecomp->collCommGetDouble();
   _globalRealignmentMass[0] = domainDecomp->collCommGetDouble();
   domainDecomp->collCommFinalize();

   for(unsigned short d = 0; d < 3; d=d+2){
      _universalRealignmentMotion[d]
         = -fraction*((_globalRealignmentBalance[d] / _globalRealignmentMass[0]) - 0.5*_globalLength[d]);
   }
}

void Domain::determineYShift( DomainDecompBase* domainDecomp, ParticleContainer* molCont,
			     double fraction){
			        // keep in mind: variables declared as static are initialized only ONCE!
   static unsigned timesCalled = 0;
   static double initialCentreOfMassY;
   double localBalance = 0.0;
   double localMass = 0.0;
   int cid;
   for(ParticleIterator tm = molCont->iteratorBegin(); tm != molCont->iteratorEnd(); ++tm)
   {
     cid = tm->componentid();
     if(_componentForYShift[cid])
     {
	double tmass = tm->gMass();

	localMass += tmass;
	localBalance += (tm->r(1) - _globalLength[1]*floor(2.0* tm->r(1)/ _globalLength[1]))*tmass; // PBC

	 }
   }
   _globalRealignmentBalance[1] = localBalance;
   _globalRealignmentMass[1] = localMass;
   // determining the global values by the use of collectiveCommunication
   domainDecomp->collCommInit(2);
   domainDecomp->collCommAppendDouble(_globalRealignmentBalance[1]);
   domainDecomp->collCommAppendDouble(_globalRealignmentMass[1]);
   domainDecomp->collCommAllreduceSum();
   _globalRealignmentBalance[1] = domainDecomp->collCommGetDouble();
   _globalRealignmentMass[1] = domainDecomp->collCommGetDouble();
   domainDecomp->collCommFinalize();

   /*
   The centre of mass of the wall (solely the wall!) is determined once at the beginning of the simulation.
   Then the realignment with respect to the y-direction is always carried out so that the wall remains at the initial posistion
   (initial y-coordinate)
   */


   if (!timesCalled){
     initialCentreOfMassY = _globalRealignmentBalance[1] / _globalRealignmentMass[1];
   }
   timesCalled++;
   //global_log->info() <<"initialCentreOfMassY = " << initialCentreOfMassY << endl;
   // the realignment motion in y-direction so that the wall is always at the bottom of the simulation box
   _universalRealignmentMotion[1] = -fraction*((_globalRealignmentBalance[1] / _globalRealignmentMass[1]) - initialCentreOfMassY);

   }

// no y-shift will be determined. edited by Michaela Heier
void Domain::noYShift( DomainDecompBase* domainDecomp, ParticleContainer* molCont,
			     double /*fraction*/){
			        // keep in mind: variables declared as static are initialized only ONCE!
   //static unsigned timesCalled = 0;
   //static double initialCentreOfMassY;
   double localBalance = 0.0;
   double localMass = 0.0;
   for(ParticleIterator tm = molCont->iteratorBegin(); tm != molCont->iteratorEnd(); ++tm)
   {

	double tmass = tm->gMass();

	localMass += tmass;
	localBalance += (tm->r(1) - _globalLength[1]*floor(2.0* tm->r(1)/ _globalLength[1]))*tmass; // PBC


   }
   _globalRealignmentBalance[1] = localBalance;
   _globalRealignmentMass[1] = localMass;
   // determining the global values by the use of collectiveCommunication
   domainDecomp->collCommInit(2);
   domainDecomp->collCommAppendDouble(_globalRealignmentBalance[1]);
   domainDecomp->collCommAppendDouble(_globalRealignmentMass[1]);
   domainDecomp->collCommAllreduceSum();
   _globalRealignmentBalance[1] = domainDecomp->collCommGetDouble();
   _globalRealignmentMass[1] = domainDecomp->collCommGetDouble();
   domainDecomp->collCommFinalize();

   /*
   The centre of mass of the wall (solely the wall!) is determined once at the beginning of the simulation.
   Then the realignment with respect to the y-direction is always carried out so that the wall remains at the initial posistion
   (initial y-coordinate)
   */

   /*
   if (!timesCalled){
     initialCentreOfMassY = _globalRealignmentBalance[1] / _globalRealignmentMass[1];
   }

   timesCalled++;
   */

   //global_log->info() <<"initialCentreOfMassY = " << initialCentreOfMassY << endl;
   // the realignment motion in y-direction so that the wall is always at the bottom of the simulation box
   //_universalRealignmentMotion[1] = -fraction*((_globalRealignmentBalance[1] / _globalRealignmentMass[1]) - initialCentreOfMassY);

   // the realignment motion is neglected
   _universalRealignmentMotion[1] = 0;

   }


void Domain::determineShift( DomainDecompBase* domainDecomp, ParticleContainer* molCont,
                             double fraction)
{
   double localBalance[3];
   for(unsigned d = 0; d < 3; d ++) localBalance[d] = 0.0; // initialising the array by zeros
   double localMass = 0.0;
   int cid;

   for(ParticleIterator tm = molCont->iteratorBegin(); tm != molCont->iteratorEnd(); ++tm)
   {
      cid = tm->componentid();
      if(_universalProfiledComponents[cid])
      {
        double tmass = tm->gMass();
        localMass += tmass;
        for(unsigned d = 0; d < 3; d++){
          localBalance[d] += tm->r(d) * tmass;
        }
      }
   }

   for(unsigned d = 0; d < 3; d++) _globalRealignmentBalance[d] = localBalance[d];
   _globalRealignmentMass[0] = localMass;
   // determining the global values by the use of collectiveCommunication
   domainDecomp->collCommInit(4);
   for(unsigned d = 0; d < 3; d++)  domainDecomp->collCommAppendDouble(_globalRealignmentBalance[d]);
   domainDecomp->collCommAppendDouble(_globalRealignmentMass[0]);
   domainDecomp->collCommAllreduceSum();
   for(unsigned d = 0; d < 3; d++)  _globalRealignmentBalance[d] = domainDecomp->collCommGetDouble();
   _globalRealignmentMass[0] = domainDecomp->collCommGetDouble();
   domainDecomp->collCommFinalize();

   for(unsigned short d = 0; d < 3; d++){
      _universalRealignmentMotion[d]
         = -fraction*((_globalRealignmentBalance[d] / _globalRealignmentMass[0]) - 0.5*_globalLength[d]);
   }
  // "quick 'n dirty hack" in order to avoid trouble with the info-line of the realign()-method, this info line uses _globalRealignmentMass[1] which does
  // occur in this method (here it is superfluos)
   _globalRealignmentMass[1] = _globalRealignmentMass[0];
}

/*
 * By Stefan Becker, see above comment on determineShift().
 * this method REQUIRES the presence of the halo
 */
void Domain::realign(
   ParticleContainer* molCont, bool _dropMove,string prefix, unsigned long _simstep, string diff_number) {
   if(!this->_localRank)
   {
#ifndef NDEBUG
      cout << "Centre of mass: (" << _globalRealignmentBalance[0]/_globalRealignmentMass[0] << " / " << _globalRealignmentBalance[1]/_globalRealignmentMass[1] << " / " << _globalRealignmentBalance[2]/_globalRealignmentMass[0] << ") "
			 << "=> adjustment: (" << _universalRealignmentMotion[0] << ", " << _universalRealignmentMotion[1] << ", " << _universalRealignmentMotion[2] << ").\n";
#endif
   
	if(_dropMove){
		outputDropMove(_universalRealignmentMotion[0],_universalRealignmentMotion[2], prefix,_simstep, diff_number);
	}
 	}
   for(ParticleIterator tm = molCont->iteratorBegin(); tm != molCont->iteratorEnd(); ++tm)
   {
     for(unsigned short d=0; d<3; d++){
       tm->move(d, _universalRealignmentMotion[d]);
     }
   }
}

/* by Stefan Becker, borrowed from Matrin Horsch
 * method cancelling the net moment
 */
void Domain::cancelMomentum(
     DomainDecompBase* domainDecomp,
     ParticleContainer* molCont
) {
   double localMomentum[3];
   for(unsigned d = 0; d < 3; d++) localMomentum[d] = 0.0;
   for(ParticleIterator tm = molCont->iteratorBegin(); tm != molCont->iteratorEnd(); ++tm)
   {
      double tmass = tm->gMass();
      for(unsigned d = 0; d < 3; d++)
         localMomentum[d] += tm->v(d) * tmass;
   }
   double globalMomentum[3];
   for(unsigned d = 0; d < 3; d++) globalMomentum[d] = localMomentum[d];
   domainDecomp->collCommInit(3);
   for(unsigned short d = 0; d<3; d++)domainDecomp->collCommAppendDouble(globalMomentum[d]);
   domainDecomp->collCommAllreduceSum();
   for(unsigned short d = 0; d < 3; d++) globalMomentum[d]  = domainDecomp->collCommGetDouble();
   for(unsigned short d = 0; d < 3; d++) globalMomentum[d] /= _globalNumMolecules;
   if(!this->_localRank)
      global_log->info() << "Average momentum: (" << globalMomentum[0] << ", " << globalMomentum[1] << ", " << globalMomentum[2] << ") => removing.\n";
   for(ParticleIterator tm = molCont->iteratorBegin(); tm != molCont->iteratorEnd(); ++tm)
   {
      double tmass = tm->gMass();
      tm->vsub(globalMomentum[0]/tmass, globalMomentum[1]/tmass, globalMomentum[2]/tmass);
   }
}

// author: Stefan Becker, method called in the case of a density profile established in cylindrical coordinates. Counterpart of outputProfile(...).
// reason for a sperate method (in addition to "outputProfile(...)"): method neatly matched to the particular needs of the (cylindrical density) profile, otherwise outpuProfile would be inflated, structure became too compilcated.
void Domain::outputCylProfile(const char* prefix, bool virialProfile){
	if(this->_localRank) return;

	   unsigned IDweight[3];
	   IDweight[0] = this->_universalNProfileUnits[1] * this->_universalNProfileUnits[2];
	   IDweight[1] = this->_universalNProfileUnits[2];
	   IDweight[2] = 1;
	   for( map<unsigned, bool>::iterator pcit = _universalProfiledComponents.begin();    // Loop ueber alle Komponenten
	   	           pcit != _universalProfiledComponents.end();
	   	           pcit++ )
	   	      {
	   	         if(!pcit->second) continue;  // ->second weist auf den key-value von map, d.h. den bool-Wert => falls falsche id, d.h. hiervon ist kein Profil zu erstellen => naechste Schleife

			 // density profile
	   	         string rhoProfName(prefix);
	   	         rhoProfName += ".rhpr";
			 // temperature profile
			 string tmpProfName(prefix);
			 tmpProfName += ".Tpr";
			 //string yVelProfname(prefix);
			 //yVelProfname+= ".vpry";
			 //string xzVelProfname(prefix);
			 //xzVelProfname+= ".vprxz";
			 //virial profile
			string VprProfName(prefix);
			VprProfName += ".Vpr";



	   	         if(!this->_universalCylindricalGeometry)
	   	         {
	   	           global_log->error() << "Incorrect call of method \"outputCylProfile()\" !";
	   	         }
				 // changed by Michaela Heier
	   	         /*ofstream* rhoProf = new ofstream(rhoProfName.c_str());
			 ofstream* tmpProf = new ofstream(tmpProfName.c_str());
			 ofstream* yVelProf = new ofstream(yVelProfname.c_str());
			 ofstream* xzVelProf = new ofstream(xzVelProfname.c_str());

	   	         if (!(*rhoProf ) || !(*tmpProf)) // geaendert durch M. Horsch, by Stefan Becker: wozu?
	   	         {
	   	            return;
	   	         }
	   	     rhoProf->precision(6);
			 tmpProf->precision(6);
			 yVelProf->precision(6);
			 xzVelProf->precision(6);
			 */
			//edited by Michaela Heier
			ofstream rhpr(rhoProfName.c_str());
			ofstream Tpr(tmpProfName.c_str());
			ofstream Vpr(VprProfName.c_str());

			rhpr.precision(6);
			Tpr.precision(6);
			if(virialProfile) Vpr.precision(6);

			//changed by Michaela Heier
	    //##########################################################################################################################################
			 // density profile: actual writing procedure

	   	         // volume of a single bin, in a0^3 (atomic units), adapted to cylindrical coordinates
	   	         //
	   	         double segmentVolume = M_PI / (this->_universalInvProfileUnit[1] * this->_universalInvProfileUnit[2] * this->_universalNProfileUnits[0]);

				rhpr << "//Segment volume: " << segmentVolume << "\n//Accumulated data sets: " << _globalAccumulatedDatasets << "\n//Local profile of the number density. Output file generated by the \"outputCylProfile\" method, located in Domain.cpp. \n";
	   	        rhpr << "//local density profile: Each matrix corresponds to a single value of \"phi_i\", measured in [rad]\n";
	   	        rhpr << "//one single matrix of the local number density rho'(phi_i;r_i',h_i') \n//       | r_i'\n// ---------------------\n//   h_i'| rho'(r_i',h_i')\n//       | \n";
	   	        rhpr << "//  T' \t sigma_ii' \t eps_ii' \t yOffset \t DELTA_phi \t DELTA_r2' \t DELTA_h' \t quantities in atomic units are denoted by an apostrophe '\n";
	   	        rhpr << this->_universalTargetTemperature[_simulation.getEnsemble()->getComponents()->size()]<<"\t"<<_simulation.getEnsemble()->getComponent(0)->getSigma(0)<<"\t"<<_simulation.getEnsemble()->getComponent(0)->getEps(0)<<"\t"; //changed by Michaela Heier
	   	        rhpr << _yOff << "\t" << 1/this->_universalInvProfileUnit[0] << "\t" << 1/this->_universalInvProfileUnit[1] << "\t" << 1/this->_universalInvProfileUnit[2]<< "\n";
	   	         // info: getSigma() und getEps() implementiert in Component.h
	   	         // end of header, start of the data-part of the density file
	   	         for(unsigned n_phi = 0; n_phi < this->_universalNProfileUnits[0]; n_phi++)
	   	         {
	   	        	 rhpr <<"> "<< (n_phi+0.5)/this->_universalInvProfileUnit[0] <<"\n 0 \t";
	   	        	 for(unsigned n_r2 = 0; n_r2 < this->_universalNProfileUnits[1]; n_r2++){
	   	        		 rhpr << 0.5*(sqrt(n_r2+1) + sqrt(n_r2))/sqrt(this->_universalInvProfileUnit[1]) <<"  \t"; // Eintragen der radialen Koordinaten r_i in Header
	   	        	 }
	   	        	 rhpr << "\n";
	   	        	for(unsigned n_h = 0; n_h < this->_universalNProfileUnits[2]; n_h++)
	   	        	{

	   	        		double hval = (n_h + 0.5) / this->_universalInvProfileUnit[2];
	   	        		rhpr << hval<< "  \t";
	   	        		for(unsigned n_r2 = 0; n_r2< this->_universalNProfileUnits[1]; n_r2++)
	   	        		{
	   	        			unsigned unID = n_phi * IDweight[0] + n_r2 * IDweight[1]
	   	        				                                    + n_h * IDweight[2];

	   	        		   	double rho_loc_cyl = this->_universalNProfile[unID] / (segmentVolume * this->_globalAccumulatedDatasets);
	   	        		   	rhpr << rho_loc_cyl << "\t";
	   	        		}
	   	        		rhpr << "\n";
	   	        	}
	   	         }
	   	         //rhoProf->close();
				 rhpr.close();
	   	         //delete rhpr;*/

		//changed by Michaela Heier
	    //##########################################################################################################################################
	   	         // temperature profile: actual writing procedure
			 Tpr << "//Local temperature profile generated by the \"outputCylProfile\" method.\n";
			 Tpr << "//Temperature expressed by 2Ekin/#DOF\n";
			 Tpr << "//one single matrix of the local temperature T(r_i,h_i) \n//      | r_i\n//---------------------\n//  h_i| T(r_i,h_i)\n//      | \n";
	   	     Tpr << "// T \t sigma_ii \t eps_ii \t yOffset \t DELTA_phi \t DELTA_r2 \t DELTA_h \n";
	   	     Tpr << this->_universalTargetTemperature[_simulation.getEnsemble()->getComponents()->size()] <<"\t"<<_simulation.getEnsemble()->getComponent(0)->getSigma(0)<<"\t"<<_simulation.getEnsemble()->getComponent(0)->getEps(0)<<"\t"; //changed by Michaela Heier
	   	     Tpr << _yOff << "\t" << 1/this->_universalInvProfileUnit[0] << "\t" << 1/this->_universalInvProfileUnit[1] << "\t" << 1/this->_universalInvProfileUnit[2]<< "\n";
	   	     Tpr <<"> "<< 0.5/this->_universalInvProfileUnit[0] <<"\n0 \t";
			 for(unsigned n_r2 = 0; n_r2 < this->_universalNProfileUnits[1]; n_r2++){
	   	        		 Tpr << sqrt(n_r2/this->_universalInvProfileUnit[1])+0.5*sqrt(this->_universalInvProfileUnit[1])<<"  \t"; // Eintragen der radialen Koordinaten r_i
	   	          }
	   	        	 Tpr << "\n";
	   	        	for(unsigned n_h = 0; n_h < this->_universalNProfileUnits[2]; n_h++)
	   	        	{
	   	        		double hval = (n_h + 0.5) / this->_universalInvProfileUnit[2];
	   	        		Tpr << hval<< "  \t";
	   	        		for(unsigned n_r2 = 0; n_r2< this->_universalNProfileUnits[1]; n_r2++)
	   	        		{
					  long double DOFc = 0.0;
					  long double twoEkinc = 0.0;
					  for(unsigned n_phi = 0; n_phi < this->_universalNProfileUnits[0]; n_phi++){
					     unsigned unID = n_phi * IDweight[0] + n_r2 * IDweight[1]
	   	        				                                    + n_h * IDweight[2];
					      DOFc += this->_universalDOFProfile[unID];
					      twoEkinc += this->_universalKineticProfile[unID];

					   }
					   if(DOFc == 0.0){
					     Tpr << 0 << "\t";
					   }
					   else{
					 Tpr << (twoEkinc/DOFc ) << "\t";

					   }
	   	        		}
	   	        		Tpr << "\n";
	   	        	}
			  //tmpProf->close();
			  //delete tmpProf;
			  Tpr.close();



			  	//added by Michaela Heier
	    //##########################################################################################################################################
			 // pressure profile: actual writing procedure

	   	         //double segmentVolume; // volume of a single bin, in a0^3 (atomic units)
	   	      	// segmentVolume = M_PI/this->_universalInvProfileUnit[1]/this->_universalInvProfileUnit[2]/this->_universalNProfileUnits[0];  // adapted to cylindrical coordinates
			  if(virialProfile){
				Vpr << "//Local profile of the pressure. Output file generated by the \"outputCylProfile\" method, located in Domain.cpp. \n";
	   	        Vpr << "//local pressure profile: Each matrix corresponds to a single value of \"phi_i\", measured in [rad]\n";
	   	        Vpr << "//one single matrix of the local number pressure p'(phi_i;r_i',h_i') \n//      | r_i'\n//---------------------\n//  h_i'| p'(r_i',h_i')\n//      | \n";
	   	        Vpr << "// T' \t sigma_ii' \t eps_ii' \t yOffset \t DELTA_phi \t DELTA_r2' \t DELTA_h' \t quantities in atomic units are denoted by an apostrophe '\n";
	   	        Vpr << this->_universalTargetTemperature[_simulation.getEnsemble()->getComponents()->size()]<<"\t"<<_simulation.getEnsemble()->getComponent(0)->getSigma(0)<<"\t"<<_simulation.getEnsemble()->getComponent(0)->getEps(0)<<"\t";
	   	        Vpr << _yOff << "\t" << 1/this->_universalInvProfileUnit[0] << "\t" << 1/this->_universalInvProfileUnit[1] << "\t" << 1/this->_universalInvProfileUnit[2]<< "\n";
	   	         // info: getSigma() und getEps() implementiert in Component.h
	   	         // end of header, start of the data-part of the density file
	   	         for(unsigned n_phi = 0; n_phi < this->_universalNProfileUnits[0]; n_phi++)
	   	         {
	   	        	 Vpr <<"> "<< (n_phi+0.5)/this->_universalInvProfileUnit[0] <<"\n0 \t";
	   	        	 for(unsigned n_r2 = 0; n_r2 < this->_universalNProfileUnits[1]; n_r2++){
	   	        		 Vpr << 0.5*(sqrt(n_r2+1) + sqrt(n_r2))/sqrt(this->_universalInvProfileUnit[1]) <<"  \t"; // Eintragen der radialen Koordinaten r_i in Header
	   	        	 }
	   	        	 Vpr << "\n";
	   	        	for(unsigned n_h = 0; n_h < this->_universalNProfileUnits[2]; n_h++)
	   	        	{

	   	        		double hval = (n_h + 0.5) / this->_universalInvProfileUnit[2];
	   	        		Vpr << hval<< "  \t";
						//long double Px = 0.0;
						//long double Py = 0.0;
						//long double Pz = 0.0;
	   	        		for(unsigned n_r2 = 0; n_r2< this->_universalNProfileUnits[1]; n_r2++)
	   	        		{
	   	        			unsigned unID = n_phi * IDweight[0] + n_r2 * IDweight[1]
	   	        				                                    + n_h * IDweight[2];
							//Px += this->_universalPXProfile[unID];
							//Py += this->_universalPYProfile[unID];
							//Pz += this->_universalPZProfile[unID];

	   	        		   	double virial_cyl =(_globalTemperatureMap[0]*this->_universalNProfile[unID]+this->_universalPYProfile[unID]) / (segmentVolume * this->_globalAccumulatedDatasets);
	   	        		   	Vpr << virial_cyl << "\t";
						}

	   	        		Vpr << "\n";
	   	        	}
	   	         }
	   	         //rhoProf->close();
				 Vpr.close();
			  }
	   	         //delete rhpr;*/
				 }

}



 void Domain::setLocalUpotCompSpecific(double UpotCspec){_localUpotCspecif = UpotCspec;}

 double Domain::getLocalUpotCompSpecific(){return _localUpotCspecif;}


double Domain::getAverageGlobalUpotCSpec() {
  global_log->debug() << "number of fluid molecules = " << getNumFluidMolecules() << "\n";
  return _globalUpotCspecif / getNumFluidMolecules();
}


void Domain::setNumFluidComponents(unsigned nc){_numFluidComponent = nc;}

unsigned Domain::getNumFluidComponents(){return _numFluidComponent;}

unsigned long Domain::getNumFluidMolecules(){
  unsigned long numFluidMolecules = 0;
  for(unsigned i = 0; i < _numFluidComponent; i++){
    Component& ci=*(global_simulation->getEnsemble()->getComponent(i));
    numFluidMolecules+=ci.getNumMolecules();
  }
  return numFluidMolecules;
}



