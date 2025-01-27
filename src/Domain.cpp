
#include <iostream>
#include <string>
#include <cmath>
#include <cstdint>

#include "Domain.h"

#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
//#include "CutoffCorrections.h"
#include "Simulation.h"
#include "ensemble/EnsembleBase.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "utils/FileUtils.h"
#include "utils/Logger.h"
#include "utils/arrayMath.h"
#include "utils/CommVar.h"



Domain::Domain(int rank) {
	_localRank = rank;
	_localUpot = 0;
	_localVirial = 0;
	_globalUpot = 0;
	_globalVirial = 0;
	_globalRho = 0;

	this->_componentToThermostatIdMap = std::map<int, int>();
	this->_localThermostatN = std::map<int, unsigned long>();
	this->_localThermostatN[-1] = 0;
	this->_localThermostatN[0] = 0;
	this->_universalThermostatN = std::map<int, unsigned long>();
	this->_universalThermostatN[-1] = 0;
	this->_universalThermostatN[0] = 0;
	this->_localRotationalDOF = std::map<int, unsigned long>();
	this->_localRotationalDOF[-1] = 0;
	this->_localRotationalDOF[0] = 0;
	this->_universalRotationalDOF = std::map<int, unsigned long>();
	this->_universalRotationalDOF[-1] = 0;
	this->_universalRotationalDOF[0] = 0;
	this->_globalLength[0] = 0;
	this->_globalLength[1] = 0;
	this->_globalLength[2] = 0;
	this->_universalBTrans = std::map<int, double>();
	this->_universalBTrans[0] = 1.0;
	this->_universalBRot = std::map<int, double>();
	this->_universalBRot[0] = 1.0;
	this->_universalTargetTemperature = std::map<int, double>();
	this->_universalTargetTemperature[0] = 1.0;
	this->_globalTemperatureMap = std::map<int, double>();
	this->_globalTemperatureMap[0] = 1.0;
	this->_local2KETrans[0] = 0.0;
	this->_local2KERot[0] = 0.0;

	this->_universalNVE = false;
	this->_globalUSteps = 0;
	this->_globalSigmaU = 0.0;
	this->_globalSigmaUU = 0.0;
	this->_componentwiseThermostat = false;
#ifdef COMPLEX_POTENTIAL_SET
	this->_universalUndirectedThermostat = std::map<int, bool>();
	this->_universalThermostatDirectedVelocity = std::map<int, std::array<double,3> >();
	this->_localThermostatDirectedVelocity = std::map<int, std::array<double,3> >();
#endif
	this->_universalSelectiveThermostatCounter = 0;
	this->_universalSelectiveThermostatWarning = 0;
	this->_universalSelectiveThermostatError = 0;

    // explosion heuristics, NOTE: turn off when using slab thermostat
    _bDoExplosionHeuristics = true;
}

void Domain::readXML(XMLfileUnits& xmlconfig) {
	/* volume */
	if ( xmlconfig.changecurrentnode( "volume" )) {
		std::string type;
		xmlconfig.getNodeValue( "@type", type );
		Log::global_log->info() << "Volume type: " << type << std::endl;
		if( type == "box" ) {
			xmlconfig.getNodeValueReduced( "lx", _globalLength[0] );
			xmlconfig.getNodeValueReduced( "ly", _globalLength[1] );
			xmlconfig.getNodeValueReduced( "lz", _globalLength[2] );
			Log::global_log->info() << "Box size: " << _globalLength[0] << ", "
				<< _globalLength[1] << ", "
				<< _globalLength[2] << std::endl;
		}
		else {
			Log::global_log->error() << "Unsupported volume type " << type << std::endl;
		}
		xmlconfig.changecurrentnode("..");
	}

	/* temperature */
	/** @todo reading temperature is performed in the Ensemble, so check and remove */
	bool bInputOk = true;
	double temperature = 0.;
	bInputOk = bInputOk && xmlconfig.getNodeValueReduced("temperature", temperature);
	if(bInputOk) {
		setGlobalTemperature(temperature);
		Log::global_log->info() << "Temperature: " << temperature << std::endl;
	}
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
	Log::global_log->debug() << "* local thermostat " << thermostat << ":  mvv = " << summv2 << std::endl;
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

double Domain::getAverageGlobalVirial() { return _globalVirial/_globalNumMolecules; }

double Domain::getAverageGlobalUpot() { return getGlobalUpot()/_globalNumMolecules; }
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

#ifdef ENABLE_PERSISTENT
	auto collComm = makeCollCommObjAllreduceAdd(domainDecomp->getCommunicator(), Upot, Virial);
	collComm.persistent();
	collComm.get(Upot, Virial);
#else
	/* FIXME stuff for the ensemble class */
	domainDecomp->collCommInit(2, 654);
	domainDecomp->collCommAppendDouble(Upot);
	domainDecomp->collCommAppendDouble(Virial);
	domainDecomp->collCommAllreduceSumAllowPrevious();
	Upot = domainDecomp->collCommGetDouble();
	Virial = domainDecomp->collCommGetDouble();
	domainDecomp->collCommFinalize();
#endif

	// Process 0 has to add the dipole correction:
	// m_UpotCorr and m_VirialCorr already contain constant (internal) dipole correction
	_globalUpot = Upot + _UpotCorr;
	_globalVirial = Virial + _VirialCorr;

	/*
	 * thermostat ID 0 represents the entire system
	 */

	std::map<int, unsigned long>::iterator thermit;
	if( _componentwiseThermostat )
	{
		Log::global_log->debug() << "* applying a component-wise thermostat" << std::endl;
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
	int thermid = 0;
	for (thermit = _universalThermostatN.begin(); thermit != _universalThermostatN.end(); thermit++, thermid++)
	{
		// number of molecules on the local process. After the reduce operation
		// num_molecules will contain the global number of molecules
		unsigned long numMolecules = _localThermostatN[thermit->first];
		double summv2 = _local2KETrans[thermit->first];
		unsigned long rotDOF = _localRotationalDOF[thermit->first];
		double sumIw2 = (rotDOF > 0)? _local2KERot[thermit->first]: 0.0;

#ifdef ENABLE_PERSISTENT
	auto collComm = makeCollCommObjAllreduceAdd(domainDecomp->getCommunicator(), summv2, sumIw2, numMolecules, rotDOF);
	collComm.persistent();
	collComm.get(summv2, sumIw2, numMolecules, rotDOF);
#else
		domainDecomp->collCommInit(4, 12+thermid);
		domainDecomp->collCommAppendDouble(summv2);
		domainDecomp->collCommAppendDouble(sumIw2);
		domainDecomp->collCommAppendUnsLong(numMolecules);
		domainDecomp->collCommAppendUnsLong(rotDOF);
		domainDecomp->collCommAllreduceSumAllowPrevious();
		summv2 = domainDecomp->collCommGetDouble();
		sumIw2 = domainDecomp->collCommGetDouble();
		numMolecules = domainDecomp->collCommGetUnsLong();
		rotDOF = domainDecomp->collCommGetUnsLong();
		domainDecomp->collCommFinalize();
#endif
		Log::global_log->debug() << "[ thermostat ID " << thermit->first << "]\tN = " << numMolecules << "\trotDOF = " << rotDOF
			<< "\tmv2 = " <<  summv2 << "\tIw2 = " << sumIw2 << std::endl;

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
			Log::global_log->warning() << "Explosion!" << std::endl;
			Log::global_log->debug() << "Selective thermostat will be applied to set " << thermit->first
				<< " (beta_trans = " << this->_universalBTrans[thermit->first]
				<< ", beta_rot = " << this->_universalBRot[thermit->first] << "!)" << std::endl;
			const double limit_energy =  KINLIMIT_PER_T * Ti;

			#if defined(_OPENMP)
			#pragma omp parallel
			#endif
			{
				for (auto tM = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tM.isValid(); ++tM) {
					const auto Utrans = tM->U_trans();
					if (Utrans > limit_energy) {
						const auto vcorr = sqrt(limit_energy / Utrans);
						Log::global_log->debug() << ": v(m" << tM->getID() << ") *= " << vcorr << std::endl;
						tM->scale_v(vcorr);
						tM->scale_F(vcorr);
					}

					const auto rot_dof = tM->component()->getRotationalDegreesOfFreedom();
					if (rot_dof > 0) {
						const auto limit_rot_energy = 3.0 * rot_dof * Ti;
						const auto Urot = tM->U_rot();
						if (Urot > limit_rot_energy) {
							const auto Dcorr = sqrt(limit_rot_energy / Urot);
							Log::global_log->debug() << "D(m" << tM->getID() << ") *= " << Dcorr << std::endl;
							tM->scale_D(Dcorr);
							tM->scale_M(Dcorr);
						}
					}
				}
			} /*_OPENMP*/

			// arbitrary values set by one of the thermodynamic guys (probs. Martin Horsch) :)
			int explosionReappearanceLimit = 4000;
			int explosionVanishGracePeriod = 40;
			int stepsSinceLastExplosion = explosionReappearanceLimit - _universalSelectiveThermostatCounter;
			// We set warning to true if the explosion is not gone after 40 steps or if it reappears within 4000 steps.
			// If it still persists after 80 steps or if it reappears twice with no more than 4000 steps between the
			// occurrences, we set error to true.
			// These counters are reduced in every step, s.t., the warning vanishes after 4000 steps without explosions.
			if (stepsSinceLastExplosion >= explosionVanishGracePeriod) {
				if( _universalSelectiveThermostatWarning > 0 )
					_universalSelectiveThermostatError = _universalSelectiveThermostatWarning;
				if( _universalSelectiveThermostatCounter > 0 )
					_universalSelectiveThermostatWarning = _universalSelectiveThermostatCounter;
				_universalSelectiveThermostatCounter = explosionReappearanceLimit;
			}
			_universalBTrans[thermit->first] = 1.0;
			_universalBRot[thermit->first] = pow(this->_universalBRot[thermit->first], 0.0091);
		}
#ifdef NDEBUG
		if( (_universalSelectiveThermostatCounter > 0) &&
				((_universalSelectiveThermostatCounter % 20) == 10) )
#endif
			/* FIXME: why difference counters? */
			Log::global_log->debug() << "counter " << _universalSelectiveThermostatCounter
				<< ",\t warning " << _universalSelectiveThermostatWarning
				<< ",\t error " << _universalSelectiveThermostatError << std::endl;

		if(collectThermostatVelocities && _universalUndirectedThermostat[thermit->first])
		{
			std::array<double, 3> sigv = _localThermostatDirectedVelocity[thermit->first];

#ifdef ENABLE_PERSISTENT
	auto collComm = makeCollCommObjAllreduceAdd(domainDecomp->getCommunicator(), sigv[0], sigv[1], sigv[2]);
	collComm.persistent();
	collComm.get(sigv[0], sigv[1], sigv[2]);
#else
			domainDecomp->collCommInit(3);
			for(int d=0; d < 3; d++) domainDecomp->collCommAppendDouble(sigv[d]);
			domainDecomp->collCommAllreduceSum();
			for(int d=0; d < 3; d++) sigv[d] = domainDecomp->collCommGetDouble();
			domainDecomp->collCommFinalize();
#endif

			_localThermostatDirectedVelocity[thermit->first].fill(0.0);

			if(numMolecules > 0)
				_universalThermostatDirectedVelocity[thermit->first] = arrayMath::mulScalar(sigv, 1.0 / static_cast<double>(numMolecules));
			else
				_universalThermostatDirectedVelocity[thermit->first].fill(0.0);

#ifndef NDEBUG
			Log::global_log->debug() << "* thermostat " << thermit->first
				<< " directed velocity: ("
				<< _universalThermostatDirectedVelocity[thermit->first][0]
				<< " / " << _universalThermostatDirectedVelocity[thermit->first][1]
				<< " / " << _universalThermostatDirectedVelocity[thermit->first][2]
				<< ")" << std::endl;
#endif
		}

#ifndef NDEBUG
		Log::global_log->debug() << "* Th" << thermit->first << " N=" << numMolecules
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
	if(this->_componentwiseThermostat)
	{
		for( std::map<int, bool>::iterator thit = _universalUndirectedThermostat.begin();
				thit != _universalUndirectedThermostat.end();
				thit ++ )
		{
			if(thit->second)
				_localThermostatDirectedVelocity[thit->first].fill(0.0);
		}

		#if defined(_OPENMP)
		#pragma omp parallel
		#endif
		{
			std::map<int, std::array<double, 3> > localThermostatDirectedVelocity_thread;

			for(auto tM = partCont->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tM.isValid(); ++tM) {
				int cid = tM->componentid();
				int thermostat = this->_componentToThermostatIdMap[cid];

				if(this->_universalUndirectedThermostat[thermostat])
					arrayMath::accumulate(localThermostatDirectedVelocity_thread[thermostat], tM->v_arr());
			}

			#if defined(_OPENMP)
			#pragma omp critical(collectVelocities1111)
			#endif
			{
				for (auto it = localThermostatDirectedVelocity_thread.begin(); it != localThermostatDirectedVelocity_thread.end(); ++it) {
					arrayMath::accumulate(_localThermostatDirectedVelocity[it->first], it->second);
				}
			}

		}
	}
	else if(this->_universalUndirectedThermostat[0])
	{
		double velX = 0.0, velY = 0.0, velZ = 0.0;

		#if defined(_OPENMP)
		#pragma omp parallel reduction(+ : velX, velY, velZ)
		#endif
		{
			for(auto tM = partCont->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tM.isValid(); ++tM) {
				velX += tM->v(0);
				velY += tM->v(1);
				velZ += tM->v(2);
			}
		}

		_localThermostatDirectedVelocity[0][0] = velX;
		_localThermostatDirectedVelocity[0][1] = velY;
		_localThermostatDirectedVelocity[0][2] = velZ;
	}
}

void Domain::calculateVelocitySums(ParticleContainer* partCont)
{
	if(this->_componentwiseThermostat)
	{
		for(auto tM = partCont->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tM.isValid(); ++tM)
		{
			int cid = tM->componentid();
			int thermostat = this->_componentToThermostatIdMap[cid];
			this->_localThermostatN[thermostat]++;
			this->_localRotationalDOF[thermostat] += tM->component()->getRotationalDegreesOfFreedom();
			if(this->_universalUndirectedThermostat[thermostat])
			{
				tM->calculate_mv2_Iw2( this->_local2KETrans[thermostat],
						this->_local2KERot[thermostat],
						this->_universalThermostatDirectedVelocity[thermostat][0],
						this->_universalThermostatDirectedVelocity[thermostat][1],
						this->_universalThermostatDirectedVelocity[thermostat][2]  );
			}
			else
			{
				tM->calculate_mv2_Iw2(_local2KETrans[thermostat], _local2KERot[thermostat]);
			}
		}
	}
	else
	{
		unsigned long N = 0, rotationalDOF = 0;
		double local2KETrans = 0.0, local2KERot = 0.0;
		#if defined(_OPENMP)
		#pragma omp parallel reduction(+:N, rotationalDOF, local2KETrans, local2KERot)
		#endif
		{

			for(auto tM = partCont->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tM.isValid(); ++tM) {
				++N;
				rotationalDOF += tM->component()->getRotationalDegreesOfFreedom();
				if(this->_universalUndirectedThermostat[0]) {
					tM->calculate_mv2_Iw2( local2KETrans,
							local2KERot,
							this->_universalThermostatDirectedVelocity[0][0],
							this->_universalThermostatDirectedVelocity[0][1],
							this->_universalThermostatDirectedVelocity[0][2]  );
				} else {
					tM->calculate_mv2_Iw2(local2KETrans, local2KERot);
				}
			}
		} /* _OPENMP */

		this->_localThermostatN[0] = N;
		this->_localRotationalDOF[0] = rotationalDOF;
		this->_local2KETrans[0] = local2KETrans;
		this->_local2KERot[0] = local2KERot;

		Log::global_log->debug() << "      * N = " << this->_localThermostatN[0]
			<< " rotDOF = " << this->_localRotationalDOF[0] << "   mv2 = "
			<< _local2KETrans[0] << " Iw2 = " << _local2KERot[0] << std::endl;
	}
}

void Domain::writeCheckpointHeader(std::string filename,
		ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, double currentTime) {
		unsigned long globalNumMolecules = this->getglobalNumMolecules(true, particleContainer, domainDecomp);
		/* Rank 0 writes file header */
		if(0 == this->_localRank) {
			std::ofstream checkpointfilestream(filename.c_str());
			checkpointfilestream << "mardyn trunk " << CHECKPOINT_FILE_VERSION;
			checkpointfilestream << "\n";  // store default format flags
			std::ios::fmtflags f( checkpointfilestream.flags() );
			checkpointfilestream << "currentTime\t"  << FORMAT_SCI_MAX_DIGITS << currentTime << "\n"; //edited by Michaela Heier
			checkpointfilestream.flags(f);  // restore default format flags
			checkpointfilestream << " Length\t" << std::setprecision(9) << _globalLength[0] << " " << _globalLength[1] << " " << _globalLength[2] << "\n";
			if(this->_componentwiseThermostat)
			{
				for( std::map<int, int>::iterator thermit = this->_componentToThermostatIdMap.begin();
						thermit != this->_componentToThermostatIdMap.end();
						thermit++ )
				{
					if(0 >= thermit->second) continue;
					checkpointfilestream << " CT\t" << 1+thermit->first
						<< "\t" << thermit->second << "\n";
				}
				for( std::map<int, double>::iterator Tit = this->_universalTargetTemperature.begin();
						Tit != this->_universalTargetTemperature.end();
						Tit++ )
				{
					if((0 >= Tit->first) || (0 >= Tit->second)) continue;
					checkpointfilestream << " ThT " << Tit->first << "\t" << Tit->second << "\n";
				}
			}
			else
			{
				checkpointfilestream << " Temperature\t" << _universalTargetTemperature[0] << std::endl;
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
				checkpointfilestream << std::setprecision(13);
				checkpointfilestream << " I\t" << this->_globalUSteps << " "
					<< this->_globalSigmaU << " " << this->_globalSigmaUU << "\n";
				checkpointfilestream << std::setprecision(8);
			}
			*/
			std::vector<Component>* components = _simulation.getEnsemble()->getComponents();
			checkpointfilestream << " NumberOfComponents\t" << components->size() << std::endl;
			for(auto pos=components->begin();pos!=components->end();++pos){
				pos->write(checkpointfilestream);
			}
			unsigned int numperline=_simulation.getEnsemble()->getComponents()->size();
			unsigned int iout=0;
			for(auto pos=_mixcoeff.begin();pos!=_mixcoeff.end();++pos){
				checkpointfilestream << *pos;
				iout++;
				// 2 parameters (xi and eta)
				if(iout/2>=numperline) {
					checkpointfilestream << std::endl;
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
			checkpointfilestream << _epsilonRF << std::endl;
			for( auto uutit = this->_universalUndirectedThermostat.begin();
					uutit != this->_universalUndirectedThermostat.end();
					uutit++ )
			{
				if(0 > uutit->first) continue;
				if(uutit->second) checkpointfilestream << " U\t" << uutit->first << "\n";
			}
			checkpointfilestream << " NumberOfMolecules\t" << globalNumMolecules << std::endl;

			checkpointfilestream << " MoleculeFormat\t" << Molecule::getWriteFormat() << std::endl;
			checkpointfilestream.close();
		}

}

void Domain::writeCheckpointHeaderXML(std::string filename, ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, double currentTime)
{
	unsigned long globalNumMolecules = this->getglobalNumMolecules(true, particleContainer, domainDecomp);
	if(0 != domainDecomp->getRank() )
		return;

	std::ofstream ofs(filename.c_str() );

	ofs << "<?xml version='1.0' encoding='UTF-8'?>" << std::endl;
	ofs << "<mardyn version=\"20100525\" >" << std::endl;
	ofs << "\t<headerinfo>" << std::endl;
	std::ios::fmtflags f( ofs.flags() );
	ofs << "\t\t<time>" << FORMAT_SCI_MAX_DIGITS_WIDTH_21 << currentTime << "</time>" << std::endl;
	ofs << "\t\t<length>" << std::endl;
	ofs << "\t\t\t<x>" << FORMAT_SCI_MAX_DIGITS_WIDTH_21 << _globalLength[0] << "</x> "
				 "<y>" << FORMAT_SCI_MAX_DIGITS_WIDTH_21 << _globalLength[1] << "</y> "
				 "<z>" << FORMAT_SCI_MAX_DIGITS_WIDTH_21 << _globalLength[2] << "</z>" << std::endl;
	ofs << "\t\t</length>" << std::endl;
	ofs.flags(f);  // restore default format flags
	ofs << "\t\t<number>" << globalNumMolecules << "</number>" << std::endl;
	ofs << "\t\t<format type=\"" << Molecule::getWriteFormat() << "\"/>" << std::endl;
	ofs << "\t</headerinfo>" << std::endl;
	ofs << "</mardyn>" << std::endl;
}

void Domain::writeCheckpoint(std::string filename,
		ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, double currentTime,
		bool useBinaryFormat) {
	domainDecomp->assertDisjunctivity(particleContainer);
#ifdef ENABLE_REDUCED_MEMORY_MODE
	Log::global_log->warning() << "The checkpoints are not adapted for RMM-mode. Velocity will be one half-timestep ahead!" << std::endl;
	Log::global_log->warning() << "See Domain::writeCheckpoint() for a suggested workaround." << std::endl;
	//TODO: desired correctness (compatibility to normal mode) should be achievable by:
	// 1. integrating positions by half a timestep forward (+ delta T / 2)
	// 2. writing the checkpoint (with currentTime + delta T ? )
	// 3. integrating positions by half a timestep backward (- delta T / 2)
#endif

	if (useBinaryFormat) {
		this->writeCheckpointHeaderXML((filename + ".header.xml"), particleContainer, domainDecomp, currentTime);
	} else {
		this->writeCheckpointHeader(filename, particleContainer, domainDecomp, currentTime);
	}
	domainDecomp->writeMoleculesToFile(filename, particleContainer, useBinaryFormat);
}


void Domain::initParameterStreams(double cutoffRadius, double cutoffRadiusLJ){
	_comp2params.initialize(*(_simulation.getEnsemble()->getComponents()), _mixcoeff, _epsilonRF, cutoffRadius, cutoffRadiusLJ);
}

void Domain::Nadd(unsigned cid, int N, int localN)
{
	Ensemble* ensemble = _simulation.getEnsemble();
	Component* component = ensemble->getComponent(cid);
	component->incNumMolecules(N);
	unsigned long rotationDegreesOfFreeedom = static_cast<unsigned long>(component->getRotationalDegreesOfFreedom());

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
	this->_globalRho = this->getglobalNumMolecules(true, nullptr, domainDecomp) /
		(this->_globalLength[0] * this->_globalLength[1] * this->_globalLength[2]);
}

void Domain::setTargetTemperature(int thermostatID, double targetT)
{
	if(thermostatID < 0)
	{
		Log::global_log->warning() << "Warning: thermostat \'" << thermostatID << "\' (T = "
			<< targetT << ") will be ignored." << std::endl;
		return;
	}

	this->_universalTargetTemperature[thermostatID] = targetT;
	if(!(this->_universalUndirectedThermostat[thermostatID] == true))
		this->_universalUndirectedThermostat[thermostatID] = false;

	/* FIXME: Substantial change in program behavior! */
	if(thermostatID == 0) {
#ifndef UNIT_TESTS
		Log::global_log->warning() << "Disabling the component wise thermostat!" << std::endl;
#endif
		disableComponentwiseThermostat();
	}
	if(thermostatID >= 1) {
		if( ! _componentwiseThermostat ) {
			/* FIXME: Substantial change in program behavior! */
			Log::global_log->warning() << "Enabling the component wise thermostat!" << std::endl;
			_componentwiseThermostat = true;
			_universalTargetTemperature.erase(0);
			_universalUndirectedThermostat.erase(0);
			this->_universalThermostatDirectedVelocity.erase(0);
			std::vector<Component>* components = _simulation.getEnsemble()->getComponents();
			for( std::vector<Component>::iterator tc = components->begin(); tc != components->end(); tc ++ ) {
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
	std::vector<Component>* components = _simulation.getEnsemble()->getComponents();
	for( std::vector<Component>::iterator tc = components->begin(); tc != components->end(); tc ++ ) {
		if(!(this->_componentToThermostatIdMap[ tc->ID() ] > 0)) {
			this->_componentToThermostatIdMap[ tc->ID() ] = -1;
		}
	}
}

void Domain::enableUndirectedThermostat(int tst)
{
	this->_universalUndirectedThermostat[tst] = true;
	this->_localThermostatDirectedVelocity[tst].fill(0.0);
	this->_universalThermostatDirectedVelocity[tst].fill(0.0);
}

std::vector<double> & Domain::getmixcoeff() { return _mixcoeff; }

double Domain::getepsilonRF() const { return _epsilonRF; }

void Domain::setepsilonRF(double erf) { _epsilonRF = erf; }

unsigned long Domain::getglobalNumMolecules(bool bUpdate, ParticleContainer* particleContainer,
											DomainDecompBase* domainDecomp) {
	if (bUpdate) {
		if (particleContainer == nullptr) {
			Log::global_log->debug() << "Domain::getglobalNumMolecules: Passed Particle Container is null! Fetching pointer "
								   "from global_simulation."
								<< std::endl;
			particleContainer = global_simulation->getMoleculeContainer();
		}
		if (domainDecomp == nullptr) {
			Log::global_log->debug() << "Domain::getglobalNumMolecules: Passed Domain Decomposition is null! Fetching "
								   "pointer from global_simulation."
								<< std::endl;
			domainDecomp = &(global_simulation->domainDecomposition());
		}
		this->updateglobalNumMolecules(particleContainer, domainDecomp);
	}
	return _globalNumMolecules;
}

void Domain::setglobalNumMolecules(unsigned long glnummol) { _globalNumMolecules = glnummol; }

void Domain::updateglobalNumMolecules(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp) {
	unsigned long oldNum = _globalNumMolecules;
	CommVar<uint64_t> numMolecules;
	numMolecules.local = particleContainer->getNumberOfParticles(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
#ifdef ENABLE_MPI

#ifdef ENABLE_PERSISTENT
	auto collComm = makeCollCommObjAllreduceAdd(domainDecomp->getCommunicator(), numMolecules.local);
	collComm.persistent();
	collComm.get(numMolecules.global);
#else
	domainDecomp->collCommInit(1);
	domainDecomp->collCommAppendUnsLong(numMolecules.local);
	domainDecomp->collCommAllreduceSum();
	numMolecules.global = domainDecomp->collCommGetUnsLong();
	domainDecomp->collCommFinalize();
#endif
#else
	numMolecules.global = numMolecules.local;
#endif
	this->setglobalNumMolecules(numMolecules.global);
	Log::global_log->debug() << "Updated global number of particles from " << oldNum << " to N_new = " << _globalNumMolecules << std::endl;
}

CommVar<uint64_t> Domain::getMaxMoleculeID() const {
	return _maxMoleculeID;
}

void Domain::updateMaxMoleculeID(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp)
{
	_maxMoleculeID.local = 0;
	for(auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit) {
		uint64_t pid = pit->getID();
		if(pid > _maxMoleculeID.local)
			_maxMoleculeID.local = pid;
	}
#ifdef ENABLE_MPI
	MPI_Allreduce(&_maxMoleculeID.local, &_maxMoleculeID.global, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
#else
	_maxMoleculeID.global = _maxMoleculeID.local;
#endif
}

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

// Profiling done in the Domain class anymore. Use SpatialProfile to add profile functionalities.
void Domain::submitDU(unsigned /*cid*/, double DU, double* r)
{
	/*unsigned xun, yun, zun;
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
	}*/
}

 void Domain::setLocalUpotCompSpecific(double UpotCspec){_localUpotCspecif = UpotCspec;}

 double Domain::getLocalUpotCompSpecific(){return _localUpotCspecif;}


double Domain::getAverageGlobalUpotCSpec() {
  Log::global_log->debug() << "number of fluid molecules = " << getNumFluidMolecules() << "\n";
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
