/*
 * DirectedPM.h
 *
 *  Created on: 25 august 2020
 *      Author: Koch
 */

#include "DirectedPM.h"

//! @brief will be called to read configuration
//!
//!
//!
//! \param xmlconfig  read from config.xml
void DirectedPM::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("Component", _component);      // Component
	xmlconfig.getNodeValue("rIncrements", _rIncrements);  // Bin Size
	xmlconfig.getNodeValue("hIncrements", _hIncrements);
	xmlconfig.getNodeValue("phiIncrements", _phiIncrements);
	// calculation which bins are in droplet and which are in vapor phase --> determination by density
	// densities are used as starting values
	xmlconfig.getNodeValue(
		"rohCutLiq",
		_rohCutLiq);  // Minimum density of liquid phase, for higher densities the bins are declared as liquid phase
	xmlconfig.getNodeValue("maxDeviation",
						   _percent);  // deviation of liquid and vapor phase density which is allowed for the
									   // calculation of the mean phase density given as a percentage
	xmlconfig.getNodeValue("heightWall", _heightWall);  // height for neglecting the adsorbate layer next to the wall
	xmlconfig.getNodeValue("heightMembrane", _heightMembrane);    // height for neglecting the influence of the membrane
	xmlconfig.getNodeValue("outputFrequency", _outputFrequency);  // Averaged time steps

	global_log->info() << "[DirectedPM] settings:" << std::endl;
	global_log->info() << "                  Component: " << _component << std::endl;
	global_log->info() << "                  r: " << _rIncrements << std::endl;
	global_log->info() << "                  h: " << _hIncrements << std::endl;
	global_log->info() << "                  phi: " << _phiIncrements << std::endl;
	global_log->info() << "                  rohCutLiq: " << _rohCutLiq << std::endl;
	global_log->info() << "                  percent: " << _percent << std::endl;
	global_log->info() << "                  heightWall: " << _heightWall << std::endl;
	global_log->info() << "                  heightMembrane: " << _heightMembrane << std::endl;
	global_log->info() << "                outputFrequency: " << _outputFrequency << std::endl;
}

void DirectedPM::beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
							  unsigned long simstep) {
	if (_enabled) {
		global_log->debug() << "[DirectedPM] before forces called" << std::endl;

		// CALCULATE SYNSIMSTEP FOR MATRIX
		double synsimstep;
		if (simstep % _outputFrequency != 0) {
			synsimstep = (simstep % _outputFrequency);
		} else {
			synsimstep = _outputFrequency;
		}

		double rohCutLiqActual = 0.;

		// BOX ASSIGNMENT HELP
		double rUN = 0.;
		double hUN = 0.;
		double phiUN = 0.;
		double xc = 0.;
		double yc = 0.;
		double zc = 0.;
		double R2 = 0.;
		double phi = 0.;
		unsigned unID = 0;
		unsigned COMPONENT;
		// COUNT ALLL PARTICLES OF THE SELCETED COMPONENT
		// double Particles = 0.;

		_iterationsSinceStart += 1;
		_xyzEkinDroplet[synsimstep] = 0.;
		_xzEkinDroplet[synsimstep] = 0.;
		_xyzEkinGas[synsimstep] = 0.;
		_xzEkinGas[synsimstep] = 0.;
		if (simstep > 0) {
			// ITERATE OVER PARTICLES AND GET VELOCITYS AND VIRIAL IN EVERY DIMENSION
			for (auto temporaryMolecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
				 temporaryMolecule.isValid(); ++temporaryMolecule) {
				COMPONENT = temporaryMolecule->getComponentLookUpID();
				if (COMPONENT == (_component - 1)) {
					// ASSIGN THE PARTICLES TO A UNID
					xc = temporaryMolecule->r(0) - _universalCentre[0];
					yc = temporaryMolecule->r(1) - _universalCentre[1];
					zc = temporaryMolecule->r(2) - _universalCentre[2];
					R2 = xc * xc + zc * zc;
					phi = atan2(zc, xc);
					if (phi < 0.0) {
						phi = phi + 2.0 * M_PI;
					}
					rUN = floor(R2 * _universalInvProfileUnit[0]);
					hUN = floor(yc * _universalInvProfileUnit[1]);
					phiUN = floor(phi * _universalInvProfileUnit[2]);

					// Check if R is inside cylinder
					if (rUN < _rIncrements) {
						if ((rUN >= 0) && (hUN >= 0) && (phiUN >= 0) && (rUN < _rIncrements) && (hUN < _hIncrements) &&
							(phiUN < _phiIncrements)) {
							unID = (hUN * _rIncrements * _phiIncrements) + (rUN * _phiIncrements) + phiUN;
						} else {
							global_log->error()
								<< "INV PROFILE UNITS " << _universalInvProfileUnit[0] << " "
								<< _universalInvProfileUnit[1] << " " << _universalInvProfileUnit[2] << "\n";
							global_log->error() << "PROFILE UNITS " << _rIncrements << " " << _hIncrements << " "
												<< _phiIncrements << "\n";
							global_log->error() << "Severe error!! Invalid profile ID (" << rUN << " / " << hUN << " / "
												<< phiUN << ").\n\n";
							global_log->error() << "Severe error!! Invalid profile unit (" << R2 << " / " << yc << " / "
												<< phi << ").\n\n";
							global_log->error()
								<< "Coordinates off center (" << xc << " / " << yc << " / " << zc << ").\n";
							global_log->error() << "unID = " << unID << "\n";
							Simulation::exit(707);
						}
						// ADD VELOCITCY AND VIRIAL TO RESPECTIVE BIN
						_localnumberOfParticles[unID] += 1.;
						_localXyzVelocities[unID][0] += temporaryMolecule->v(0);
						_localXyzVelocities[unID][1] += temporaryMolecule->v(1);
						_localXyzVelocities[unID][2] += temporaryMolecule->v(2);
						_localXyzVelocities2[unID][0] += temporaryMolecule->v(0) * temporaryMolecule->v(0);
						_localXyzVelocities2[unID][1] += temporaryMolecule->v(1) * temporaryMolecule->v(1);
						_localXyzVelocities2[unID][2] += temporaryMolecule->v(2) * temporaryMolecule->v(2);
						_localDirYVelocity2[unID] +=
							(temporaryMolecule->v(1) - _directedVelocityOld) *
							(temporaryMolecule->v(1) -
							 _directedVelocityOld);  // velocity2 without directed velocity of Droplet
						_localXyzVi[unID][0] += temporaryMolecule->Vi(0);
						_localXyzVi[unID][1] += temporaryMolecule->Vi(1);
						_localXyzVi[unID][2] += temporaryMolecule->Vi(2);
						_particles += 1.;
					}
				}
			}
		}

		if (synsimstep == _outputFrequency) {
			// PARALLELIZATION
			double unIDs = _rIncrements * _hIncrements * _phiIncrements;
			domainDecomp->collCommInit(unIDs * 11);

			for (int uID = 1; uID <= (_rIncrements * _hIncrements * _phiIncrements); ++uID) {
				domainDecomp->collCommAppendDouble(_localnumberOfParticles[uID]);
				domainDecomp->collCommAppendDouble(_localXyzVelocities[uID][0]);
				domainDecomp->collCommAppendDouble(_localXyzVelocities[uID][1]);
				domainDecomp->collCommAppendDouble(_localXyzVelocities[uID][2]);
				domainDecomp->collCommAppendDouble(_localXyzVelocities2[uID][0]);
				domainDecomp->collCommAppendDouble(_localXyzVelocities2[uID][1]);
				domainDecomp->collCommAppendDouble(_localXyzVelocities2[uID][2]);
				domainDecomp->collCommAppendDouble(_localDirYVelocity2[uID]);
				domainDecomp->collCommAppendDouble(_localXyzVi[uID][0]);
				domainDecomp->collCommAppendDouble(_localXyzVi[uID][1]);
				domainDecomp->collCommAppendDouble(_localXyzVi[uID][2]);
			}
			domainDecomp->collCommAllreduceSum();
			for (int uID = 1; uID <= (_rIncrements * _hIncrements * _phiIncrements); ++uID) {
				_globalnumberOfParticles[uID] = domainDecomp->collCommGetDouble();
				_globalXyzVelocities[uID][0] = domainDecomp->collCommGetDouble();
				_globalXyzVelocities[uID][1] = domainDecomp->collCommGetDouble();
				_globalXyzVelocities[uID][2] = domainDecomp->collCommGetDouble();
				_globalXyzVelocities2[uID][0] = domainDecomp->collCommGetDouble();
				_globalXyzVelocities2[uID][1] = domainDecomp->collCommGetDouble();
				_globalXyzVelocities2[uID][2] = domainDecomp->collCommGetDouble();
				_globalDirYVelocity2[uID] = domainDecomp->collCommGetDouble();
				_globalXyzVi[uID][0] = domainDecomp->collCommGetDouble();
				_globalXyzVi[uID][1] = domainDecomp->collCommGetDouble();
				_globalXyzVi[uID][2] = domainDecomp->collCommGetDouble();
			}
			domainDecomp->collCommFinalize();

			// VELOCITYS AND PARTICLE NUMBERS OF EVERY BOX IN THE DROPLET
			double velocDroplet[3];
			velocDroplet[0] = 0.;
			velocDroplet[1] = 0.;
			velocDroplet[2] = 0.;
			double velocDroplet2[3];
			velocDroplet2[0] = 0.;
			velocDroplet2[1] = 0.;
			velocDroplet2[2] = 0.;
			double viDroplet[3];
			viDroplet[0] = 0.;
			viDroplet[1] = 0.;
			viDroplet[2] = 0.;
			double viGas[3];
			viGas[0] = 0.;
			viGas[1] = 0.;
			viGas[2] = 0.;
			double particleDropletTrue = 0.;
			double particleDropletFalse = 0.;
			double sumBoxesLiquid = 0.;
			double sumBoxesVapor = 0.;
			double sumVirialDroplet = 0.;
			double sumVirialGas = 0.;
			double sumRohLiq = 0.;
			double UnIDMin = _rIncrements * _hIncrements * _phiIncrements;
			double UnIDMax = 0.;
			double hUnMin = 0.;
			double hUnMax = 0.;
			double hMin = 0.;
			double hMax = 0.;
			double r2Max = 0.;
			double rUnMax = 0.;
			double sumBoxNotGas = 0.;
			// DEFINE ROHCUTLIQ
			unsigned ModuloStart = _iterationsSinceStart % _outputFrequency;

			_iterationsSinceStart = 0;
			if (_firstTime == true) {
				rohCutLiqActual = _rohCutLiq * _percent;

				if (ModuloStart == 0) {
					_firstTime = false;
				}
			} else {
				rohCutLiqActual = _rohCutLiqNew * _percent;
			}
			// COUNT ALL VELOCITYS AND PARTICLESNUMBERS OVER EVERY BOX AND EVERY SYNSIMSTEP TOGETHER
			for (unsigned unID = 1; unID <= (_rIncrements * _hIncrements * _phiIncrements); ++unID) {
				// CALCULATE DESITYS OUT OF THE PARTICLENUMBERS AND THE BOX-DIMENSIONS
				_densityBox[unID] = (_globalnumberOfParticles[unID] / _outputFrequency) / _volumebox;
				_temperatureBox[unID] =
					2. / (3. * _globalnumberOfParticles[unID]) * 0.5 *
					(_globalXyzVelocities2[unID][0] + _globalXyzVelocities2[unID][1] + _globalXyzVelocities2[unID][2]);
				_temperatureBoxXZ[unID] = 2. / (2. * _globalnumberOfParticles[unID]) * 0.5 *
										  (_globalXyzVelocities2[unID][0] + _globalXyzVelocities2[unID][2]);
				_EkinBox[unID] = 0.5 * (_globalXyzVelocities2[unID][0] + _globalXyzVelocities2[unID][1] +
										_globalXyzVelocities2[unID][2]);
				_virialBox[unID] =
					_temperatureBox[unID] * (_globalnumberOfParticles[unID] / _outputFrequency) / _volumebox +
					(1. / _volumebox) * (_globalXyzVi[unID][0] + _globalXyzVi[unID][1] + _globalXyzVi[unID][2]) /
						(3 * _outputFrequency);
				// pressurexDROPLET = TxyzDROPLET * (particleDropletTrue/ _outputFrequency) / (_volumebox *
				// sumBoxesLiquid) + 1. / (_volumebox * sumBoxesLiquid) * (viDroplet[0] / _outputFrequency);
				// COUNT ALL VALUES IF BOX IS IN THE DROPLET
				if ((_densityBox[unID] > rohCutLiqActual) && (_permissibleRange[unID] == 1)) {
					_temperatureBox[unID] =
						2. / (3. * _globalnumberOfParticles[unID]) * 0.5 *
						(_globalXyzVelocities2[unID][0] + _globalDirYVelocity2[unID] + _globalXyzVelocities2[unID][2]);
					_EkinBox[unID] = 0.5 * (_globalXyzVelocities2[unID][0] + _globalDirYVelocity2[unID] +
											_globalXyzVelocities2[unID][2]);
					_virialBox[unID] =
						_temperatureBox[unID] * (_globalnumberOfParticles[unID] / _outputFrequency) / _volumebox +
						(1. / _volumebox) * (_globalXyzVi[unID][0] + _globalXyzVi[unID][1] + _globalXyzVi[unID][2]) /
							(3 * _outputFrequency);
					particleDropletTrue += (_globalnumberOfParticles[unID]);
					velocDroplet[0] += ((_globalXyzVelocities[unID][0]));
					velocDroplet[1] += ((_globalXyzVelocities[unID][1]));
					velocDroplet[2] += ((_globalXyzVelocities[unID][2]));
					viDroplet[0] += _globalXyzVi[unID][0];
					viDroplet[1] += _globalXyzVi[unID][1];
					viDroplet[2] += _globalXyzVi[unID][2];
					_xyzEkinDroplet[_outputFrequency] +=
						0.5 *
						(_globalXyzVelocities2[unID][0] + _globalDirYVelocity2[unID] + _globalXyzVelocities2[unID][2]);
					_xzEkinDroplet[_outputFrequency] +=
						0.5 * (_globalXyzVelocities2[unID][0] + _globalXyzVelocities2[unID][2]);
					sumRohLiq += _densityBox[unID];
					sumBoxesLiquid += 1.;
					if (unID < UnIDMin) {
						UnIDMin = unID;
					}
					if (unID > UnIDMax) {
						UnIDMax = unID;
					}
				}
			}
			// CUT OFF CYLINDER WITH DROPLET
			hUnMin = floor(UnIDMin / _rIncrements);
			hUnMax = floor(UnIDMax / _rIncrements);
			hMin = hUnMin / _universalInvProfileUnit[1];
			hMax = hUnMax / _universalInvProfileUnit[1];

			hMin = hMin - 5.;
			hMax = hMax + 5.;

			hUnMin = floor(hMin * _universalInvProfileUnit[1]);
			hUnMax = floor(hMax * _universalInvProfileUnit[1]);

			// GET RADIUS OF DROPLET WITH HIGHEST AND LOWEST BOX FROM DROPLET
			r2Max = (hMax - hMin) / 2;
			r2Max = r2Max * r2Max;
			rUnMax = ceil(r2Max * _universalInvProfileUnit[0]);

			for (int h = hUnMin; h <= hUnMax; h++) {
				for (int r = 0; r <= rUnMax; r++) {
					_permissibleRange[h * _rIncrements + r] = 0.;
					sumBoxNotGas += 1;
				}
			}

			// COUNT ALL VALUES IF BOX IS IN THE GAS PHASE
			for (unsigned unID = 1; unID <= (_rIncrements * _hIncrements * _phiIncrements); ++unID) {
				if (_permissibleRange[unID] == 1) {
					particleDropletFalse += (_globalnumberOfParticles[unID]);
					viGas[0] += _globalXyzVi[unID][0];
					viGas[1] += _globalXyzVi[unID][1];
					viGas[2] += _globalXyzVi[unID][2];
					_xyzEkinGas[_outputFrequency] +=
						0.5 * (_globalXyzVelocities2[unID][0] + _globalXyzVelocities2[unID][1] +
							   _globalXyzVelocities2[unID][2]);
					_xzEkinGas[_outputFrequency] +=
						0.5 * (_globalXyzVelocities2[unID][0] + _globalXyzVelocities2[unID][2]);
					sumBoxesVapor += 1.;
				}
				// SET ALL LOCALS TO ZERO
				_localnumberOfParticles[unID] = 0.;
				_localXyzVelocities[unID][0] = 0.;
				_localXyzVelocities[unID][1] = 0.;
				_localXyzVelocities[unID][2] = 0.;
				_localXyzVelocities2[unID][0] = 0.;
				_localXyzVelocities2[unID][1] = 0.;
				_localXyzVelocities2[unID][2] = 0.;
				_localDirYVelocity2[unID] = 0.;
				_localXyzVi[unID][0] = 0.;
				_localXyzVi[unID][1] = 0.;
				_localXyzVi[unID][2] = 0.;
			}
			// RESET PERMISSIBLE RANGE
			for (int h = hUnMin; h <= hUnMax; h++) {
				for (int r = 0; r <= rUnMax; r++) {
					_permissibleRange[h * _rIncrements + r] = 1.;
				}
			}
			// GET NEW DENSITY CUTOFF FOR LIQUIDPHASE
			_rohCutLiqNew = (particleDropletTrue / _outputFrequency) / (sumBoxesLiquid * _volumebox);
			// GET THE TOTAL-VELOCITYS OF THE BOXES IN THE DROPLET FROM EVERY DIRECTION
			velocDroplet[0] = velocDroplet[0] / particleDropletTrue;
			velocDroplet[1] = velocDroplet[1] / particleDropletTrue;
			velocDroplet[2] = velocDroplet[2] / particleDropletTrue;
			_directedVelocityOld = velocDroplet[1];

			// SAFE ACTUAL SIMSTEP AND VELOCITYS FOR OUTPUT
			_simstepArray[_counter] = simstep;
			_velocDroplet[_counter][0] = velocDroplet[0];
			_velocDroplet[_counter][1] = velocDroplet[1];
			_velocDroplet[_counter][2] = velocDroplet[2];

			// CALCULATE GLOBAL TEMPERATURE --> Hierf�r Ekin �ber alle Bins aufsummieren und durch particelnumberdroplet
			// teilen
			double TxyzDROPLET = 2. / (3. * particleDropletTrue) * _xyzEkinDroplet[_outputFrequency];
			double TxzDROPLET = 2. / (2. * particleDropletTrue) * _xzEkinDroplet[_outputFrequency];
			double TxyzGAS = 2. / (3. * particleDropletFalse) * _xyzEkinGas[_outputFrequency];
			double TxzGAS = 2. / (2. * particleDropletFalse) * _xzEkinGas[_outputFrequency];

			// CALCULATE PRESSURE
			double pressurexDROPLET =
				TxyzDROPLET * (particleDropletTrue / _outputFrequency) / (_volumebox * sumBoxesLiquid) +
				1. / (_volumebox * sumBoxesLiquid) * (viDroplet[0] / _outputFrequency);
			double pressureyDROPLET =
				TxyzDROPLET * (particleDropletTrue / _outputFrequency) / (_volumebox * sumBoxesLiquid) +
				1. / (_volumebox * sumBoxesLiquid) * (viDroplet[1] / _outputFrequency);
			double pressurezDROPLET =
				TxyzDROPLET * (particleDropletTrue / _outputFrequency) / (_volumebox * sumBoxesLiquid) +
				1. / (_volumebox * sumBoxesLiquid) * (viDroplet[2] / _outputFrequency);
			double pressurexGAS = TxyzGAS * (particleDropletFalse / _outputFrequency) / (_volumebox * sumBoxesVapor) +
								  1. / (_volumebox * sumBoxesVapor) * (viGas[0] / _outputFrequency);
			double pressureyGAS = TxyzGAS * (particleDropletFalse / _outputFrequency) / (_volumebox * sumBoxesVapor) +
								  1. / (_volumebox * sumBoxesVapor) * (viGas[1] / _outputFrequency);
			double pressurezGAS = TxyzGAS * (particleDropletFalse / _outputFrequency) / (_volumebox * sumBoxesVapor) +
								  1. / (_volumebox * sumBoxesVapor) * (viGas[2] / _outputFrequency);

			_particles = 0.;
			_counter += 1;

			int mpi_rank = domainDecomp->getRank();
			if (mpi_rank == 0) {
				// WRITE ALL SIMSTEPS, gerichtete Geschwindigkeit, dichteGas, dichteLiq, druckGas, druckLiq, TGas, TLiq,
				// EkinGas und EkinLiq IN AN OUTPUT FILE
				_DPMGlobalStream.open("Global_output_DPM_MK.txt", std::ios::app);
				if (simstep / _outputFrequency > 1) {
					_DPMGlobalStream << simstep << " \t\t" << _directedVelocityOld << " \t\t "
									 << ((particleDropletFalse / _outputFrequency) / (_volumebox * sumBoxesVapor))
									 << " \t\t " << _rohCutLiqNew << " \t\t "
									 << ((pressurexGAS + pressureyGAS + pressurezGAS) / 3) << " \t\t "
									 << ((pressurexDROPLET + pressureyDROPLET + pressurezDROPLET) / 3) << " \t\t "
									 << TxyzGAS << " \t\t " << TxyzDROPLET << " \t\t "
									 << _xyzEkinGas[_outputFrequency] / (3 * particleDropletFalse) << " \t\t "
									 << _xyzEkinDroplet[_outputFrequency] / (3 * particleDropletTrue) << " \t\t "
									 << TxzGAS << " \t\t " << TxzDROPLET << " \t\t "
									 << _xzEkinGas[_outputFrequency] / (2 * particleDropletFalse) << " \t\t "
									 << _xzEkinDroplet[_outputFrequency] / (2 * particleDropletTrue) << endl;
				}
				_DPMGlobalStream.close();

				// 2D DENSITY OUTPUT
				DPMStreamDensity.open("drop_MK_DirectedPM_" + to_string(simstep) + ".NDpr", std::ios::out);
				DPMStreamDensity.precision(6);
				// Write Header
				DPMStreamDensity << "//Segment volume: " << _volumebox
								 << "\n//Accumulated data sets: " << _outputFrequency
								 << "\n//Local profile of the number density. Output file generated by the "
									"\"DensityProfile\" method, plugins/profiles. \n";
				DPMStreamDensity << "//local density profile: Y - Z || X-projection\n";
				DPMStreamDensity << "// \t dr \t dh \t dphi \n";
				DPMStreamDensity << "\t" << 1 / _universalInvProfileUnit[0] << "\t" << 1 / _universalInvProfileUnit[1]
								 << "\t" << 1 / _universalInvProfileUnit[2] << "\n";
				DPMStreamDensity << "0 \t";
				for (unsigned r = 0; r < _rIncrements; r++) {
					DPMStreamDensity << 0.5 * (sqrt(r + 1) + sqrt(r)) / sqrt(_universalInvProfileUnit[0])
									 << " \t";  // Eintragen der R-Koordinate in Header korrigiert
				}
				//_DPMStreamDensity << "END R2";
				DPMStreamDensity << "\n";
				// Y - axis label
				for (unsigned h = 0; h < _hIncrements; h++) {
					double hval = (h + 0.5) / _universalInvProfileUnit[1];
					DPMStreamDensity << hval << "  \t";
					for (unsigned phiID = 0; phiID < _phiIncrements; phiID++) {
						for (unsigned r = 0; r < _rIncrements; r++) {
							auto ID = (long)(h * _rIncrements * _phiIncrements + r * _phiIncrements + phiID);
							DPMStreamDensity << _densityBox[ID] << "\t";
						}
					}
					DPMStreamDensity << "\n";
				}
				DPMStreamDensity.close();

				// 2D Temperature OUTPUT
				DPMStreamTemperature.open("drop_MK_DirectedPM_" + to_string(simstep) + ".Temppr", std::ios::out);
				DPMStreamTemperature.precision(6);
				// Write Header
				DPMStreamTemperature << "//Segment volume: " << _volumebox
									 << "\n//Accumulated data sets: " << _outputFrequency
									 << "\n//Local profile of the number temperature. \n";
				DPMStreamTemperature << "//Temperature expressed by 2Ekin/#DOF\n";
				DPMStreamTemperature << "//local temperature profile: Y - Z || X-projection\n";
				DPMStreamTemperature << "// \t dr \t dh \t dphi \n";
				DPMStreamTemperature << "\t" << 1 / _universalInvProfileUnit[0] << "\t"
									 << 1 / _universalInvProfileUnit[1] << "\t" << 1 / _universalInvProfileUnit[2]
									 << "\n";
				DPMStreamTemperature << "0 \t";
				for (unsigned r = 0; r < _rIncrements; r++) {
					DPMStreamTemperature << 0.5 * (sqrt(r + 1) + sqrt(r)) / sqrt(_universalInvProfileUnit[0])
										 << " \t";  // Eintragen der R-Koordinate in Header korrigiert
				}
				DPMStreamTemperature << "\n";
				// Y - axis label
				for (unsigned h = 0; h < _hIncrements; h++) {
					double hval = (h + 0.5) / _universalInvProfileUnit[1];
					DPMStreamTemperature << hval << "  \t";
					for (unsigned phiID = 0; phiID < _phiIncrements; phiID++) {
						for (unsigned r = 0; r < _rIncrements; r++) {
							auto ID = (long)(h * _rIncrements * _phiIncrements + r * _phiIncrements + phiID);
							if (isnan(_temperatureBox[ID])) {
								_temperatureBox[ID] = 0.;
							}
							DPMStreamTemperature << _temperatureBox[ID] << "\t";
						}
					}
					DPMStreamTemperature << "\n";
				}
				DPMStreamTemperature.close();

				// 2D TemperatureXZ OUTPUT
				DPMStreamTemperatureXZ.open("drop_MK_DirectedPM_" + to_string(simstep) + ".TempprXZ", std::ios::out);
				DPMStreamTemperatureXZ.precision(6);
				// Write Header
				DPMStreamTemperatureXZ << "//Segment volume: " << _volumebox
									   << "\n//Accumulated data sets: " << _outputFrequency
									   << "\n//Local profile of the number temperature. \n";
				DPMStreamTemperatureXZ
					<< "//Temperature expressed by 2Ekin/#DOF\t MIT XZ-WERTEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
				DPMStreamTemperatureXZ << "//local temperature profile: Y - Z || X-projection\n";
				DPMStreamTemperatureXZ << "// \t dr \t dh \t dphi \n";
				DPMStreamTemperatureXZ << "\t" << 1 / _universalInvProfileUnit[0] << "\t"
									   << 1 / _universalInvProfileUnit[1] << "\t" << 1 / _universalInvProfileUnit[2]
									   << "\n";
				DPMStreamTemperatureXZ << "0 \t";
				for (unsigned r = 0; r < _rIncrements; r++) {
					DPMStreamTemperatureXZ << 0.5 * (sqrt(r + 1) + sqrt(r)) / sqrt(_universalInvProfileUnit[0])
										   << " \t";  // Eintragen der R-Koordinate in Header korrigiert
				}
				DPMStreamTemperatureXZ << "\n";
				// Y - axis label
				for (unsigned h = 0; h < _hIncrements; h++) {
					double hval = (h + 0.5) / _universalInvProfileUnit[1];
					DPMStreamTemperatureXZ << hval << "  \t";
					for (unsigned phiID = 0; phiID < _phiIncrements; phiID++) {
						for (unsigned r = 0; r < _rIncrements; r++) {
							auto ID = (long)(h * _rIncrements * _phiIncrements + r * _phiIncrements + phiID);
							if (isnan(_temperatureBoxXZ[ID])) {
								_temperatureBoxXZ[ID] = 0.;
							}
							DPMStreamTemperatureXZ << _temperatureBoxXZ[ID] << "\t";
						}
					}
					DPMStreamTemperatureXZ << "\n";
				}
				DPMStreamTemperatureXZ.close();

				// 2D Ekin OUTPUT
				DPMStreamEkin.open("drop_MK_DirectedPM_" + to_string(simstep) + ".Ekin", std::ios::out);
				DPMStreamEkin.precision(6);
				// Write Header
				DPMStreamEkin << "//Segment volume: " << _volumebox << "\n//Accumulated data sets: " << _outputFrequency
							  << "\n//Local profile of the number Ekin. \n";
				DPMStreamEkin << "//local Ekin profile: Y - Z || X-projection\n";
				DPMStreamEkin << "// \t dr \t dh \t dphi \n";
				DPMStreamEkin << "\t" << 1 / _universalInvProfileUnit[0] << "\t" << 1 / _universalInvProfileUnit[1]
							  << "\t" << 1 / _universalInvProfileUnit[2] << "\n";
				DPMStreamEkin << "0 \t";
				for (unsigned r = 0; r < _rIncrements; r++) {
					DPMStreamEkin << 0.5 * (sqrt(r + 1) + sqrt(r)) / sqrt(_universalInvProfileUnit[0])
								  << " \t";  // Eintragen der R-Koordinate in Header korrigiert
				}
				DPMStreamEkin << "\n";
				// Y - axis label
				for (unsigned h = 0; h < _hIncrements; h++) {
					double hval = (h + 0.5) / _universalInvProfileUnit[1];
					DPMStreamEkin << hval << "  \t";
					for (unsigned phiID = 0; phiID < _phiIncrements; phiID++) {
						for (unsigned r = 0; r < _rIncrements; r++) {
							auto ID = (long)(h * _rIncrements * _phiIncrements + r * _phiIncrements + phiID);
							DPMStreamEkin << _EkinBox[ID] << "\t";
						}
					}
					DPMStreamEkin << "\n";
				}
				DPMStreamEkin.close();

				// 2D VIRIAL OUTPUT
				DPMStreamVirial.open("drop_MK_DirectedPM_" + to_string(simstep) + ".Vipr", std::ios::out);
				DPMStreamVirial.precision(6);
				// Write Header
				DPMStreamVirial << "//Segment volume: " << _volumebox
								<< "\n//Accumulated data sets: " << _outputFrequency
								<< "\n//Local profile of the partial pressures p = (px+py+pz)/3). \n";
				DPMStreamVirial << "// \t dr \t dh \t dphi \n";
				DPMStreamVirial << "\t" << 1 / _universalInvProfileUnit[0] << "\t" << 1 / _universalInvProfileUnit[1]
								<< "\t" << 1 / _universalInvProfileUnit[2] << "\n";
				DPMStreamVirial << "0 \t";
				for (unsigned r = 0; r < _rIncrements; r++) {
					DPMStreamVirial << 0.5 * (sqrt(r + 1) + sqrt(r)) / sqrt(_universalInvProfileUnit[0])
									<< " \t";  // Eintragen der R-Koordinate in Header korrigiert
				}
				DPMStreamVirial << "\n";
				// Y - axis label
				for (unsigned h = 0; h < _hIncrements; h++) {
					double hval = (h + 0.5) / _universalInvProfileUnit[1];
					DPMStreamVirial << hval << "  \t";
					for (unsigned phiID = 0; phiID < _phiIncrements; phiID++) {
						for (unsigned r = 0; r < _rIncrements; r++) {
							auto ID = (long)(h * _rIncrements * _phiIncrements + r * _phiIncrements + phiID);
							if (isnan(_virialBox[ID])) {
								_virialBox[ID] = 0.;
							}
							DPMStreamVirial << _virialBox[ID] << "\t";
						}
					}
					DPMStreamVirial << "\n";
				}
				DPMStreamVirial.close();
			}
		}
	}
}

void DirectedPM::endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
						 unsigned long simstep) {
	if (_enabled) {
		if (simstep > 0) {
			// CALCULATE SYNSIMSTEP FOR MATRIX
			double synsimstep;
			if ((simstep) % _outputFrequency != 0) {
				synsimstep = (simstep % _outputFrequency);
			} else {
				synsimstep = _outputFrequency;
			}

			// synsimstep als globale Variable deklarieren, dann ist Berechnung oben dr�ber unn�tig
			// WRITE AN OUTPUT FILE TO CHECK THE DENSITYS AND THE POSITION OF THE DROP IN THE MIDDLE X-LAYER (ONLY FOR
			// TESTS)
			if (synsimstep == _outputFrequency) {
				for (int j = 1; j <= (_rIncrements * _hIncrements * _phiIncrements); j++) {
					for (int i = 0; i < 3; i++) {
						_globalXyzVelocities[j][i] = 0.;
						_globalXyzVelocities2[j][i] = 0.;
						_globalXyzVi[j][i] = 0.;
					}
					_globalDirYVelocity2[j] = 0.;
					_globalnumberOfParticles[j] = 0.;
					_densityBox[j] = 0.;
					_temperatureBox[j] = 0.;
					_EkinBox[j] = 0.;
					_virialBox[j] = 0.;
				}
			}
		}
	}
}
