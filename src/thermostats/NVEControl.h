/*
 * NVEControl.h
 *
 *  Created on: 05.07.2024
 *      Author: jniemann
 * 
 * maintaines total energy to match initial total energy level.
 * component-wise not supported. only manipulates vTrans, not rotation. (suited e.g. for single-centered LJ)
 * 
 * 
 */

#pragma once 

#include "particleContainer/ParticleContainer.h"
	
	
	/** @brief Read in XML configuration for NVEControl and all its included objects.
	 *
	 * The following XML object structure is handled by this method:
	 * \code{.xml}
		<thermostats>
			<thermostat type="NVEControl">
				<control>
					<start>UNSIGNED_LONG</start>           <!-- Timestep turning thermostat ON -->
					<frequency>UNSIGNED_LONG</frequency>   <!-- Thermosatt is active every <frequency>-th time step -->
					<stop>UNSIGNED_LONG</stop>             <!-- Timestep turning thermostat OFF -->
					<!-- SO FAR NOT DOING ANYTHING <writefreq>5000</writefreq>         <!-- SO FAR NOT DOING ANYTHING Log thermostat scaling factors betaTrans --> 
					<fileprefix>betalog</fileprefix>    <!-- Prefix of log file -->
				</control>
			</thermostat>
		</thermostats>
	   \endcode
	 */



class NVEControl {
public:
	NVEControl();
	~NVEControl();
	void setBetaTrans(double beta);
	double getBetaTrans() { return _globalBetaTrans; }
	void apply(ParticleContainer *moleculeContainer);
	void readXML(XMLfileUnits& xmlconfig);

private:
	double _globalBetaTrans;

	unsigned long int _Nglobal;
	double _Utarget;

	unsigned long _startTime;
	unsigned long _endTime;
	int _frequency;
	int _writeFrequency;
	std::string _fileprefix;
};