/*
 * DirectedPM.h
 *
 *  Created on: 25 august 2020
 *      Author: Koch
 */

#ifndef MARDYN_TRUNK_DIRECTEDPM_H
#define MARDYN_TRUNK_DIRECTEDPM_H

//class DirectedPMTest;
#include "PluginBase.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"

/** @brief
* Plugin: can be enabled via config.xml <br>
*
* Calculates velocity of a moving droplet Also determines density, temperature, and pressure<br>
* of the droplet and the surrounding vapor phase. For this purpose the phase boundary is determined<br>
* via the local density in a cylindrical coordinate system using a binning scheme<br>
* Temperature is then calculated by subtracting the overall droplet velocity from the individual particle velocities<br>
* to yield only the thermal velocities. Alternatively, temperature is also calculated by only considering velocities in x-<br>
* and z-direction. The droplet is assumed to move in y-direction.<br>
* \code{.xml}
* <plugin>
*           <Component>1</Component>
*			<rIncrements>50</rIncrements>
*           <hIncrements>50</hIncrements>
*	        <phiIncrements>1</phiIncrements>
*	        <rohCutLiq>0.5</rohCutLiq>
*	        <maxDeviation>5</maxDeviation>
*	        <heightWall>5</heightWall>
*	        <heightMembrane>170</heightMembrane>
*           <outputFrequency>1000</outputFrequency>
*			
* </plugin>
* \endcode
*/

class DirectedPM : public PluginBase{

private:
    //friend DirectedPMTest;
	
	bool _enabled = true;

	int _measureinterval =1 ;
	double _directedVelocityOld;

    unsigned _outputFrequency;
    unsigned _component;
	double _rIncrements;
	double _hIncrements;
	double _phiIncrements;
    double _boxLength[3];
    double _binSize[3];
    double _rohCutLiq;
    double _heightWall;
    double _heightMembrane;
    double _percent;
    int _yLowestBox;
    int _yHighestBox;
    double _yLow;
    double _yHigh;
    double _volumebox;
    int _counter;
    
    double _particles;
    double _universalInvProfileUnit[3];
    double _universalCentre[3];
    double _minXZ;
    double _R2max;

    bool _firstTime;
    double _rohCutLiqNew;

    unsigned _iterationsSinceStart;

    std::map<unsigned, std::map<unsigned, double>> _xyzEkin;
    std::map<unsigned, std::map<unsigned, double>> _xzEkin;

    std::map<unsigned, double> _localnumberOfParticles;
    std::map<unsigned, std::map<unsigned, double>> _localXyzVi;
    std::map<unsigned, std::map<unsigned, double>> _localXyzVelocities;
    std::map<unsigned, std::map<unsigned, double>> _localXyzVelocities2;
    std::map<unsigned, double> _localDirYVelocity2;    
    
    std::map<unsigned, double> _globalnumberOfParticles;
    std::map<unsigned, std::map<unsigned, double>> _globalXyzVi;
    std::map<unsigned, std::map<unsigned, double>> _globalXyzVelocities;
    std::map<unsigned, std::map<unsigned, double>> _globalXyzVelocities2;
    std::map<unsigned, double> _globalDirYVelocity2;

    std::map<unsigned, double> _xyzEkinDroplet;
    std::map<unsigned, double> _xyzEkinGas;
    std::map<unsigned, double> _xzEkinDroplet;
    std::map<unsigned, double> _xzEkinGas;

    std::map<unsigned, double> _permissibleRange;

    std::map<unsigned, double> _densityBox;
    std::map<unsigned, double> _temperatureBox;
    std::map<unsigned, double> _temperatureBoxXZ;
    std::map<unsigned, double> _EkinBox;
    std::map<unsigned, double> _virialBox; 

    std::ofstream  DPMStreamDensity;
    std::ofstream  DPMStreamVirial;
    std::ofstream  DPMStreamTemperature;
    std::ofstream  DPMStreamTemperatureXZ;
    std::ofstream  DPMStreamEkin;

    std::map<unsigned long, double> _simstepArray;
    std::map<unsigned, std::map<unsigned, double>> _velocDroplet;
    std::ofstream _DPMGlobalStream;

public:


    void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override {
        global_log->debug() << "DirectPM enabled" << std::endl;
        _firstTime = true;
        for (unsigned d = 0; d < 3; d++) {
            _boxLength[d] = domain->getGlobalLength(d);
        }
        _iterationsSinceStart = 0;
        _binSize[0] = _boxLength[0] / _rIncrements;
        _binSize[1] = _boxLength[1] / _hIncrements;
        _binSize[2] = _boxLength[2] / _phiIncrements;
        _yLowestBox = ceil(_heightWall / _binSize[1]);
        _yHighestBox = floor(_heightMembrane / _binSize[1]);
        _yLow = _yLowestBox * _binSize[1];
        _yHigh = _yHighestBox * _binSize[1];

        _directedVelocityOld = 0.;
        _counter = 0;
        _particles = 0.;
        //GET CYLINDRICAL COORDINATES
        _minXZ = _boxLength[0];
        if (_boxLength[2] < _minXZ) {
            _minXZ = _boxLength[2];
        }
        _R2max = 0.24 * _minXZ * _minXZ;
        _universalInvProfileUnit[0] = _rIncrements / (_R2max);                    // delta_R2
        _universalInvProfileUnit[1] = _hIncrements / _boxLength[1];              // delta_H
        _universalInvProfileUnit[2] = _phiIncrements / (2 * M_PI);                 // delta_Phi
        _universalCentre[0] = 0.5 * _boxLength[0];
        _universalCentre[1] = 0.;
        _universalCentre[2] = 0.5 * _boxLength[2];
        _volumebox = M_PI / (_universalInvProfileUnit[0] * _universalInvProfileUnit[1] * _phiIncrements);

        //SET PERMISSIBLE RANGE OF EVERY BOX TO ONE
        for (int l = 0; l <= ((_rIncrements * _hIncrements * _phiIncrements)); l++) {
            _permissibleRange[l] = 1.;
        }
        // SET PERMISSIBLE RANGE OF EVERY BOX BELOW ADSORPTION LAYER TO ZERO
        for (int i = 0; i <= _yLowestBox; i++) {
            for (int j = 0; j <= _rIncrements; j++) {
                for (int k = 0; k <= _phiIncrements; k++) {
                    _permissibleRange[(i * _rIncrements * _phiIncrements) + (j * _phiIncrements) + k] = 0.;
                }
            }
        }
        // SET PERMISSIBLE RANGE OF EVERY BOX ABOVE MEMBRANE TO ZERO
        for (int m = _yHighestBox; m <= (_hIncrements); m++) {
            for (int n = 0; n <= _rIncrements; n++) {
                for (int o = 0; o <= _phiIncrements; o++) {
                    _permissibleRange[(m * _rIncrements * _phiIncrements) + (n * _phiIncrements) + o] = 0.;
                }
            }
        }
        // WRITE ALL SIMSTEPS, gerichtete Geschwindigkeit, dichteGas, dichteLiq, druckGas, druckLiq, TGas, TLiq, EkinGas und EkinLiq IN AN OUTPUT FILE
        _DPMGlobalStream.open("Global_output_DPM_MK.txt", std::ios::out);
        _DPMGlobalStream << "Ausgabe der globalen Größen gerichtete Geschwindigkeit, dichteGas, dichteLiq, druckGas, druckLiq, TGas, TLiq, EkinxyzGas, EkinxyzLiq, TGasXZ, TLiqXZ,EkinxzGas und EkinxzLiq," << endl;
        _DPMGlobalStream << endl;
        _DPMGlobalStream << "Timestept \t\t  gerichtete Geschw. \t\t  dichteGas \t\t  dichteLiq \t\t  druckGas \t\t  druckLiq \t\t  TGas \t\t  TLiq \t\t  EkinxyzGas \t\t  EkinxyzLiq\t\t  TGasXZ\t\t TLiqXZ\t\t  EkinxzGas \t\t  EkinxzLiq\t\t" << endl;
        _DPMGlobalStream.close();
    }
	void readXML (XMLfileUnits& xmlconfig) override;

    void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;

    void endStep(
            ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
            Domain *domain, unsigned long simstep) override;


    void finish(ParticleContainer *particleContainer,
                DomainDecompBase *domainDecomp, Domain *domain) override {};

    std::string getPluginName()override {return std::string("DirectedPM");}

    static PluginBase* createInstance(){return new DirectedPM();}

};


#endif //MARDYN_TRUNK_DIRECTEDPM_H
