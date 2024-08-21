/*
 * SphericalSampling.cpp
 *
 *  Created on: Aug 2024
 *      Author: JakNiem
 */


#pragma once

#include "LongRangeCorrection.h"
#include "Simulation.h"



#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class Spherical:public LongRangeCorrection {
public:
	Spherical(double cutoffT, double cutoffLJ, Domain* domain,  DomainDecompBase* domainDecomposition, ParticleContainer* particleContainer, Simulation* simulation);
	virtual ~Spherical();

	virtual void init();
	virtual void readXML(XMLfileUnits& xmlconfig);
	virtual void calculateLongRange();
	virtual void writeProfiles(DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep);


private:
	template<typename T>
	void resizeExactly(std::vector<T>& v, unsigned int numElements) const {
		v.reserve(numElements);
		v.resize(numElements);
	}

	template<typename T>
	void initCommVarVector(CommVar<std::vector<T>>& commVar, unsigned int len) const {
		commVar.local.resize(len);
		commVar.global.resize(len);
	}
	template<typename T>
	void initCommVarVector(CommVar<std::vector<T>>& commVar) const {
		commVar.local.resize(_lenVector);
		commVar.global.resize(_lenVector);
		std::fill(commVar.global.begin(), commVar.global.end(), static_cast<T>(0.)); 
		std::fill(commVar.local.begin(), commVar.local.end(), static_cast<T>(0.)); 
	}



	void calculateBulkDensities();
	void calculateVirtualDensity();

    double _cutoffLJ;
	Domain* _domain;
	DomainDecompBase* _domainDecomposition;
	ParticleContainer* _particleContainer;

    // Control: general
    unsigned int _nShells {60};
    unsigned int _lenVector {60+1};
    bool _isBubble {false};
	bool _disableLRC {false};
    double _T {std::nan("0")};      // no default --> error message if not set
	unsigned int _calcFreq {100};
	unsigned int _writeFreq {1000}; //has to be a multiple of calculation Frequency
    
    std::string _outputPrefix {"sphericalLRC_"};


    // Auxiliary variables
    std::array<double, 3> _globalBoxLength;
    std::array<double, 3> _globalCenter;
	double _distMax;
	// double _distMaxSquared;
	// double _shellWidth;
    unsigned long _globalNumMols;


	//SHELL PROPERTIES:
    std::vector<double> _shellVolume;             
    std::vector<double> _shellLowerBound; //_shellUpperBound[i] = _shellLowerBound[i+1]            

	//DENSITY VECTORS:
	std::vector<unsigned long> _numMolecules_accum;  //for correction caclulation (denstity)           // Number of molecules in bin
	std::vector<double> _density_avg; 	//avg over time (TODO: so far: all time average -- should be changed in the future)
	std::vector<double> _density_avg_fitted; 
	std::vector<double> _virtual_density;   // _density_avg_fitted - homogenous density

	//SHELLWISE CORRECTION TERMS:
	std::vector<double> _FcorrectionShell;
	std::vector<double> _UcorrectionShell;
	std::vector<double> _VirNcorrectionShell;
	std::vector<double> _VirTcorrectionShell;

	//SAMPLING FOR OUTPUT:
	std::vector<double> 	   _VirN_accum;
	std::vector<double> 	   _VirT_accum;
	std::vector<double> 	   _VirX_accum;
	std::vector<double> 	   _VirY_accum;
	std::vector<double> 	   _VirZ_accum;
	std::vector<double> 	   _velocityN_accum;
	std::vector<unsigned long> _numMolecules_accum_output;   



    void initVectors() // Resize vectors & set to 0
	{  
		//TODO: alle vektoren hier vorhanden?
		_shellVolume.resize(_lenVector);
		_shellLowerBound.resize(_lenVector);
		_numMolecules_accum.resize(_lenVector);
		_density_avg.resize(_lenVector);
		_density_avg_fitted.resize(_lenVector);
		_virtual_density.resize(_lenVector);
		_FcorrectionShell.resize(_lenVector);
		_UcorrectionShell.resize(_lenVector);
		_VirNcorrectionShell.resize(_lenVector);
		_VirTcorrectionShell.resize(_lenVector);
		std::fill(_shellVolume.begin(), _shellVolume.end(), 0.); 				 
		std::fill(_shellLowerBound.begin(), _shellLowerBound.end(), 0.);		 
		std::fill(_numMolecules_accum.begin(), _numMolecules_accum.end(), 0ul);  
		std::fill(_density_avg.begin(), _density_avg.end(), 0.);				 
		std::fill(_density_avg_fitted.begin(), _density_avg_fitted.end(), 0.); 	
		std::fill(_virtual_density.begin(), _virtual_density.end(), 0.);		 
		std::fill(_FcorrectionShell.begin(), _FcorrectionShell.end(), 0.);		 
		std::fill(_UcorrectionShell.begin(), _UcorrectionShell.end(), 0.);		 
		std::fill(_VirNcorrectionShell.begin(), _VirNcorrectionShell.end(), 0.);		 
		std::fill(_VirTcorrectionShell.begin(), _VirTcorrectionShell.end(), 0.);		 
		
		//sampling for output:		
		_VirN_accum.resize(_lenVector);
		_VirT_accum.resize(_lenVector);
		_VirX_accum.resize(_lenVector);
		_VirY_accum.resize(_lenVector);
		_VirZ_accum.resize(_lenVector);
		_velocityN_accum.resize(_lenVector);
		_numMolecules_accum_output.resize(_lenVector);
		std::fill(_VirN_accum.begin(), _VirN_accum.end(), 0.);		 
		std::fill(_VirT_accum.begin(), _VirT_accum.end(), 0.);		 
		std::fill(_VirX_accum.begin(), _VirX_accum.end(), 0.);		 
		std::fill(_VirY_accum.begin(), _VirY_accum.end(), 0.);		 
		std::fill(_VirZ_accum.begin(), _VirZ_accum.end(), 0.);		 
		std::fill(_velocityN_accum.begin(), _velocityN_accum.end(), 0.);		 
		std::fill(_numMolecules_accum_output.begin(), _numMolecules_accum_output.end(), 0.);		 
	}; 

	void resetDensityVectors()
	{
		std::fill(_numMolecules_accum.begin(), _numMolecules_accum.end(), 0ul); // korrekte 0 einfüllen: 0.0f, 0ul, ...
		std::fill(_density_avg.begin(), _density_avg.end(), 0.); // korrekte 0 einfüllen: 0.0f, 0ul, ...
		std::fill(_density_avg_fitted.begin(), _density_avg_fitted.end(), 0.); // korrekte 0 einfüllen: 0.0f, 0ul, ...
		std::fill(_virtual_density.begin(), _virtual_density.end(), 0.); // korrekte 0 einfüllen: 0.0f, 0ul, ...
	};


    unsigned long _numsamples {0};  // Number of samples for correction (density)


	struct BulkBoundaries
	{
		int inside_from;
		int inside_to;
		int outside_from;
		int outside_to;
	} _bulkBoundaries;

	double _rho_in {0.};
	double _rho_out {0.};


	//for testing:
	// void calculation_01_directCalc();
	void calculation_2_getShellPropertiesByLoopingOverAllMolecules();
	void calculation_3_getShellPropertiesFromRepresentativeParticles();
	void calculation_4_isabelReverseEngineered();

	double TICCu(int n, double rc, double sigma2);
	double TICSu(int n, double rc, double sigma2, double tau);
	double TISSu(int n, double rc, double sigma2, double tau1, double tau2);
	double TICCp(int n, double rc, double sigma2);
	double TICSp(int n, double rc, double sigma2, double tau);
	double TISSp(int n, double rc, double sigma2, double tau1, double tau2);

};