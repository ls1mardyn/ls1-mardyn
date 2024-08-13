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
	void resizeCommVarVector(CommVar<std::vector<T>>& commVar, unsigned int len) const {
		commVar.local.resize(len);
		commVar.global.resize(len);
	}
	template<typename T>
	void resizeCommVarVector(CommVar<std::vector<T>>& commVar) const {
		commVar.local.resize(_lenVector);
		commVar.global.resize(_lenVector);
	}


    void resizeVectors();  // Change size of accumulation vectors
    void resetVectors();   // Set accumulation vectors to zero

	void calculateBulkDensities();
	void calculateTanhProfile();



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

    std::vector<double> _shellVolume;             
    std::vector<double> _shellLowerBound; //_shellUpperBound[i] = _shellLowerBound[i+1]            

	std::vector<unsigned long> _numMolecules_accum;             // Number of molecules in bin
	std::vector<double> _density_avg; 	//avg over time (TODO: so far: all time average -- should be changed in the future)
	std::vector<double> _density_avg_fitted; 
	std::vector<double> _virtual_density;   // _density_avg_fitted - homogenous density


    unsigned long _numsamples {0};                   // Number of samples; can vary from bin to bin as some bins could be empty
	
	struct BulkBoundaries
	{
		int inside_from;
		int inside_to;
		int outside_from;
		int outside_to;
	} _bulkBoundaries;

	double _rho_in {0.};
	double _rho_out {0.};
};