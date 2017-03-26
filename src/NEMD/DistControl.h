/*
 * DistControl.h
 *
 *  Created on: 16.03.2015
 *      Author: mheinen
 */

#ifndef DISTCONTROL_H_
#define DISTCONTROL_H_

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "utils/ObserverBase.h"
#include "utils/Region.h"
#include "molecules/MoleculeForwardDeclaration.h"

using namespace std;

class Domain;
class ParticleContainer;
class DomainDecompBase;

enum DistControlUpdateMethods
{
	DCUM_UNKNOWN = 0,
	DCUM_DENSITY_PROFILE = 1,
	DCUM_DENSITY_PROFILE_DERIVATION = 2,
//	DCUM_FORCE_PROFILE = 3,
};

enum DistControlInitMethods
{
	DCIM_UNKNOWN = 0,
	DCIM_START_CONFIGURATION = 1,
	DCIM_MIDPOINT_VALUES = 2,
	DCIM_READ_FROM_FILE = 3,
};

class DistControl : public ControlInstance, public SubjectBase
{
public:
    DistControl(DomainDecompBase* domainDecomp, Domain* domain, unsigned int nUpdateFreq, unsigned int nWriteFreqProfiles);
    ~DistControl();

    std::string GetShortName() {return "DiC";}

    // set subdivision
	void SetSubdivision(unsigned int nNumSlabs) {_nNumShells = nNumSlabs; _nSubdivisionOpt = SDOPT_BY_NUM_SLABS;}
    void SetSubdivision(double dSlabWidth) {_dShellWidth = dSlabWidth; _nSubdivisionOpt = SDOPT_BY_SLAB_WIDTH;}
    void PrepareSubdivision();  // need to be called before PrepareDataStructures()

    // data structures
    void PrepareDataStructures();

    // init
    void InitPositions(double dInterfaceMidLeft, double dInterfaceMidRight);

    double GetInterfaceMidLeft() {return _dInterfaceMidLeft;}
    double GetInterfaceMidRight() {return _dInterfaceMidRight;}
    unsigned int GetUpdateFreq() {return _nUpdateFreq;}
    unsigned int GetWriteFreqProfiles() {return _nWriteFreqProfiles;}

    // set update/init method
    void SetUpdateMethod(const std::string& strMethod, const std::stringstream& sstr);
    void SetInitMethod(const int& nMethod, const std::stringstream& sstr) {_nMethodInit = nMethod; _sstrInit << sstr.rdbuf();}
    int GetInitMethod() {return _nMethodInit;}

    void Init(ParticleContainer* particleContainer);
    void WriteHeader();
    void WriteData(unsigned long simstep);
    void WriteDataProfiles(unsigned long simstep);

    // place method inside loop over molecule container
    void SampleProfiles(Molecule* mol);

    void UpdatePositionsInit(ParticleContainer* particleContainer);  // call in Simulation::prepare_start()
    void UpdatePositions(unsigned long simstep);
    void AlignSystemCenterOfMass(Molecule* mol, unsigned long simstep);

    // SubjectBase methods
	virtual void registerObserver(ObserverBase* observer);
	virtual void deregisterObserver(ObserverBase* observer);
	virtual void informObserver();

private:
    // place methods after the loop
	void CalcProfiles();
    void EstimateInterfaceMidpoint();  // called by UpdatePositions
    void EstimateInterfaceMidpointsByForce();
    void ResetLocalValues();

    // data structures
    void InitDataStructurePointers();
    void AllocateDataStructures();
	void InitDataStructures();

	// processing profiles
	void SmoothProfiles(const unsigned int& nNeighbourVals);
	void DerivateProfiles(const unsigned int& nNeighbourVals);

private:
    double _dInterfaceMidLeft;
    double _dInterfaceMidRight;

    unsigned long* _nNumMoleculesLocal;
    unsigned long* _nNumMoleculesGlobal;
    double* _dMidpointPositions;
	double* _dForceSumLocal;
	double* _dForceSumGlobal;
    double* _dDensityProfile;
    double* _dDensityProfileSmoothed;
    double* _dDensityProfileSmoothedDerivation;
	double* _dForceProfile;
    double* _dForceProfileSmoothed;
    unsigned int _nNumShells;
    double _dShellWidth;
    double _dInvertShellWidth;
    double _dShellVolume;
    unsigned int _nUpdateFreq;
    unsigned int _nWriteFreq;
    unsigned int _nWriteFreqProfiles;

	// update method
    int _nMethod;
	double _dVaporDensity;
	unsigned short _nNeighbourValsSmooth;
	unsigned short _nNeighbourValsDerivate;

    int _nMethodInit;
    std::stringstream _sstrInit;
	std::string _strFilenameInit;
	unsigned long _nRestartTimestep;

    // write data
    std::string _strFilename;
    std::string _strFilenameProfilesPrefix;

    // simtype
    int _nSimType;

    // observer
	std::vector<ObserverBase*> _observer;

	int _nSubdivisionOpt;

};  // class DistControl


#endif /* DISTCONTROL_H_ */
