/*
 * DistControl.cpp
 *
 *  Created on: 16.03.2015
 *      Author: mheinen
 */

#include "DistControl.h"
#include "molecules/Molecule.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "utils/Region.h"
#include "utils/DynAlloc.h"
#include "utils/Math.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric>
#include <cstdlib>
#include <limits>

using namespace std;


DistControl::DistControl(DomainDecompBase* domainDecomp, Domain* domain, unsigned int nUpdateFreq, unsigned int nWriteFreqProfiles)
: ControlInstance(domain, domainDecomp)
{
	// update frequency
    _nUpdateFreq = nUpdateFreq;

    // write data
    _nWriteFreq = _nUpdateFreq;
    _nWriteFreqProfiles = nWriteFreqProfiles;

    _strFilename = "DistControl.dat";
    _strFilenameProfilesPrefix = "DistControlProfiles";

    // init method variables
    _sstrInit.clear();
	_strFilenameInit = "unknown";
	_nRestartTimestep = 0;

	// init data structure pointers
	this->InitDataStructurePointers();
}

DistControl::~DistControl()
{

}

void DistControl::PrepareSubdivision()
{
	Domain* domain = this->GetDomain();
	double dWidth = domain->getGlobalLength(1);

	switch(_nSubdivisionOpt)
	{
	case SDOPT_BY_NUM_SLABS:
		_dShellWidth = dWidth / ( (double)(_nNumShells) );
		break;
	case SDOPT_BY_SLAB_WIDTH:
		_nNumShells = round(dWidth / _dShellWidth);
		_dShellWidth = dWidth / ( (double)(_nNumShells) );
		break;
	case SDOPT_UNKNOWN:
	default:
		global_log->error() << "ERROR in DistControl::PrepareSubdivision(): Neither _dShellWidth nor _nNumShells was set correctly! Programm exit..." << endl;
		exit(-1);
	}

	_dInvertShellWidth = 1. / _dShellWidth;
    _dShellVolume = _dShellWidth * domain->getGlobalLength(0) * domain->getGlobalLength(2);

//	cout << "DistControl::_nNumShells = " << _nNumShells << endl;
}

void DistControl::InitDataStructurePointers()
{
    // number of molecules
	_nNumMoleculesLocal = NULL;
	_nNumMoleculesGlobal = NULL;

	// density profile
	_dDensityProfile = NULL;
	_dDensityProfileSmoothed = NULL;
	_dDensityProfileSmoothedDerivation = NULL;

    // force profile
	_dForceSumLocal = NULL;
	_dForceSumGlobal = NULL;
	_dForceProfile = NULL;
	_dForceProfileSmoothed = NULL;

	// profile midpoint positions
	_dMidpointPositions = NULL;
}

void DistControl::AllocateDataStructures()
{
    // number of molecules
	AllocateUnsLongArray(_nNumMoleculesLocal, _nNumShells);
	AllocateUnsLongArray(_nNumMoleculesGlobal, _nNumShells);

	// density profile
	AllocateDoubleArray(_dDensityProfile, _nNumShells);
	AllocateDoubleArray(_dDensityProfileSmoothed, _nNumShells);
	AllocateDoubleArray(_dDensityProfileSmoothedDerivation, _nNumShells);

    // force profile
	AllocateDoubleArray(_dForceSumLocal, _nNumShells);
	AllocateDoubleArray(_dForceSumGlobal, _nNumShells);
	AllocateDoubleArray(_dForceProfile, _nNumShells);
	AllocateDoubleArray(_dForceProfileSmoothed, _nNumShells);

	// profile midpoint positions
	AllocateDoubleArray(_dMidpointPositions, _nNumShells);
}

void DistControl::InitDataStructures()
{
    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
		// number of molecules
		_nNumMoleculesLocal[s] = 0;
		_nNumMoleculesGlobal[s] = 0;

		// density profile
		_dDensityProfile[s] = 0.;
		_dDensityProfileSmoothed[s] = 0.;
		_dDensityProfileSmoothedDerivation[s] = 0.;

		// force profile
		_dForceSumLocal[s] = 0.;
		_dForceSumGlobal[s] = 0.;
		_dForceProfile[s] = 0.;
		_dForceProfileSmoothed[s] = 0.;

		// profile midpoint positions
		_dMidpointPositions[s] = (0.5 + s)*_dShellWidth;
    }
}

void DistControl::PrepareDataStructures()
{
	this->AllocateDataStructures();
	this->InitDataStructures();
}

// update method
void DistControl::SetUpdateMethod(const std::string& strMethod, const std::stringstream& sstr)
{
	std::stringstream sstrcopy;
	sstrcopy << sstr.rdbuf();
	bool bWrongNumArg = false;
	std::string strExpected = "";

	if(strMethod == "density")
	{
		_nMethod = DCUM_DENSITY_PROFILE;
		sstrcopy >> _dVaporDensity;
		bWrongNumArg = (sstrcopy.fail() || !(sstrcopy.eof() ) );
		strExpected = "Expected: 'DistControl method density <double>'. ";
	}
	else if(strMethod == "denderiv")
	{
		_nMethod = DCUM_DENSITY_PROFILE_DERIVATION;
		sstrcopy >> _nNeighbourValsSmooth;
		sstrcopy >> _nNeighbourValsDerivate;
		bWrongNumArg = (sstrcopy.fail() || !(sstrcopy.eof() ) );
		strExpected = "Expected: 'DistControl method denderiv <unsigned short> <unsigned short>'. ";
	}
	else
	{
		strExpected = "Expected: 'DistControl method density|denderiv ... '. ";
		cout << "ERROR in DistControl. " << strExpected << "Programm exit..." << endl;
		exit(-1);
	}

	// check number of arguments
	if( true == bWrongNumArg)
	{
		cout << "ERROR in DistControl. " << strExpected << "Programm exit..." << endl;
		exit(-1);
	}
}

void DistControl::SampleProfiles(Molecule* mol)
{
    unsigned int nPosIndex;
    unsigned int nIndexMax = _nNumShells - 1;

    // calc position index
    double dPos = mol->r(1);

    nPosIndex = (unsigned int) floor(dPos * _dInvertShellWidth);

    // ignore outer (halo) molecules
    if(nPosIndex > nIndexMax)  // negative values will be ignored to: cast to unsigned int --> high value
        return;

    _nNumMoleculesLocal[nPosIndex]++;
    _dForceSumLocal[nPosIndex] += mol->F(1);
}

void DistControl::CalcProfiles()
{
	#ifdef ENABLE_MPI

    MPI_Allreduce( _nNumMoleculesLocal, _nNumMoleculesGlobal, _nNumShells, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce( _dForceSumLocal, _dForceSumGlobal, _nNumShells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
        _nNumMoleculesGlobal[s] = _nNumMoleculesLocal[s];
        _dForceSumGlobal[s] = _dForceSumLocal[s];
    }
#endif

    // calc profiles
    double dInvertShellVolume = 1. / _dShellVolume;
    double dInvertSampleTimesteps = 1. / ( (double)(_nUpdateFreq) );

    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
        _dDensityProfile[s] = _nNumMoleculesGlobal[s] * dInvertSampleTimesteps * dInvertShellVolume;

        if(_nNumMoleculesGlobal[s] > 0) 
            _dForceProfile[s] = _dForceSumGlobal[s] / (double)(_nNumMoleculesGlobal[s]) * dInvertSampleTimesteps;
        else
            _dForceProfile[s] = 0.;
    }

	// smooth profiles
	this->SmoothProfiles(_nNeighbourValsSmooth);

    // density derivation
    this->DerivateProfiles(_nNeighbourValsDerivate);
}

void DistControl::EstimateInterfaceMidpointsByForce()
{
    // find min/max in drho/dy profile
    double dMin = 0;
    double dMax = 0;
    unsigned int nIndexMin = 0;
    unsigned int nIndexMax = 0;

    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
//        if(_dDensityProfileSmoothed[s] > _dVaporDensity)
//        {
		double drhody = _dDensityProfileSmoothedDerivation[s];

		if(drhody > dMax)
		{
			dMax = drhody;
			nIndexMax = s;
		}

		if(drhody < dMin)
		{
			dMin = drhody;
			nIndexMin = s;
		}
//        }
    }

/*
    // find min/max in force profile
    double dFmin = 0;
    double dFmax = 0;
    unsigned int nIndexMin = 0;
    unsigned int nIndexMax = 0;

    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
        if(_dDensityProfileSmoothed[s] > _dVaporDensity)
        {
            double dF = _dForceProfileSmoothed[s];

            if(dF > dFmax)
            {
                dFmax = dF;
                nIndexMax = s;
            }

            if(dF < dFmin)
            {
                dFmin = dF;
                nIndexMin = s;
            }
        }
    }
*/

    _dInterfaceMidLeft  = _dMidpointPositions[nIndexMax];
    _dInterfaceMidRight = _dMidpointPositions[nIndexMin];
}

void DistControl::EstimateInterfaceMidpoint()
{

/*
    // mheinen_2015-03-17 --> DEBUG_TEST_CASE
    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
        _dDensityProfile[s] = 0.7;
    }

    _dDensityProfile[0] = 0.1;
    _dDensityProfile[1] = 0.1;
    _dDensityProfile[2] = 0.1;
    _dDensityProfile[3] = 0.1;
    _dDensityProfile[4] = 0.2;
    _dDensityProfile[5] = 0.3;
    _dDensityProfile[6] = 0.4;
    _dDensityProfile[7] = 0.5;
    _dDensityProfile[8] = 0.6;

    _dDensityProfile[211] = 0.6;
    _dDensityProfile[212] = 0.5;
    _dDensityProfile[213] = 0.4;
    _dDensityProfile[214] = 0.3;
    _dDensityProfile[215] = 0.2;
    _dDensityProfile[216] = 0.1;
    _dDensityProfile[217] = 0.1;
    _dDensityProfile[218] = 0.1;
    _dDensityProfile[219] = 0.1;


    // result:
    // ----------
    // nIndexSlabGreater = 7
    // x1 = 6.5
    // dx1m = 0
    // _dInterfaceMidLeft = 6.5
    // _dInterfaceMidRight = 213.5

    // <-- DEBUG_TEST_CASE
*/


    // determine left midpoint
    {
        // determine min max
         double dMin = _dDensityProfile[0];
         double dMax = _dDensityProfile[0];

//         int nIndexMin = 0;
//         int nIndexMax = 0;

         for(unsigned int s = 0; s < _nNumShells/2; ++s)
         {
             double dDensity = _dDensityProfile[s];

             // find min
             if(dDensity < dMin)
             {
                 dMin = dDensity;
//                 nIndexMin = s;
             }

             // find max
             if(dDensity > dMax)
             {
                 dMax = dDensity;
//                 nIndexMax = s;
             }
         }

         // dN/dShellWidth = (Nsoll - Nlower) / dInterpol
         //
         // => dInterpol = (Nsoll - Nlower) / dN * dShellWidth
         //
         // dMidPos = dLeftBound + nNumSlabs_lower * dShellWidth + dInterpol

         // 1. ym bestimmen
         double ym = dMin + (dMax - dMin) / 2;

//         cout << "dMin = " << dMin << endl;
//         cout << "dMax = " << dMax << endl;
//         cout << "ym = " << ym << endl;

         // 2. Slab mit gerade niedrigerer Anzahl finden --> von links angefangen ersten Slab größer numMax finden --> Index -1 rechnen
         int nIndexSlabGreater = 0;

         for(unsigned int s = 0; s < _nNumShells; ++s)
         {
             if(_dDensityProfile[s] > ym)
             {
                 nIndexSlabGreater = s;
                 break;
             }
         }

         if(nIndexSlabGreater < 1)
         {
             cout << "ERROR in MEXRegion::EstimateInterfaceMidpoint(): could not find valid index in density profile" << endl;
             return;
         }

//         cout << "nIndexSlabGreater = " << nIndexSlabGreater << endl;

         //3. Interpolierte Distanz dInterpol berechnen, die von dem Slab mit N < Nsoll ausgehend in Richtung N > Nsoll gegangen werden muss
         double dx12, dx1m, y1, y2, dy12, dy1m;

         y1 = _dDensityProfile[nIndexSlabGreater-1];
         y2 = _dDensityProfile[nIndexSlabGreater];

         dy1m = ym - y1;

         dy12 = y2 - y1;
         dx12 = _dShellWidth;

         // dy12 / dx12 = dy1m / dx1m => dx1m = dy1m / dy12 * dx12
         dx1m = dy1m / dy12 * dx12;

         // 4. Berechnung der Mittelpunktsposition
         // nNumSlabslower == nIndexSlabLower + 1 == nIndexSlabGreater
         double xm, x1;
         double dOuterBoundarySampleZone = 0.;  //this->GetLowerCorner()[1];

         x1 = dOuterBoundarySampleZone + nIndexSlabGreater * _dShellWidth - _dShellWidth*0.5;

//         cout << "x1 = " << x1 << endl;
//         cout << "dx1m = " << dx1m << endl;

         xm = x1 + dx1m;
         _dInterfaceMidLeft = xm;

//         cout << "_dInterfaceMidLeft = " << _dInterfaceMidLeft << endl;
     }


    // determine right midpoint
    {
        // determine min max
        double dMin = _dDensityProfile[0];
        double dMax = _dDensityProfile[0];

//         int nIndexMin = 0;
//         int nIndexMax = 0;

         for(unsigned int s = _nNumShells/2; s < _nNumShells; ++s)
         {
             double dDensity = _dDensityProfile[s];

             // find min
             if(dDensity < dMin)
             {
                 dMin = dDensity;
//                 nIndexMin = s;
             }

             // find max
             if(dDensity > dMax)
             {
                 dMax = dDensity;
//                 nIndexMax = s;
             }
         }

         // dN/dShellWidth = (Nsoll - Nlower) / dInterpol
         //
         // => dInterpol = (Nsoll - Nlower) / dN * dShellWidth
         //
         // dMidPos = dLeftBound + nNumSlabs_lower * dShellWidth + dInterpol

         // 1. ym bestimmen
         double ym = dMin + (dMax - dMin) / 2;

//         cout << "ym = " << ym << endl;

         // 2. Slab mit gerade niedrigerer Anzahl finden --> von links angefangen ersten Slab größer numMax finden --> Index -1 rechnen
         unsigned int nIndexSlabGreater = 0;

         for(int s = _nNumShells-1; s >= 0; --s)
         {
             if(_dDensityProfile[s] > ym)
             {
                 nIndexSlabGreater = s;
                 break;
             }
         }

         if(nIndexSlabGreater < 1)
         {
             cout << "ERROR in MEXRegion::EstimateInterfaceMidpoint(): could not find valid index in density profile" << endl;
             return;
         }

         //3. Interpolierte Distanz dInterpol berechnen, die von dem Slab mit N < Nsoll ausgehend in Richtung N > Nsoll gegangen werden muss
         double dx12, dx1m, y1, y2, dy12, dy1m;

         y1 = _dDensityProfile[nIndexSlabGreater+1];
         y2 = _dDensityProfile[nIndexSlabGreater];

         dy1m = ym - y1;

         dy12 = y2 - y1;
         dx12 = _dShellWidth;

         // dy12 / dx12 = dy1m / dx1m => dx1m = dy1m / dy12 * dx12
         dx1m = dy1m / dy12 * dx12;

         // 4. Berechnung der Mittelpunktsposition
         // nNumSlabslower == nIndexSlabLower + 1 == nIndexSlabGreater
         Domain* domain = this->GetDomain();
         double xm, x1;
         double dOuterBoundarySampleZone = domain->getGlobalLength(1);  //this->GetUpperCorner()[1];
         unsigned int numShellsLower = _nNumShells-1 - nIndexSlabGreater;

//         cout << "numShellsLower = " << numShellsLower << endl;

         x1 = dOuterBoundarySampleZone - numShellsLower * _dShellWidth + _dShellWidth*0.5;

//         cout << "x1 = " << x1 << endl;
//         cout << "dx1m = " << dx1m << endl;

         xm = x1 - dx1m;
         _dInterfaceMidRight = xm;

//         cout << "_dInterfaceMidRight = " << _dInterfaceMidRight << endl;
    }
}


// init
void DistControl::InitPositions(double dInterfaceMidLeft, double dInterfaceMidRight)
{
    // set interface midpoints manually
    _dInterfaceMidLeft  = dInterfaceMidLeft;
    _dInterfaceMidRight = dInterfaceMidRight;
}

void DistControl::UpdatePositionsInit(ParticleContainer* particleContainer)
{
	switch(_nMethodInit )
	{
	case DCIM_START_CONFIGURATION:

		for( ParticleIterator tM  = particleContainer->iteratorBegin();
			 tM != particleContainer->iteratorEnd();
			 ++tM)
		{
			// sample density profile
			this->SampleProfiles(&(*tM));
		}

		// determine interface midpoints and update region positions
		this->UpdatePositions(this->GetUpdateFreq() );
		break;
	case DCIM_MIDPOINT_VALUES:

		_sstrInit >> _dInterfaceMidLeft;
		_sstrInit >> _dInterfaceMidRight;
		break;
	case DCIM_READ_FROM_FILE:
		{
		_sstrInit >> _strFilenameInit;
		_sstrInit >> _nRestartTimestep;
//		cout << "_strFilenameInit = " << _strFilenameInit << endl;
//		cout << "_nRestartTimestep = " << _nRestartTimestep << endl;

		ifstream filein(_strFilenameInit.c_str(), ios::in);

		string strLine, strToken;
		string strTokens[20];

		while (getline (filein, strLine))
		{
			stringstream sstr;
			sstr << strLine;

			sstr >> strToken;

			if( (unsigned long) (atoi(strToken.c_str() ) ) == _nRestartTimestep)
			{
				int i=1;
				while (sstr >> strTokens[i])
					i++;
			}
		}

		_dInterfaceMidLeft  = atof(strTokens[1].c_str() );
		_dInterfaceMidRight = atof(strTokens[2].c_str() );

		filein.close();
		break;
		}
	case DCIM_UNKNOWN:
	default:
		cout << "DistControl: Wrong Init Method! Programm exit..." << endl;
		exit(-1);
	}

#ifndef NDEBUG
	cout << "DistControl::_dInterfaceMidLeft = " << _dInterfaceMidLeft << endl;
	cout << "DistControl::_dInterfaceMidRight = " << _dInterfaceMidRight << endl;
#endif

    // update positions
    this->informObserver();

    // reset local values
    this->ResetLocalValues();
}

void DistControl::UpdatePositions(unsigned long simstep)
{
    // update with respect to update frequency
    if(simstep % _nUpdateFreq != 0 || simstep == 0 || simstep == 1)  // TODO init timestep
        return;

    // calc profiles
    this->CalcProfiles();

    // update midpoint coordinates with respect to desired method
    switch(_nMethod)
    {
    case DCUM_DENSITY_PROFILE:
		this->EstimateInterfaceMidpoint();
		break;
    case DCUM_DENSITY_PROFILE_DERIVATION:
    	this->EstimateInterfaceMidpointsByForce();
    	break;
    case DCUM_UNKNOWN:
    default:
		cout << "DistControl::UpdatePositions: Corrupted code!!! Programm exit..." << endl;
    	exit(-1);
    }

    // update positions
    this->informObserver();

    // reset local values
    this->ResetLocalValues();
}


void DistControl::ResetLocalValues()
{
    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
        _nNumMoleculesLocal[s] = 0;
        _dForceSumLocal[s] = 0.;
    }
}

void DistControl::WriteData(unsigned long simstep)
{
	// domain decomposition
	DomainDecompBase* domainDecomp = this->GetDomainDecomposition();

    // write out data
    stringstream outputstream;

    // write data
    if(simstep % _nWriteFreq == 0)
    {
        // calc global values
//        this->CalcGlobalValues(domainDecomp);

#ifdef ENABLE_MPI
    int rank = domainDecomp->getRank();
    // int numprocs = domainDecomp->getNumProcs();
    if (rank == 0)
    {
#endif
        outputstream << std::setw(20) << simstep;
        outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dInterfaceMidLeft;
        outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dInterfaceMidRight;

        // observer data
    	std::vector<ObserverBase*>::iterator it;

    	for(it=_observer.begin(); it!=_observer.end(); it++)
    	{
    		CuboidRegionObs* region = static_cast<CuboidRegionObs*>(*it);
    		if(NULL != region)
    		{
    	        outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << region->GetLowerCorner(1);
    	        outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << region->GetUpperCorner(1);
    		}
    	}

        outputstream << endl;

        ofstream fileout(_strFilename.c_str(), ios::out|ios::app);
        fileout << outputstream.str();
        fileout.close();
#ifdef ENABLE_MPI
    }
#endif
    }
}

void DistControl::WriteHeader()
{
	// domain decomposition
	DomainDecompBase* domainDecomp = this->GetDomainDecomposition();

    // write header
    stringstream outputstream;

#ifdef ENABLE_MPI
    int rank = domainDecomp->getRank();
    // int numprocs = domainDecomp->getNumProcs();
    if (rank == 0)
    {
#endif

    outputstream << "             simstep";
    outputstream << "                 midLeft";
    outputstream << "                midRight";

    // observer data
	std::vector<ObserverBase*>::iterator it;

	for(it=_observer.begin(); it!=_observer.end(); it++)
	{
		CuboidRegionObs* region = static_cast<CuboidRegionObs*>(*it);
		if(NULL != region)
		{
			outputstream << "                 " << region->GetParent()->GetShortName() << "l" << "[" << region->GetID() << "]";
			outputstream << "                 " << region->GetParent()->GetShortName() << "r" << "[" << region->GetID() << "]";
		}
	}

    outputstream << endl;

    ofstream fileout(_strFilename.c_str(), ios::out);
    fileout << outputstream.str();
    fileout.close();

#ifdef ENABLE_MPI
    }
#endif
}


void DistControl::WriteDataProfiles(unsigned long simstep)
{
	// domain decomposition
	DomainDecompBase* domainDecomp = this->GetDomainDecomposition();

    // write out data
    stringstream outputstream;
    stringstream filenamestream;

    filenamestream << _strFilenameProfilesPrefix << "_TS";
    filenamestream.fill('0');
    filenamestream.width(9);
    filenamestream << right << simstep;
    filenamestream << ".dat";

    // write data
    if(simstep % _nWriteFreqProfiles == 0)
    {
        // calc global values
//        this->CalcGlobalValues(domainDecomp);

#ifdef ENABLE_MPI
    int rank = domainDecomp->getRank();
    // int numprocs = domainDecomp->getNumProcs();
    if (rank == 0)
    {
#endif

        // write header
        outputstream << "                   coord";
        outputstream << "                     rho";
        outputstream << "              rho_smooth";
        outputstream << "                 drho/dy";
        outputstream << "                      Fy";
        outputstream << "               Fy_smooth";
        outputstream << endl;

        for(unsigned int s=0; s<_nNumShells; s++)
        {
            outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dMidpointPositions[s];
            outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dDensityProfile[s];
            outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dDensityProfileSmoothed[s];
            outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dDensityProfileSmoothedDerivation[s];
            outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dForceProfile[s];
            outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dForceProfileSmoothed[s];
            outputstream << endl;
        }

        ofstream fileout(filenamestream.str().c_str(), ios::out|ios::out);
        fileout << outputstream.str();
        fileout.close();
#ifdef ENABLE_MPI
    }
#endif
    }
}


void DistControl::Init(ParticleContainer* particleContainer)
{
    // write output file with header
    this->WriteHeader();
}


void DistControl::AlignSystemCenterOfMass(Molecule* mol, unsigned long simstep)
{
	// domain
	Domain* domain = this->GetDomain();

    // update with respect to update frequency
    if(simstep % _nUpdateFreq != 0 || simstep == 0 || simstep == 1)  // TODO init timestep
        return;

    double dDeltaLeftY  = _dInterfaceMidLeft;  // (minus 0)
    double dDeltaRightY = domain->getGlobalLength(1) - _dInterfaceMidRight;
    double dDeltaY = 0.5 * (dDeltaRightY - dDeltaLeftY);

    double dNewPosition = mol->r(1) + dDeltaY;

/*
 * seems not to work: simulation crashes!!!
 *
    double dBoxLengthY = domain->getGlobalLength(1);

    if (dNewPosition > dBoxLengthY )
        dNewPosition -= dBoxLengthY;
    else if (dNewPosition < 0. )
        dNewPosition += dBoxLengthY;
*/

    mol->setr(1, dNewPosition);

#ifdef DEBUG
    cout << "_dInterfaceMidLeft = " << _dInterfaceMidLeft << endl;
    cout << "_dInterfaceMidRight = " << _dInterfaceMidRight << endl;
    cout << "dDeltaLeftY = " << dDeltaLeftY << endl;
    cout << "dDeltaRightY = " << dDeltaRightY << endl;
    cout << "dDeltaY = " << dDeltaY << endl;
#endif
}

void DistControl::registerObserver(ObserverBase* observer)
{
	_observer.push_back(observer);
}

void DistControl::deregisterObserver(ObserverBase* observer)
{
	std::vector<ObserverBase*>::iterator it;

	for(it=_observer.begin(); it!=_observer.end(); it++)
	{
		if( (*it) == observer)
		{
			_observer.erase(it);
		}
	}
}

void DistControl::informObserver()
{
	std::vector<ObserverBase*>::iterator it;

	for(it=_observer.begin(); it!=_observer.end(); it++)
	{
		(*it)->set(_dInterfaceMidLeft, _dInterfaceMidRight);
	}
}

void DistControl::SmoothProfiles(const unsigned int& nNeighbourVals)
{
	// smooth profiles
	double dNumSmoothedVals = (double)(2.*nNeighbourVals+1);

	unsigned int l = 0;
	unsigned int r = nNeighbourVals;

	for(unsigned int s = 0; s < _nNumShells; ++s)
	{
		double dDensitySum = 0.;
		double dForceSum   = 0.;

		for(unsigned int t = s-l; t <= s+r; ++t)
		{
			dDensitySum += _dDensityProfile[t];
			dForceSum   += _dForceProfile[t];
		}
		_dDensityProfileSmoothed[s] = dDensitySum / dNumSmoothedVals;
		_dForceProfileSmoothed[s]   = dForceSum   / dNumSmoothedVals;

		if(l < nNeighbourVals)
			l++;

		if(s >= (_nNumShells-1 - nNeighbourVals) )
			r--;
	}
}

void DistControl::DerivateProfiles(const unsigned int& nNeighbourVals)
{
    unsigned int l = 0;
    unsigned int r = nNeighbourVals;

/*
 * derivate by two points
 *
    for(unsigned int s = 0; s < _nNumShells; ++s)
    {
        double dy = _dDensityProfileSmoothed[s+r] - _dDensityProfileSmoothed[s-l];
        double dx = (r+l)*_dShellWidth;

        _dDensityProfileSmoothedDerivation[s] = dy / dx;

        if(l < nNumNeighbourVals)
            l++;

        if(s >= (_nNumShells-1 - nNumNeighbourVals) )
            r--;
    }
*/

	vector<double> x;
	vector<double> y;

	for(unsigned int s = 0; s < _nNumShells; ++s)
    {
		x.clear();
		y.clear();

		for(unsigned int i = s-l; i<=s+r; ++i)
    	{
			x.push_back(_dMidpointPositions[i] );
			y.push_back(_dDensityProfileSmoothed[i] );
    	}
		double beta1, beta2;
        LinearRegression(x, y, beta1, beta2);

        _dDensityProfileSmoothedDerivation[s] = beta2;

        if(l < nNeighbourVals)
            l++;

        if(s >= (_nNumShells-1 - nNeighbourVals) )
            r--;
    }
}












