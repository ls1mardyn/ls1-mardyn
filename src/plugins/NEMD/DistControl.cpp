/*
 * DistControl.cpp
 *
 *  Created on: 16.03.2015
 *      Author: mheinen
 */

#include "DistControl.h"
#include "molecules/Molecule.h"
#include "Domain.h"
#include "Simulation.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "utils/Region.h"
#include "utils/DynAlloc.h"
#include "utils/Math.h"
#include "utils/xmlfileUnits.h"
#include "plugins/Mirror.h"  // TODO: Not so nice that DistControl has to know Mirror plugin explicitly

#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric>
#include <limits>
#include <cstdlib>
#include <cstdint>
#include <cmath>



DistControl::DistControl()
		: ControlInstance(),
		_dInterfaceMidLeft(0.),
		_dInterfaceMidRight(0.),
		_nNumComponents(0),
		_nTargetCompID(0),
		_nNumValuesScalar(0),
		_nMethod(DCUM_UNKNOWN),
		_dVaporDensity(0.),
		_nNeighbourValsSmooth(0),
		_nNeighbourValsDerivate(0),
		_nMethodInit(DCIM_UNKNOWN),
		_strFilenameInit("unknown"),
		_nRestartTimestep(0),
		_strFilename("unknown"),
		_strFilenameProfilesPrefix("unknown"),
		_nSubdivisionOpt(SDOPT_UNKNOWN)
{
	// number of components
	Domain* domain = global_simulation->getDomain();
	_nNumComponents = domain->getNumberOfComponents()+1;  // +1: because 0 stands for sum of all components;
}

DistControl::~DistControl()
{
}

void DistControl::readXML(XMLfileUnits& xmlconfig)
{
	// control
	_controlFreqs.update = 5000;
	_controlFreqs.write.profiles = 50000;
	xmlconfig.getNodeValue("control/update", _controlFreqs.update);
	xmlconfig.getNodeValue("control/writeprofiles", _controlFreqs.write.profiles);
	_controlFreqs.write.data = _controlFreqs.update;

	// filenames
	_strFilename = "DistControl.dat";
	_strFilenameProfilesPrefix = "DistControlProfiles";
	xmlconfig.getNodeValue("filenames/control",  _strFilename);
	xmlconfig.getNodeValue("filenames/profiles", _strFilenameProfilesPrefix);

	Log::global_log->info() << "[DistControl] Writing control data to file: " << _strFilename << std::endl;
	Log::global_log->info() << "[DistControl] Writing profile data to files with prefix: " << _strFilenameProfilesPrefix << std::endl;

	// subdivision of system
	uint32_t nSubdivisionType = SDOPT_UNKNOWN;
	std::string strSubdivisionType;
	if( !xmlconfig.getNodeValue("subdivision@type", strSubdivisionType) )
	{
		Log::global_log->error() << "[DistControl] Missing attribute \"subdivision@type\"! Programm exit..." << std::endl;
		exit(-1);
	}
	if("number" == strSubdivisionType)
	{
		unsigned int nNumSlabs = 0;
		if( !xmlconfig.getNodeValue("subdivision/number", nNumSlabs) )
		{
			Log::global_log->error() << "[DistControl] Missing element \"subdivision/number\"! Programm exit..." << std::endl;
			exit(-1);
		}
		else
			this->SetSubdivision(nNumSlabs);
	}
	else if("width" == strSubdivisionType)
	{
		double dSlabWidth = 0.;
		if( !xmlconfig.getNodeValue("subdivision/width", dSlabWidth) )
		{
			Log::global_log->error() << "[DistControl] Missing element \"subdivision/width\"! Programm exit..." << std::endl;
			exit(-1);
		}
		else
			this->SetSubdivision(dSlabWidth);
	}
	else
	{
		Log::global_log->error() << "[DistControl] Wrong attribute \"subdivision@type\". Expected: type=\"number|width\"! Programm exit..." << std::endl;
		exit(-1);
	}

	// init method
	_nMethodInit = DCIM_UNKNOWN;
	_dInterfaceMidLeft  = 10.0;
	_dInterfaceMidRight = 20.0;
	_strFilenameInit = "unknown";
	_nRestartTimestep = 1;
	std::string strInitMethodType;

	xmlconfig.getNodeValue("init@type", strInitMethodType);
	if("startconfig" == strInitMethodType)
	{
		_nMethodInit = DCIM_START_CONFIGURATION;
		Log::global_log->info() << "[DistControl] Init method 'startconfig', dertermining interface midpoints from start configuration." << std::endl;
	}
	else if("values" == strInitMethodType)
	{
		_nMethodInit = DCIM_MIDPOINT_VALUES;
		bool bInputIsValid = true;
		bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("init/values/left",  _dInterfaceMidLeft);
		bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("init/values/right", _dInterfaceMidRight);
		if(true == bInputIsValid)
		{
			Log::global_log->info() << "[DistControl] Init method 'values' => interface midpoint left: " << _dInterfaceMidLeft << ", "
					"right: " << _dInterfaceMidRight << "." << std::endl;
		}
		else
		{
			Log::global_log->error() << "[DistControl] Missing elements \"init/values/left\" or \"init/values/right\" or both! Programm exit..." << std::endl;
			exit(-1);
		}
	}
	else if("file" == strInitMethodType)
	{
		_nMethodInit = DCIM_READ_FROM_FILE;
		bool bInputIsValid = true;
		bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("init/file",    _strFilenameInit);
		bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("init/simstep", _nRestartTimestep);
		if(true == bInputIsValid)
		{
			Log::global_log->info() << "[DistControl] Init method 'file', reading from file: " << _strFilenameInit << ", "
					"goto line with simstep == " << _nRestartTimestep << "." << std::endl;
		}
		else
		{
			Log::global_log->error() << "[DistControl] Missing elements \"init/file\" or \"init/simstep\" or both! Programm exit..." << std::endl;
			exit(-1);
		}
	}
	else
	{
		Log::global_log->error() << "[DistControl] Wrong attribute \"init@type\", type = " << strInitMethodType << ", "
				"expected: type=\"startconfig|values|file\"! Programm exit..." << std::endl;
		exit(-1);
	}

	// update method
	_nMethod = DCUM_UNKNOWN;
	_nTargetCompID = 0;
	_dVaporDensity = 0.01;
	_nNeighbourValsSmooth   = 2;
	_nNeighbourValsDerivate = 2;
	std::string strUpdateMethodType;

	xmlconfig.getNodeValue("method@type", strUpdateMethodType);
	if("density" == strUpdateMethodType)
	{
		_nMethod = DCUM_DENSITY_PROFILE;
		bool bInputIsValid = true;
		uint32_t nTargetCompID = 0;
		bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("method/componentID", nTargetCompID);
		_nTargetCompID = (uint16_t)(nTargetCompID);
		bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("method/density",     _dVaporDensity);
		if(true == bInputIsValid)
		{
			Log::global_log->info() << "[DistControl] Update method 'density', using constant value for vapor density rho_vap == " << _dVaporDensity << ", "
					"target componentID: " << _nTargetCompID << "." << std::endl;
		}
		else
		{
			Log::global_log->error() << "[DistControl] Missing elements \"method/componentID\" or \"method/density\" or both! Programm exit..." << std::endl;
			exit(-1);
		}
	}
	else if("denderiv" == strUpdateMethodType)
	{
		_nMethod = DCUM_DENSITY_PROFILE_DERIVATION;
		bool bInputIsValid = true;
		uint32_t nTargetCompID = 0;
		uint32_t nNeighbourValsSmooth = 0;
		uint32_t nNeighbourValsDerivate = 0;
		bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("method/componentID", nTargetCompID);
		bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("method/neighbourvals[@algorithm='smooth']",   nNeighbourValsSmooth);
		bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("method/neighbourvals[@algorithm='derivate']", nNeighbourValsDerivate);
		_nTargetCompID = (uint16_t)(nTargetCompID);
		_nNeighbourValsSmooth = (uint16_t)(nNeighbourValsSmooth);
		_nNeighbourValsDerivate = (uint16_t)(nNeighbourValsDerivate);
		if(true == bInputIsValid)
		{
			Log::global_log->info() << "[DistControl] Update method 'denderiv', using " << _nNeighbourValsSmooth << " neigbour values for smoothing "
					" and " << _nNeighbourValsDerivate << " neigbour values for derivation of the density profile, target componentID: " << _nTargetCompID << "." << std::endl;
		}
		else
		{
			Log::global_log->error() << "[DistControl] Missing elements \"method/componentID\" or \"method/density\" or both! Programm exit..." << std::endl;
			exit(-1);
		}
	}
	else
	{
		Log::global_log->error() << "[DistControl] Wrong attribute \"method@type\", type = " << strUpdateMethodType << ", "
				"expected: type=\"density|denderiv\"! Programm exit..." << std::endl;
		exit(-1);
	}
}

void DistControl::init(ParticleContainer *particleContainer,
		  DomainDecompBase *domainDecomp, Domain *domain)
{
	this->Init(particleContainer);
	this->UpdatePositionsInit(particleContainer);
	this->WriteHeader();
	this->WriteData(0);
}

void DistControl::beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		unsigned long simstep)
{
	for(auto it = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it) {
		// sample density profile
		this->SampleProfiles(&(*it));
	}

	// determine interface midpoints and update region positions
	this->UpdatePositions(simstep);

	// write data
	this->WriteData(simstep);
	this->WriteDataProfiles(simstep);
}

void DistControl::PrepareSubdivision()
{
	Domain* domain = global_simulation->getDomain();
	double dWidth = domain->getGlobalLength(1);

	switch(_nSubdivisionOpt)
	{
	case SDOPT_BY_NUM_SLABS:
		_binParams.width = dWidth / ( (double)(_binParams.count) );
		break;
	case SDOPT_BY_SLAB_WIDTH:
		_binParams.count = round(dWidth / _binParams.width);
		_binParams.width = dWidth / ( (double)(_binParams.count) );
		break;
	case SDOPT_UNKNOWN:
	default:
		Log::global_log->error() << "[DistControl] PrepareSubdivision(): Neither _binParams.width nor _binParams.count was set correctly! Programm exit..." << std::endl;
		exit(-1);
	}

	_binParams.invWidth = 1. / _binParams.width;
	_binParams.volume = _binParams.width * domain->getGlobalLength(0) * domain->getGlobalLength(2);
}

void DistControl::InitDataStructures()
{
	_nNumValuesScalar = static_cast<uint64_t>(_binParams.count) * static_cast<uint64_t>(_nNumComponents);

	// init offset array
	_nOffsets.resize(_nNumComponents);
	for(uint64_t cid = 0ul; cid < _nNumComponents; ++cid) {
		_nOffsets.at(cid) = cid * static_cast<uint64_t>(_binParams.count);
	}

	// profile midpoint positions
	_dMidpointPositions.resize(_binParams.count);
	for(auto s=0u; s<_binParams.count; ++s)
		_dMidpointPositions.at(s) = (0.5 + s)*_binParams.width;

	// resize
	_nNumMolecules.local.resize(_nNumValuesScalar);
	_nNumMolecules.global.resize(_nNumValuesScalar);
	_dDensityProfile.resize(_nNumValuesScalar);
	_dDensityProfileSmoothed.resize(_nNumValuesScalar);
	_dDensityProfileSmoothedDerivation.resize(_nNumValuesScalar);
	_dForceSum.local.resize(_nNumValuesScalar);
	_dForceSum.global.resize(_nNumValuesScalar);
	_dForceProfile.resize(_nNumValuesScalar);
	_dForceProfileSmoothed.resize(_nNumValuesScalar);

	// init
	for(auto s=0u; s<_nNumValuesScalar; ++s) {
		// number of molecules
		_nNumMolecules.local.at(s) = 0;
		_nNumMolecules.global.at(s) = 0;

		// density profile
		_dDensityProfile.at(s) = 0.;
		_dDensityProfileSmoothed.at(s) = 0.;
		_dDensityProfileSmoothedDerivation.at(s) = 0.;

		// force profile
		_dForceSum.local.at(s) = 0.;
		_dForceSum.global.at(s) = 0.;
		_dForceProfile.at(s) = 0.;
		_dForceProfileSmoothed.at(s) = 0.;
	}
}

void DistControl::PrepareDataStructures()
{


	this->InitDataStructures();
}

void DistControl::SampleProfiles(Molecule* mol)
{
	unsigned int nPosIndex;
	unsigned int nIndexMax = _binParams.count - 1;

	// calc position index
	double dPos = mol->r(1);

	nPosIndex = (unsigned int) floor(dPos * _binParams.invWidth);

	// ignore outer (halo) molecules
	if(nPosIndex > nIndexMax)  // negative values will be ignored to: cast to unsigned int --> high value
		return;

	_nNumMolecules.local.at(nPosIndex)++;
	_dForceSum.local.at(nPosIndex) += mol->F(1);

	unsigned short cid = mol->componentid() + 1;  // cid == 0: sum of all components
	unsigned long nOffset = _nOffsets[cid] + nPosIndex;
	_nNumMolecules.local.at(nOffset)++;
	_dForceSum.local.at(nOffset) += mol->F(1);
}

void DistControl::CalcProfiles()
{
#ifdef ENABLE_MPI

	MPI_Allreduce( _nNumMolecules.local.data(), _nNumMolecules.global.data(), _nNumValuesScalar, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce( _dForceSum.local.data(), _dForceSum.global.data(), _nNumValuesScalar, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
	for(auto s=0u; s<_nNumValuesScalar; ++s)
	{
		_nNumMolecules.global.at(s) = _nNumMolecules.local.at(s);
		_dForceSum.global.at(s) = _dForceSum.local.at(s);
	}
#endif

	// calc profiles
	double dInvertShellVolume = 1. / _binParams.volume;
	double dInvertSampleTimesteps = 1. / ( (double)(_controlFreqs.update) );

	for(auto s=0u; s<_nNumValuesScalar; ++s)
	{
		unsigned long nNumMolecules = _nNumMolecules.global.at(s);
		_dDensityProfile[s] = nNumMolecules * dInvertSampleTimesteps * dInvertShellVolume;

		if(nNumMolecules > 0)
			_dForceProfile[s] = _dForceSum.global.at(s) / ( (double)(nNumMolecules) );  // * dInvertSampleTimesteps; <-- wrong!?? TODO: check
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
//	for(unsigned short cid=0; cid<_nNumComponents; ++cid)
//		std::cout << "_nOffsets[cid] = " << _nOffsets[cid] << std::endl;

	unsigned int nIndexMin = 0;
	unsigned int nIndexMax = 0;
//	std::cout << "_nTargetCompID = " << _nTargetCompID << std::endl;
//	std::cout << "_nOffsets[_nTargetCompID] = " << _nOffsets[_nTargetCompID] << std::endl;
	double* dProfile = _dDensityProfileSmoothedDerivation.data() + _nOffsets.at(_nTargetCompID);
	double dMin = dProfile[0];
	double dMax = dProfile[0];

	// find min/max values and their indexes (positions) in profile
	for(auto s=1u; s<_binParams.count; ++s)
	{
		double dTmp = dProfile[s];
//		if(dTmp < _dVaporDensity)
//			continue;

		if(dTmp > dMax) {
			dMax = dTmp;
			nIndexMax = s;
		}

		if(dTmp < dMin) {
			dMin = dTmp;
			nIndexMin = s;
		}
	}
	_dInterfaceMidLeft  = _dMidpointPositions[nIndexMax];
	_dInterfaceMidRight = _dMidpointPositions[nIndexMin];
}

void DistControl::EstimateInterfaceMidpoint()
{

/*
	// mheinen_2015-03-17 --> DEBUG_TEST_CASE
	for(unsigned int s = 0; s < _binParams.count; ++s)
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
		double dMin = _dDensityProfile.at(0);
		double dMax = _dDensityProfile.at(0);

//		int nIndexMin = 0;
//		int nIndexMax = 0;

		for(auto s=0u; s<_binParams.count/2; ++s)
		{
			 double dDensity = _dDensityProfile.at(s);

			// find min
			if(dDensity < dMin)
			{
				dMin = dDensity;
//				nIndexMin = s;
			}

			// find max
			if(dDensity > dMax)
			{
				dMax = dDensity;
//				nIndexMax = s;
			}
		}

		// dN/dShellWidth = (Nsoll - Nlower) / dInterpol
		//
		// => dInterpol = (Nsoll - Nlower) / dN * dShellWidth
		//
		// dMidPos = dLeftBound + nNumSlabs_lower * dShellWidth + dInterpol

		// 1. ym bestimmen
		double ym = dMin + (dMax - dMin) / 2;

//		std::cout << "dMin = " << dMin << std::endl;
//		std::cout << "dMax = " << dMax << std::endl;
//		std::cout << "ym = " << ym << std::endl;

		// 2. Slab mit gerade niedrigerer Anzahl finden --> von links angefangen ersten Slab größer numMax finden --> Index -1 rechnen
		int nIndexSlabGreater = 0;

		for(auto s=0u; s<_binParams.count; ++s)
		{
			if(_dDensityProfile.at(s) > ym)
			{
				nIndexSlabGreater = s;
				break;
			}
		}

		if(nIndexSlabGreater < 1)
		{
			Log::global_log->error() << "[DistControl] EstimateInterfaceMidpoint(): could not find valid index in density profile" << std::endl;
			return;
		}

//		std::cout << "nIndexSlabGreater = " << nIndexSlabGreater << std::endl;

		//3. Interpolierte Distanz dInterpol berechnen, die von dem Slab mit N < Nsoll ausgehend in Richtung N > Nsoll gegangen werden muss
		double dx12, dx1m, y1, y2, dy12, dy1m;

		y1 = _dDensityProfile.at(nIndexSlabGreater-1);
		y2 = _dDensityProfile.at(nIndexSlabGreater);

		dy1m = ym - y1;

		dy12 = y2 - y1;
		dx12 = _binParams.width;

		// dy12 / dx12 = dy1m / dx1m => dx1m = dy1m / dy12 * dx12
		dx1m = dy1m / dy12 * dx12;

		// 4. Berechnung der Mittelpunktsposition
		// nNumSlabslower == nIndexSlabLower + 1 == nIndexSlabGreater
		double xm, x1;
		double dOuterBoundarySampleZone = 0.;  //this->GetLowerCorner()[1];

		x1 = dOuterBoundarySampleZone + nIndexSlabGreater * _binParams.width - _binParams.width*0.5;

//		std::cout << "x1 = " << x1 << std::endl;
//		std::cout << "dx1m = " << dx1m << std::endl;

		xm = x1 + dx1m;
		_dInterfaceMidLeft = xm;

//		std::cout << "_dInterfaceMidLeft = " << _dInterfaceMidLeft << std::endl;
	}


	// determine right midpoint
	{
		// determine min max
		double dMin = _dDensityProfile.at(0);
		double dMax = _dDensityProfile.at(0);

//		int nIndexMin = 0;
//		int nIndexMax = 0;

		for(auto s=_binParams.count/2; s < _binParams.count; ++s) {
			double dDensity = _dDensityProfile.at(s);

			// find min
			if(dDensity < dMin) {
				dMin = dDensity;
//				nIndexMin = s;
			}

			// find max
			if(dDensity > dMax) {
				dMax = dDensity;
//				nIndexMax = s;
			}
		}

		// dN/dShellWidth = (Nsoll - Nlower) / dInterpol
		//
		// => dInterpol = (Nsoll - Nlower) / dN * dShellWidth
		//
		// dMidPos = dLeftBound + nNumSlabs_lower * dShellWidth + dInterpol

		// 1. ym bestimmen
		double ym = dMin + (dMax - dMin) / 2;

//		std::cout << "ym = " << ym << std::endl;

		// 2. Slab mit gerade niedrigerer Anzahl finden --> von links angefangen ersten Slab größer numMax finden --> Index -1 rechnen
		unsigned int nIndexSlabGreater = 0;

		for(int64_t s= static_cast<int64_t>(_binParams.count) - 1; s >= 0; --s)
		{
			if(_dDensityProfile.at(s) > ym)
			{
				nIndexSlabGreater = s;
				break;
			}
		}

		if(nIndexSlabGreater < 1)
		{
			Log::global_log->error() << "[DistControl] EstimateInterfaceMidpoint(): could not find valid index in density profile" << std::endl;
			return;
		}

		//3. Interpolierte Distanz dInterpol berechnen, die von dem Slab mit N < Nsoll ausgehend in Richtung N > Nsoll gegangen werden muss
		double dx12, dx1m, y1, y2, dy12, dy1m;

		y1 = _dDensityProfile.at(nIndexSlabGreater+1);
		y2 = _dDensityProfile.at(nIndexSlabGreater);

		dy1m = ym - y1;

		dy12 = y2 - y1;
		dx12 = _binParams.width;

		// dy12 / dx12 = dy1m / dx1m => dx1m = dy1m / dy12 * dx12
		dx1m = dy1m / dy12 * dx12;

		// 4. Berechnung der Mittelpunktsposition
		// nNumSlabslower == nIndexSlabLower + 1 == nIndexSlabGreater
		Domain* domain = global_simulation->getDomain();
		double xm, x1;
		double dOuterBoundarySampleZone = domain->getGlobalLength(1);  //this->GetUpperCorner()[1];
		unsigned int numShellsLower = _binParams.count-1 - nIndexSlabGreater;

//		std::cout << "numShellsLower = " << numShellsLower << std::endl;

		x1 = dOuterBoundarySampleZone - numShellsLower * _binParams.width + _binParams.width*0.5;

//		std::cout << "x1 = " << x1 << std::endl;
//		std::cout << "dx1m = " << dx1m << std::endl;

		xm = x1 - dx1m;
		_dInterfaceMidRight = xm;

//		std::cout << "_dInterfaceMidRight = " << _dInterfaceMidRight << std::endl;
	}
}

// init
void DistControl::InitPositions(const double& dInterfaceMidLeft, const double& dInterfaceMidRight)
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

		for (auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit)
		{
			// sample density profile
			this->SampleProfiles(&(*pit));
		}

		// determine interface midpoints and update region positions
		this->UpdatePositions(this->GetUpdateFreq() );
		break;
	case DCIM_MIDPOINT_VALUES:
		break;
	case DCIM_READ_FROM_FILE:
		{
//		std::cout << "_strFilenameInit = " << _strFilenameInit << std::endl;
//		std::cout << "_nRestartTimestep = " << _nRestartTimestep << std::endl;

		std::ifstream filein(_strFilenameInit.c_str(), std::ios::in);

		std::string strLine, strToken;
		std::string strTokens[20];

		while (getline (filein, strLine))
		{
			std::stringstream sstr;
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
		Log::global_log->error() << "[DistControl] Wrong Init Method! Programm exit..." << std::endl;
		exit(-1);
	}

#ifndef NDEBUG
	Log::global_log->error() << "[DistControl] _dInterfaceMidLeft = " << _dInterfaceMidLeft << std::endl;
	Log::global_log->error() << "[DistControl] _dInterfaceMidRight = " << _dInterfaceMidRight << std::endl;
#endif

	// update positions
	this->informObserver();

	// reset local values
	this->ResetLocalValues();
}

void DistControl::UpdatePositions(const uint64_t& simstep)
{
	// update with respect to update frequency
	if(simstep % _controlFreqs.update != 0 || simstep == 0 || simstep == 1)  // TODO init timestep
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
		Log::global_log->error() << "[DistControl] UpdatePositions() Corrupted code!!! Programm exit..." << std::endl;
		exit(-1);
	}

	// update positions
	this->informObserver();

	// reset local values
	this->ResetLocalValues();
}


void DistControl::ResetLocalValues()
{
	for(auto s=0u; s<_nNumValuesScalar; ++s) {
		_nNumMolecules.local.at(s) = 0;
		_dForceSum.local.at(s) = 0.;
	}
}

void DistControl::WriteData(const uint64_t& simstep)
{
	// domain decomposition
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();

	// write out data
	std::stringstream outputstream;

	// write data
	if(simstep % _controlFreqs.write.data == 0)
	{
		// calc global values
//		this->CalcGlobalValues(domainDecomp);

#ifdef ENABLE_MPI
	int rank = domainDecomp.getRank();
	// int numprocs = domainDecomp.getNumProcs();
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
			CuboidRegionObs* region = dynamic_cast<CuboidRegionObs*>(*it);
			if(nullptr != region) {
				outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << region->GetLowerCorner()[1];
				outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << region->GetUpperCorner()[1];
			}
			Mirror* mirror = dynamic_cast<Mirror*>(*it);
			if(nullptr != mirror) {
				outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << mirror->getPosition();
			}
		}

		outputstream << std::endl;

		std::ofstream fileout(_strFilename.c_str(), std::ios::out|std::ios::app);
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
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();

	// write header
	std::stringstream outputstream;

#ifdef ENABLE_MPI
	int rank = domainDecomp.getRank();
	// int numprocs = domainDecomp.getNumProcs();
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
		CuboidRegionObs* region = dynamic_cast<CuboidRegionObs*>(*it);
		if(nullptr != region) {
			outputstream << std::setw(20) << region->GetParent()->getShortName() << "l" << "[" << region->GetID() << "]";
			outputstream << std::setw(20) << region->GetParent()->getShortName() << "r" << "[" << region->GetID() << "]";
		}

		ControlInstance* ctrlInst = dynamic_cast<ControlInstance*>(*it);
		if(nullptr != ctrlInst)
			outputstream << std::setw (24) << ctrlInst->getShortName();
	}
	outputstream << std::endl;

	std::ofstream fileout(_strFilename.c_str(), std::ios::out);
	fileout << outputstream.str();
	fileout.close();

#ifdef ENABLE_MPI
	}
#endif
}


void DistControl::WriteDataProfiles(const uint64_t& simstep)
{
	if(simstep % _controlFreqs.write.profiles != 0)
		return;

	// domain decomposition
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();

	std::stringstream outputstream;
	std::stringstream filenamestream;

	filenamestream << _strFilenameProfilesPrefix << "_TS";
	filenamestream.fill('0');
	filenamestream.width(9);
	filenamestream << std::right << simstep;
	filenamestream << ".dat";

#ifdef ENABLE_MPI
	int rank = domainDecomp.getRank();
	// int numprocs = domainDecomp.getNumProcs();
	if (rank != 0)
		return;
#endif

	// write header
	outputstream << "                   coord";
	for(auto cid=0; cid<_nNumComponents; ++cid)
	{
		outputstream << "                  rho[" << cid << "]";
		outputstream << "           rho_smooth[" << cid << "]";
		outputstream << "              drho/dy[" << cid << "]";
		outputstream << "                   Fy[" << cid << "]";
		outputstream << "            Fy_smooth[" << cid << "]";
	}
	outputstream << std::endl;
	// write data
	for(auto s=0u; s<_binParams.count; ++s)
	{
		outputstream << FORMAT_SCI_MAX_DIGITS << _dMidpointPositions.at(s);
		for(auto cid=0; cid<_nNumComponents; ++cid)
		{
			uint64_t nIndex = _nOffsets[cid]+s;
			outputstream << FORMAT_SCI_MAX_DIGITS << _dDensityProfile.at(nIndex);
			outputstream << FORMAT_SCI_MAX_DIGITS << _dDensityProfileSmoothed.at(nIndex);
			outputstream << FORMAT_SCI_MAX_DIGITS << _dDensityProfileSmoothedDerivation.at(nIndex);
			outputstream << FORMAT_SCI_MAX_DIGITS << _dForceProfile.at(nIndex);
			outputstream << FORMAT_SCI_MAX_DIGITS << _dForceProfileSmoothed.at(nIndex);
		}
		outputstream << std::endl;
	}
	std::ofstream fileout(filenamestream.str().c_str(), std::ios::out|std::ios::out);
	fileout << outputstream.str();
	fileout.close();
}


void DistControl::Init(ParticleContainer* particleContainer)
{
	this->PrepareSubdivision();
	this->PrepareDataStructures();
}

void DistControl::registerObserver(ObserverBase* observer)
{
	_observer.push_back(observer);
}

void DistControl::unregisterObserver(ObserverBase* observer)
{
	std::vector<ObserverBase*>::iterator it;
	for(it=_observer.begin(); it!=_observer.end(); ++it)
		if(*it == observer)
			_observer.erase(it);
}

void DistControl::informObserver()
{
	for(auto it:_observer)
		it->update(this);
}

void DistControl::SmoothProfile(double* dData, double* dSmoothData, const uint64_t& nNumVals, const uint32_t& nNeighbourVals)
{
	uint32_t li = 0;
	uint32_t ri = nNeighbourVals;

	for(auto s=0u; s<nNumVals; ++s)
	{
		double dSum = 0.;
		uint16_t nNumValsSmooth = 0;
		for(auto ti = s-li; ti <= s+ri; ++ti) {
			dSum += dData[ti];
			nNumValsSmooth++;
		}
		double dInvNumValsSmooth = 1. / ( (double)(nNumValsSmooth) );
		dSmoothData[s] = dSum * dInvNumValsSmooth;

		if(li < nNeighbourVals)
			li++;

		if(s >= (nNumVals-1 - nNeighbourVals) )
			ri--;
	}
}

void DistControl::SmoothProfiles(const uint32_t& nNeighbourVals)
{
	for(auto cid=0u; cid<_nNumComponents; ++cid) {
		uint64_t nOffset = _nOffsets.at(cid);
		SmoothProfile(_dDensityProfile.data()+nOffset, _dDensityProfileSmoothed.data()+nOffset, _binParams.count, nNeighbourVals);
		SmoothProfile(_dForceProfile.data()+nOffset, _dForceProfileSmoothed.data()+nOffset, _binParams.count, nNeighbourVals);
	}
}

void DistControl::DerivateProfile(double* dDataX, double* dDataY, double* dDerivDataY, const uint64_t& nNumVals, const uint32_t& nNeighbourVals)
{
	uint32_t li = 0;
	uint32_t ri = nNeighbourVals;

/*
 * derivate by two points
 *
	for(unsigned int s = 0; s < _binParams.count; ++s)
	{
		double dy = _dDensityProfileSmoothed[s+r] - _dDensityProfileSmoothed[s-l];
		double dx = (r+l)*_binParams.width;

		_dDensityProfileSmoothedDerivation[s] = dy / dx;

		if(l < nNumNeighbourVals)
			l++;

		if(s >= (_binParams.count-1 - nNumNeighbourVals) )
			r--;
	}
*/

	std::vector<double> x;
	std::vector<double> y;

	for(auto s=0u; s<nNumVals; ++s) {
		x.clear();
		y.clear();

		for(auto i=s-li; i<=s+ri; ++i) {
			x.push_back(dDataX[i] );
			y.push_back(dDataY[i] );
		}
		double beta1, beta2;
		LinearRegression(x, y, beta1, beta2);

		dDerivDataY[s] = beta2;

		if(li < nNeighbourVals)
			li++;

		if(s >= (_binParams.count-1 - nNeighbourVals) )
			ri--;
	}
}

void DistControl::DerivateProfiles(const unsigned int& nNeighbourVals)
{
	for(auto cid=0; cid<_nNumComponents; ++cid) {
		uint64_t nOffset = _nOffsets.at(cid);
		this->DerivateProfile(_dMidpointPositions.data(), _dDensityProfileSmoothed.data()+nOffset, _dDensityProfileSmoothedDerivation.data()+nOffset,
				_binParams.count, nNeighbourVals);
	}
}










