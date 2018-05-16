/*
 * RegionSampling.cpp
 *
 *  Created on: 18.03.2015
 *      Author: mheinen
 */

#include "RegionSampling.h"
#include "Domain.h"
#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "utils/FileUtils.h"
#include "utils/xmlfileUnits.h"
#include "NEMD/DistControl.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <cmath>
#include <cstdint>
#include <algorithm>  // std::fill
#include <array>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;


// init static ID --> instance counting
unsigned short SampleRegion::_nStaticID = 0;

void SampleRegion::get_v(double* q, Molecule* mol)
{
	for(uint16_t d=0; d<3; ++d)
		q[d] = mol->v(d);
}

void SampleRegion::get_F(double* q, Molecule* mol)
{
	for(uint16_t d=0; d<3; ++d)
		q[d] = mol->F(d);
}

void SampleRegion::get_v2(double& q, Molecule* mol)
{
	q = mol->v2();
}

void SampleRegion::get_F2(double& q, Molecule* mol)
{
	q = mol->F2();
}

SampleRegion::SampleRegion( RegionSampling* parent, double dLowerCorner[3], double dUpperCorner[3] )
	: CuboidRegionObs(parent, dLowerCorner, dUpperCorner),
	_bDiscretisationDoneProfiles(false),
	_SamplingEnabledProfiles(false),
	_nSubdivisionOptProfiles(SDOPT_UNKNOWN),
	_bDiscretisationDoneVDF(false),
	_SamplingEnabledVDF(false),
	_nSubdivisionOptVDF(SDOPT_UNKNOWN),
	_fnamePrefixVDF("VDF"),
	_bDiscretisationDoneFieldYR(false),
	_SamplingEnabledFieldYR(false),
	_nSubdivisionOptFieldYR_Y(SDOPT_UNKNOWN),
	_nSubdivisionOptFieldYR_R(SDOPT_UNKNOWN)
{
	_nID = ++_nStaticID;
	_nSubdivisionOpt = SDOPT_UNKNOWN;

	_numComponents = global_simulation->getEnsemble()->getComponents()->size()+1;  // cid == 0: all components
	_vecMass.resize(_numComponents);
	_vecMass.at(0) = 0.;
	for(uint8_t cid=1; cid<_numComponents; ++cid)
	{
		_vecMass.at(cid) = global_simulation->getEnsemble()->getComponent(cid-1)->m();
//		cout << "cid = " << cid << ": mass = " << _vecMass.at(cid) << endl;
	}

	// Init component specific parameters for VDF sampling
	this->InitComponentSpecificParamsVDF();

	_fptr = get_v;
	_f2ptr = get_v2;
}

SampleRegion::~SampleRegion()
{
	// profile sampling
	if(true == _SamplingEnabledProfiles)
	{
		delete [] _dBinMidpointsProfiles;
		// Offsets
		for(uint8_t dir = 0; dir<3; ++dir)
			delete [] _nOffsetScalar[dir];

		for(uint8_t dim = 0; dim<3; ++dim)
		{
			for(uint8_t dir = 0; dir<3; ++dir)
				delete [] _nOffsetVector[dim][dir];
		}
		// Scalar quantities
		delete [] _nNumMoleculesLocal;
		delete [] _nNumMoleculesGlobal;
		delete [] _nRotDOFLocal;
		delete [] _nRotDOFGlobal;
		delete [] _d2EkinRotLocal;
		delete [] _d2EkinRotGlobal;
		// output profiles
		delete [] _dDensity;
		delete [] _d2EkinTotal;
		delete [] _d2EkinTrans;
		delete [] _d2EkinDrift;
		delete [] _d2EkinRot;
		delete [] _d2EkinT;
		delete [] _dTemperature;
		delete [] _dTemperatureTrans;
		delete [] _dTemperatureRot;
		// Vector quantities
		delete [] _dVelocityLocal;
		delete [] _dVelocityGlobal;
		delete [] _dSquaredVelocityLocal;
		delete [] _dSquaredVelocityGlobal;
		delete [] _dForceLocal;
		delete [] _dForceGlobal;
		// output profiles
		delete [] _dForce;
		delete [] _dDriftVelocity;
		delete [] _d2EkinTransComp;
		delete [] _d2EkinDriftComp;
		delete [] _dTemperatureComp;
	}

	// VDF sampling
	if(true == _SamplingEnabledVDF)
	{
		// local
		delete [] _VDF_pjy_abs_local;
		delete [] _VDF_pjy_pvx_local;
		delete [] _VDF_pjy_pvy_local;
		delete [] _VDF_pjy_pvz_local;
		delete [] _VDF_pjy_nvx_local;
		delete [] _VDF_pjy_nvz_local;

		delete [] _VDF_njy_abs_local;
		delete [] _VDF_njy_pvx_local;
		delete [] _VDF_njy_pvz_local;
		delete [] _VDF_njy_nvx_local;
		delete [] _VDF_njy_nvy_local;
		delete [] _VDF_njy_nvz_local;

		// global
		delete [] _VDF_pjy_abs_global;
		delete [] _VDF_pjy_pvx_global;
		delete [] _VDF_pjy_pvy_global;
		delete [] _VDF_pjy_pvz_global;
		delete [] _VDF_pjy_nvx_global;
		delete [] _VDF_pjy_nvz_global;

		delete [] _VDF_njy_abs_global;
		delete [] _VDF_njy_pvx_global;
		delete [] _VDF_njy_pvz_global;
		delete [] _VDF_njy_nvx_global;
		delete [] _VDF_njy_nvy_global;
		delete [] _VDF_njy_nvz_global;
	}

	// fieldYR sampling
	if(true == _SamplingEnabledFieldYR)
	{
		delete [] _dShellVolumesFieldYR;
		delete [] _dBinMidpointsFieldYR;
		delete [] _dShellMidpointsFieldYR;
		// Offsets
		for(uint8_t dim=0; dim<3; ++dim)
		{
			for(uint8_t sec=0; sec<3; ++sec)
			{
				for(uint32_t cid = 0; cid<_numComponents; ++cid)
					delete [] _nOffsetFieldYR[dim][sec][cid];
				delete [] _nOffsetFieldYR[dim][sec];
			}
			delete [] _nOffsetFieldYR[dim];
		}
		delete [] _nOffsetFieldYR;

		// Scalar quantities
		delete [] _nNumMoleculesFieldYRLocal;
		delete [] _nNumMoleculesFieldYRGlobal;
		// output profiles
		delete [] _dDensityFieldYR;
		delete [] _dInvShellVolumesFieldYR;
	}
}

void SampleRegion::InitComponentSpecificParamsVDF()
{
	_vecComponentSpecificParamsVDF.resize(_numComponents);
	for(auto&& csp : _vecComponentSpecificParamsVDF)
	{
		csp.bSamplingEnabled = false;
		csp.dVeloMax = 0.0;
		csp.dVelocityClassWidth = 0.0;
		csp.dInvVelocityClassWidth = 0.0;
		csp.numVelocityClasses = 0;
		csp.numVals = 0;
		csp.nOffsetDataStructure = 0;
		csp.dDiscreteVelocityValues.clear();
	}
}

void SampleRegion::ShowComponentSpecificParamsVDF()
{
	global_log->info() << ">>>> ComponentSpecificParamsVDF" << endl;
	global_log->info() << "-----------------------------------------------------" << endl;
	uint32_t cid = 0;
	for(auto&& csp : _vecComponentSpecificParamsVDF)
	{
		global_log->info() << "bSamplingEnabled = "       << csp.bSamplingEnabled << endl;
		global_log->info() << "dVeloMax = "               << csp.dVeloMax << endl;
		global_log->info() << "dVelocityClassWidth = "    << csp.dVelocityClassWidth << endl;
		global_log->info() << "dInvVelocityClassWidth = " << csp.dInvVelocityClassWidth << endl;
		global_log->info() << "numVelocityClasses = "     << csp.numVelocityClasses << endl;
		global_log->info() << "numVals = "                << csp.numVals << endl;
		global_log->info() << "nOffsetDataStructure = "   << csp.nOffsetDataStructure << endl;
		global_log->info() << "-----------------------------------------------------" << endl;
	}
	global_log->info() << "<<<< ComponentSpecificParamsVDF" << endl;
}

void SampleRegion::readXML(XMLfileUnits& xmlconfig)
{
	// sampling modules
	uint32_t nSamplingModuleID = 0;
	uint32_t numSamplingModules = 0;
	XMLfile::Query query = xmlconfig.query("sampling");
	numSamplingModules = query.card();
	global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]: Number of sampling modules: " << numSamplingModules << endl;
	if(numSamplingModules < 1) {
		global_log->error() << "RegionSampling->region["<<this->GetID()-1<<"]: No sampling module parameters specified. Program exit ..." << endl;
		Simulation::exit(-1);
	}

	XMLfile::Query::const_iterator outputSamplingIter;
	for( outputSamplingIter = query.begin(); outputSamplingIter != query.end(); outputSamplingIter++ )
	{
		xmlconfig.changecurrentnode( outputSamplingIter );
		std::string strSamplingModuleType = "unknown";
		xmlconfig.getNodeValue("@type", strSamplingModuleType);

		if("profiles" == strSamplingModuleType)
		{
			// enable profile sampling
			_SamplingEnabledProfiles = true;
			// control
			xmlconfig.getNodeValue("control/start", _initSamplingProfiles);
			xmlconfig.getNodeValue("control/frequency", _writeFrequencyProfiles);
			xmlconfig.getNodeValue("control/stop", _stopSamplingProfiles);
			this->SetParamProfiles();
			global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Start sampling from simstep: " << _initSamplingProfiles << endl;
			global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Sample with frequency: " << _writeFrequencyProfiles << endl;
			global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Stop sampling at simstep: " << _stopSamplingProfiles << endl;

			// subdivision of region
			std::string strSubdivisionType;
			if( !xmlconfig.getNodeValue("subdivision@type", strSubdivisionType) )
			{
				global_log->error() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Missing attribute subdivision@type! Program exit..." << endl;
				Simulation::exit(-1);
			}
			if("number" == strSubdivisionType)
			{
				unsigned int nNumSlabs = 0;
				if( !xmlconfig.getNodeValue("subdivision/number", nNumSlabs) )
				{
					global_log->error() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Missing element subdivision/number! Program exit..." << endl;
					Simulation::exit(-1);
				}
				else
				{
					this->SetSubdivisionProfiles(nNumSlabs);
					global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): subdivision by '"<< strSubdivisionType << "': " << nNumSlabs << endl;
				}
			}
			else if("width" == strSubdivisionType)
			{
				double dSlabWidth = 0.;
				if( !xmlconfig.getNodeValue("subdivision/width", dSlabWidth) )
				{
					global_log->error() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Missing element subdivision/width! Program exit..." << endl;
					Simulation::exit(-1);
				}
				else
				{
					this->SetSubdivisionProfiles(dSlabWidth);
					global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): subdivision by '"<< strSubdivisionType << "': " << dSlabWidth << endl;
				}
			}
			else
			{
				global_log->error() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Wrong attribute subdivision@type. Expected type=\"number|width\"! Program exit..." << endl;
				Simulation::exit(-1);
			}
		}
		else if("VDF" == strSamplingModuleType || "FDF" == strSamplingModuleType)
		{
			if("VDF" == strSamplingModuleType)
			{
				_fptr = get_v;
				_f2ptr = get_v2;
				_fnamePrefixVDF = "VDF";
			}
			else if("FDF" == strSamplingModuleType)
			{
				_fptr = get_F;
				_f2ptr = get_F2;
				_fnamePrefixVDF = "FDF";
			}

			// enable VDF sampling
			_SamplingEnabledVDF = true;
			// control
			xmlconfig.getNodeValue("control/start", _initSamplingVDF);
			xmlconfig.getNodeValue("control/frequency", _writeFrequencyVDF);
			xmlconfig.getNodeValue("control/stop", _stopSamplingVDF);

			// velocity dicretization
			{
				string oldpath = xmlconfig.getcurrentnodepath();
				xmlconfig.changecurrentnode("discretizations");
				uint32_t numDiscretizations = 0;
				XMLfile::Query query_vd = xmlconfig.query("discretization");
				numDiscretizations = query_vd.card();
				global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]: Number of velocity discretizations (components for that VDF should be sampled): " << numDiscretizations << endl;
				if(numDiscretizations < 1) {
					global_log->error() << "RegionSampling->region["<<this->GetID()-1<<"]: No velocity discretizations specified for VDF sampling. Program exit ..." << endl;
					Simulation::exit(-1);
				}
				XMLfile::Query::const_iterator nodeIter;
				for( nodeIter = query_vd.begin(); nodeIter != query_vd.end(); nodeIter++ )
				{
					xmlconfig.changecurrentnode(nodeIter);
					uint32_t cid = 0;
					if(false == (xmlconfig.getNodeValue("@cid", cid) && cid < _numComponents) )
						global_log->error() << "RegionSampling->region["<<this->GetID()-1<<"]: VDF velocity discretization corrupted. Program exit ..." << endl;

					ComponentSpecificParamsVDF& csp = _vecComponentSpecificParamsVDF.at(cid);
					csp.bSamplingEnabled = true;
					xmlconfig.getNodeValue("numclasses", csp.numVelocityClasses);
					xmlconfig.getNodeValue("maxvalue", csp.dVeloMax);
				}
				xmlconfig.changecurrentnode(oldpath);
			}

			global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Start sampling from simstep: " << _initSamplingVDF << endl;
			global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Sample with frequency: " << _writeFrequencyVDF << endl;
			global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Stop sampling at simstep: " << _stopSamplingVDF << endl;

			// subdivision of region
			std::string strSubdivisionType;
			if( !xmlconfig.getNodeValue("subdivision@type", strSubdivisionType) )
			{
				global_log->error() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Missing attribute subdivision@type! Program exit..." << endl;
				Simulation::exit(-1);
			}
			if("number" == strSubdivisionType)
			{
				unsigned int nNumSlabs = 0;
				if( !xmlconfig.getNodeValue("subdivision/number", nNumSlabs) )
				{
					global_log->error() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Missing element subdivision/number! Program exit..." << endl;
					Simulation::exit(-1);
				}
				else
				{
					this->SetSubdivisionVDF(nNumSlabs);
					global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): subdivision by '"<< strSubdivisionType << "': " << nNumSlabs << endl;
				}
			}
			else if("width" == strSubdivisionType)
			{
				double dSlabWidth = 0.;
				if( !xmlconfig.getNodeValue("subdivision/width", dSlabWidth) )
				{
					global_log->error() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Missing element subdivision/width! Program exit..." << endl;
					Simulation::exit(-1);
				}
				else
				{
					this->SetSubdivisionVDF(dSlabWidth);
					global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): subdivision by '"<< strSubdivisionType << "': " << dSlabWidth << endl;
				}
			}
			else
			{
				global_log->error() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Wrong attribute subdivision@type. Expected type=\"number|width\"! Program exit..." << endl;
				Simulation::exit(-1);
			}
		}
		else if("fieldYR" == strSamplingModuleType)
		{
			// enable fieldYR sampling
			_SamplingEnabledFieldYR = true;
			bool bInputIsValid = true;

			// output file
			bInputIsValid = true;
			std::string strFileTypeFieldYR;
			_nFileTypeFieldYR = RSFT_BINARY;
			_strFilePrefixFieldYR = "fieldYR";
			bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("outputfile/@type", strFileTypeFieldYR);
			if(true == bInputIsValid)
			{
				if("ASCII" == strFileTypeFieldYR)
					_nFileTypeFieldYR = RSFT_ASCII;
				else if("binary" == strFileTypeFieldYR)
					_nFileTypeFieldYR = RSFT_BINARY;
				else
					bInputIsValid = false;
			}
			bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("outputfile/prefix", _strFilePrefixFieldYR);
			if(true == bInputIsValid)
				global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Writing y-r-field data to '" << strFileTypeFieldYR << "'-files "
				"staring with prefix: " << _strFilePrefixFieldYR << endl;
			else
			{
				global_log->error() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Parameters of element: 'outputfile' corrupted! Program exit..." << endl;
				Simulation::exit(-1);
			}

			// control
			bInputIsValid = true;
			bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("control/start", _initSamplingFieldYR);
			bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("control/frequency", _writeFrequencyFieldYR);
			bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("control/stop", _stopSamplingFieldYR);
			if(true == bInputIsValid)
			{
				global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Start sampling from simstep: " << _initSamplingFieldYR << endl;
				global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Sample with frequency: " << _writeFrequencyFieldYR << endl;
				global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Stop sampling at simstep: " << _stopSamplingFieldYR << endl;
			}
			else
			{
				global_log->error() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Parameters of element: 'control' corrupted! Program exit..." << endl;
				Simulation::exit(-1);
			}

			// subdivision of region
			uint32_t numSubdivisions = 0;
			XMLfile::Query query_sd = xmlconfig.query("subdivision");
			numSubdivisions = query_sd.card();
			if(numSubdivisions != 2) {
				global_log->error() << "RegionSampling->region["<<this->GetID()-1<<"]: Found " << numSubdivisions << " 'subdivision' elements, "
						"expected: 2. Program exit ..." << endl;
				Simulation::exit(-1);
			}
			string oldpath = xmlconfig.getcurrentnodepath();
			XMLfile::Query::const_iterator outputSubdivisionIter;
			for( outputSubdivisionIter = query_sd.begin(); outputSubdivisionIter != query_sd.end(); outputSubdivisionIter++ )
			{
				xmlconfig.changecurrentnode( outputSubdivisionIter );
				bInputIsValid = true;
				std::string strSubdivisionType;
				std::string strSubdivisionDim;
				bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("@dim", strSubdivisionDim);
				if("y" == strSubdivisionDim)
				{
					bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("@type", strSubdivisionType);
					if("number" == strSubdivisionType) {
						_nSubdivisionOptFieldYR_Y = SDOPT_BY_NUM_SLABS;
						bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("number", _nNumBinsFieldYR);
						global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): subdivision 'y' by '"<< strSubdivisionType << "': " << _nNumBinsFieldYR << endl;
					}
					else if("width" == strSubdivisionType) {
						_nSubdivisionOptFieldYR_Y = SDOPT_BY_SLAB_WIDTH;
						bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("width", _dBinWidthInitFieldYR);
						global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): subdivision 'y' by '"<< strSubdivisionType << "': " << _dBinWidthInitFieldYR << endl;
					}
					else
						bInputIsValid = false;
				}
				else if("r" == strSubdivisionDim)
				{
					bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("@type", strSubdivisionType);
					if("number" == strSubdivisionType) {
						_nSubdivisionOptFieldYR_R = SDOPT_BY_NUM_SLABS;
						bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("number", _nNumShellsFieldYR);
						global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): subdivision 'r' by '"<< strSubdivisionType << "': " << _nNumShellsFieldYR << endl;
					}
					else if("width" == strSubdivisionType) {
						_nSubdivisionOptFieldYR_R = SDOPT_BY_SLAB_WIDTH;
						bInputIsValid = bInputIsValid && xmlconfig.getNodeValue("width", _dShellWidthInitFieldYR);
						global_log->info() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): subdivision 'r' by '"<< strSubdivisionType << "': " << _dShellWidthInitFieldYR << endl;
					}
					else
						bInputIsValid = false;
				}
				else
					bInputIsValid = false;

				if(false == bInputIsValid)
				{
					global_log->error() << "RegionSampling->region["<<this->GetID()-1<<"]->sampling('"<<strSamplingModuleType<<"'): Parameters for elements: 'subdivision' corrupted! Program exit..." << endl;
					Simulation::exit(-1);
				}
			}  // for( outputSubdivisionIter = query.begin(); outputSubdivisionIter; outputSubdivisionIter++ )
			xmlconfig.changecurrentnode(oldpath);
		}
		else
		{
			global_log->error() << "RegionSampling->region["<<this->GetID()-1<<"]: Wrong attribute 'sampling@type', expected type='profiles|VDF|fieldYR'! Program exit..." << endl;
			Simulation::exit(-1);
		}
	}  // for( outputSamplingIter = query.begin(); outputSamplingIter; outputSamplingIter++ )
}

void SampleRegion::PrepareSubdivisionProfiles()
{
	if(false == _SamplingEnabledProfiles)
		return;

	double dWidth = this->GetWidth(1);

	switch(_nSubdivisionOptProfiles)
	{
	case SDOPT_BY_NUM_SLABS:
		_dBinWidthProfilesInit = this->GetWidth(1) / ( (double)(_nNumBinsProfiles) );
		_dBinWidthProfiles = _dBinWidthProfilesInit;
		break;
	case SDOPT_BY_SLAB_WIDTH:
		_nNumBinsProfiles = round(dWidth / _dBinWidthProfilesInit);
		_dBinWidthProfiles = dWidth / ( (double)(_nNumBinsProfiles) );
		break;
	case SDOPT_UNKNOWN:
	default:
		global_log->error() << "tec::ControlRegion::PrepareSubdivisionProfiles(): Unknown subdivision type! Program exit..." << endl;
		exit(-1);
	}
}

void SampleRegion::PrepareSubdivisionVDF()
{
	if(false == _SamplingEnabledVDF)
		return;

	double dWidth = this->GetWidth(1);

	switch(_nSubdivisionOptVDF)
	{
	case SDOPT_BY_NUM_SLABS:
		_dBinWidthVDFInit = this->GetWidth(1) / ( (double)(_numBinsVDF) );
		_dBinWidthVDF = _dBinWidthVDFInit;
		break;
	case SDOPT_BY_SLAB_WIDTH:
		_numBinsVDF = round(dWidth / _dBinWidthVDFInit);
		_dBinWidthVDF = dWidth / ( (double)(_numBinsVDF) );
		break;
	case SDOPT_UNKNOWN:
	default:
		global_log->error() << "ERROR in SampleRegion::PrepareSubdivisionVDF(): Unknown subdivision type! Program exit..." << endl;
		exit(-1);
	}
}

void SampleRegion::PrepareSubdivisionFieldYR()
{
	if(false == _SamplingEnabledFieldYR)
		return;

	double dWidth = this->GetWidth(1);

	switch(_nSubdivisionOptFieldYR_Y)
	{
	case SDOPT_BY_NUM_SLABS:
		_dBinWidthInitFieldYR = dWidth / ( (double)(_nNumBinsFieldYR) );
		_dBinWidthFieldYR = _dBinWidthInitFieldYR;
		break;
	case SDOPT_BY_SLAB_WIDTH:
		_nNumBinsFieldYR = round(dWidth / _dBinWidthInitFieldYR);
		_dBinWidthFieldYR = dWidth / ( (double)(_nNumBinsFieldYR) );
		break;
	case SDOPT_UNKNOWN:
	default:
		global_log->error() << "SampleRegion::PrepareSubdivisionFieldYR(): Unknown subdivision type! Program exit..." << endl;
		Simulation::exit(-1);
	}

	dWidth = (this->GetWidth(0) < this->GetWidth(2) ) ? this->GetWidth(0) : this->GetWidth(2);
	dWidth *= 0.5;

	switch(_nSubdivisionOptFieldYR_R)
	{
	case SDOPT_BY_NUM_SLABS:
		_dShellWidthInitFieldYR = dWidth / ( (double)(_nNumShellsFieldYR) );
		_dShellWidthFieldYR = _dShellWidthInitFieldYR;
		break;
	case SDOPT_BY_SLAB_WIDTH:
		_nNumShellsFieldYR = round(dWidth / _dShellWidthInitFieldYR);
		_dShellWidthFieldYR = dWidth / ( (double)(_nNumShellsFieldYR) );
		break;
	case SDOPT_UNKNOWN:
	default:
		global_log->error() << "SampleRegion::PrepareSubdivisionFieldYR(): Unknown subdivision type! Program exit..." << endl;
		Simulation::exit(-1);
	}
}

void SampleRegion::InitSamplingProfiles(int nDimension)
{
	if(false == _SamplingEnabledProfiles)
		return;

	// Bin width
	double dNumBinsProfiles = (double) _nNumBinsProfiles;
	_dBinWidthProfilesInit = this->GetWidth(nDimension) / dNumBinsProfiles;
	_dBinWidthProfiles = _dBinWidthProfilesInit;

	// Bin volume
	double dArea;
	Domain* domain = this->GetParent()->GetDomain();

	switch(nDimension)
	{
	case RS_DIMENSION_X:
		dArea = domain->getGlobalLength(1) * domain->getGlobalLength(2);
		break;

	case RS_DIMENSION_Y:
		dArea = domain->getGlobalLength(0) * domain->getGlobalLength(2);
		break;

	case RS_DIMENSION_Z:
		dArea = domain->getGlobalLength(0) * domain->getGlobalLength(1);
		break;
	default:
		dArea = domain->getGlobalLength(0) * domain->getGlobalLength(2);
	}
	_dBinVolumeProfiles = _dBinWidthProfiles * dArea;
	_dInvertBinVolumeProfiles = 1. / _dBinVolumeProfiles;
	_dInvertBinVolSamplesProfiles = _dInvertBinVolumeProfiles * _dInvertNumSamplesProfiles;


	// discrete values: Bin midpoints, velocity values
	_dBinMidpointsProfiles = new double[_nNumBinsProfiles];

	_nNumValsScalar = _nNumBinsProfiles * _numComponents * 3;  // * 3: directions: all(+/-) | only (+) | only (-)
	_nNumValsVector = _nNumValsScalar * 3;                        // * 3: x, y, z-component

#ifndef NDEBUG
	cout << "_nNumBinsProfiles = " << _nNumBinsProfiles << endl;
	cout << "_numComponents    = " << _numComponents   << endl;
	cout << "_nNumValsScalar   = " << _nNumValsScalar   << endl;
	cout << "_nNumValsVector   = " << _nNumValsVector   << endl;
#endif

	// Offsets
	_nOffsetScalar = new unsigned long*[3];
	for(unsigned int dir = 0; dir<3; ++dir)
		_nOffsetScalar[dir] = new unsigned long[_numComponents];

	_nOffsetVector = new unsigned long**[3];
	for(unsigned int dim = 0; dim<3; ++dim)
	{
		_nOffsetVector[dim] = new unsigned long*[3];
		for(unsigned int dir = 0; dir<3; ++dir)
			_nOffsetVector[dim][dir] = new unsigned long[_numComponents];
	}

	unsigned long nOffset;

	// Scalar quantities
	nOffset = 0;
	for(unsigned int dir = 0; dir<3; ++dir){
		for(unsigned int cid = 0; cid<_numComponents; ++cid){
			_nOffsetScalar[dir][cid] = nOffset;
			nOffset += _nNumBinsProfiles;
		}
	}
	// Vector quantities
	nOffset = 0;
	for(unsigned int dim = 0; dim<3; ++dim){
		for(unsigned int dir = 0; dir<3; ++dir){
			for(unsigned int cid = 0; cid<_numComponents; ++cid){
				_nOffsetVector[dim][dir][cid] = nOffset;
				nOffset += _nNumBinsProfiles;
			}
		}
	}

	// Scalar quantities
	// [direction all|+|-][component][position]
	_nNumMoleculesLocal  = new unsigned long[_nNumValsScalar];
	_nNumMoleculesGlobal = new unsigned long[_nNumValsScalar];
	_nRotDOFLocal  = new unsigned long[_nNumValsScalar];
	_nRotDOFGlobal = new unsigned long[_nNumValsScalar];
	_d2EkinRotLocal  = new double[_nNumValsScalar];
	_d2EkinRotGlobal = new double[_nNumValsScalar];

	// output profiles
	_dDensity = new double[_nNumValsScalar];
	_d2EkinTotal = new double[_nNumValsScalar];
	_d2EkinTrans = new double[_nNumValsScalar];
	_d2EkinDrift = new double[_nNumValsScalar];
	_d2EkinRot   = new double[_nNumValsScalar];
	_d2EkinT     = new double[_nNumValsScalar];
	_dTemperature      = new double[_nNumValsScalar];
	_dTemperatureTrans = new double[_nNumValsScalar];
	_dTemperatureRot   = new double[_nNumValsScalar];

	// Vector quantities
	// [direction all|+|-][component][position][dimension x|y|z]
	_dVelocityLocal  = new double[_nNumValsVector];
	_dVelocityGlobal = new double[_nNumValsVector];
	_dSquaredVelocityLocal  = new double[_nNumValsVector];
	_dSquaredVelocityGlobal = new double[_nNumValsVector];
	_dForceLocal  = new double[_nNumValsVector];
	_dForceGlobal = new double[_nNumValsVector];

	// output profiles
	_dForce = new double[_nNumValsVector];
	_dDriftVelocity   = new double[_nNumValsVector];
	_d2EkinTransComp   = new double[_nNumValsVector];
	_d2EkinDriftComp   = new double[_nNumValsVector];
	_dTemperatureComp = new double[_nNumValsVector];

	// init sampling data structures
	this->ResetLocalValuesProfiles();

	// discretisation
	this->DoDiscretisationProfiles(RS_DIMENSION_Y);
}


void SampleRegion::InitSamplingVDF(int nDimension)
{
	if(false == _SamplingEnabledVDF)
		return;

	// Init component specific parameters
	uint32_t nOffsetDataStructure = 0;
	for(auto&& csp : _vecComponentSpecificParamsVDF)
	{
		if(false == csp.bSamplingEnabled)
			continue;
		csp.nOffsetDataStructure = nOffsetDataStructure;
		csp.dVelocityClassWidth = csp.dVeloMax / ((double)(csp.numVelocityClasses));
		csp.dInvVelocityClassWidth = 1./csp.dVelocityClassWidth;
		csp.numVals = _numBinsVDF * csp.numVelocityClasses;
		nOffsetDataStructure += csp.numVals;
		csp.dDiscreteVelocityValues.resize(csp.numVelocityClasses);
		std::fill (csp.dDiscreteVelocityValues.begin(),csp.dDiscreteVelocityValues.end(),0);
	}
	_numValsVDF = nOffsetDataStructure;

	// discrete values: Bin midpoints, velocity values
	_dInvBinWidthVDF = 1./_dBinWidthVDF;
	_dBinMidpointsVDF.resize(_numBinsVDF);
	std::fill (_dBinMidpointsVDF.begin(),_dBinMidpointsVDF.end(),0);

	// local
	_VDF_pjy_abs_local = new uint64_t[_numValsVDF];
	_VDF_pjy_pvx_local = new uint64_t[_numValsVDF];
	_VDF_pjy_pvy_local = new uint64_t[_numValsVDF];
	_VDF_pjy_pvz_local = new uint64_t[_numValsVDF];
	_VDF_pjy_nvx_local = new uint64_t[_numValsVDF];
	_VDF_pjy_nvz_local = new uint64_t[_numValsVDF];

	_VDF_njy_abs_local = new uint64_t[_numValsVDF];
	_VDF_njy_pvx_local = new uint64_t[_numValsVDF];
	_VDF_njy_pvz_local = new uint64_t[_numValsVDF];
	_VDF_njy_nvx_local = new uint64_t[_numValsVDF];
	_VDF_njy_nvy_local = new uint64_t[_numValsVDF];
	_VDF_njy_nvz_local = new uint64_t[_numValsVDF];

	// global
	_VDF_pjy_abs_global = new uint64_t[_numValsVDF];
	_VDF_pjy_pvx_global = new uint64_t[_numValsVDF];
	_VDF_pjy_pvy_global = new uint64_t[_numValsVDF];
	_VDF_pjy_pvz_global = new uint64_t[_numValsVDF];
	_VDF_pjy_nvx_global = new uint64_t[_numValsVDF];
	_VDF_pjy_nvz_global = new uint64_t[_numValsVDF];

	_VDF_njy_abs_global = new uint64_t[_numValsVDF];
	_VDF_njy_pvx_global = new uint64_t[_numValsVDF];
	_VDF_njy_pvz_global = new uint64_t[_numValsVDF];
	_VDF_njy_nvx_global = new uint64_t[_numValsVDF];
	_VDF_njy_nvy_global = new uint64_t[_numValsVDF];
	_VDF_njy_nvz_global = new uint64_t[_numValsVDF];

	// store pointers of local data structures (velocity component only) in 2D array
	_dataPtrs.at(0).at(0) = _VDF_njy_nvx_local;
	_dataPtrs.at(0).at(1) = _VDF_njy_pvx_local;
	_dataPtrs.at(0).at(2) = _VDF_pjy_nvx_local;
	_dataPtrs.at(0).at(3) = _VDF_pjy_pvx_local;
	_dataPtrs.at(1).at(0) = _VDF_njy_nvy_local;
	_dataPtrs.at(1).at(1) = NULL;
	_dataPtrs.at(1).at(2) = NULL;
	_dataPtrs.at(1).at(3) = _VDF_pjy_pvy_local;
	_dataPtrs.at(2).at(0) = _VDF_njy_nvz_local;
	_dataPtrs.at(2).at(1) = _VDF_njy_pvz_local;
	_dataPtrs.at(2).at(2) = _VDF_pjy_nvz_local;
	_dataPtrs.at(2).at(3) = _VDF_pjy_pvz_local;


	// init local values
	this->ResetLocalValuesVDF();
/*
	cout << "_initSamplingVDF = " << _initSamplingVDF << endl;
	cout << "_writeFrequencyVDF = " << _writeFrequencyVDF << endl;
	cout << "_numBinsVDF = " << _numBinsVDF << endl;
	cout << "_numVelocityClassesVDF = " << _numVelocityClassesVDF << endl;
*/

	// discrete velocity values
	this->DoDiscretisationVDF(RS_DIMENSION_Y);

#ifndef NDEBUG
	// show component specific parameters
	this->ShowComponentSpecificParamsVDF();
#endif
}


void SampleRegion::InitSamplingFieldYR(int nDimension)
{
	if(false == _SamplingEnabledFieldYR)
		return;

	// equal shell volumes
	double ra = _dShellWidthInitFieldYR;
	double V = (ra*ra) * M_PI * _dBinWidthFieldYR;
	_dShellVolumeFieldYR = V;
	_dInvShellVolumeFieldYR = 2. / V;  // 2 because of upper/lower section
	double dWidth = (this->GetWidth(0) < this->GetWidth(2) ) ? this->GetWidth(0) : this->GetWidth(2);
	dWidth *= 0.5;
	_nNumShellsFieldYR = ceil( dWidth*dWidth/(ra*ra) );
	_dShellWidthFieldYR = _dShellWidthInitFieldYR;
	_dShellWidthSquaredFieldYR = _dShellWidthFieldYR*_dShellWidthFieldYR;

	// shell volumes
	_dShellVolumesFieldYR = new double[_nNumShellsFieldYR];

	/*
	 * TODO: Implement both variants
	 *
	// equidistant shells
	for(uint32_t si=0; si<_nNumShellsFieldYR; ++si)
	{
		double ri = si * _dShellWidthFieldYR;
		double ra = ri + _dShellWidthFieldYR;
		double V = (ra*ra - ri*ri) * M_PI * _dBinWidthFieldYR;
		_dShellVolumesFieldYR[si] = V;
	}
	*/

	// discrete values: Bin midpoints, Shell midpoints
	_dBinMidpointsFieldYR = new double[_nNumBinsFieldYR];
	_dShellMidpointsFieldYR = new double[_nNumShellsFieldYR];

	_nNumValsFieldYR = _numComponents * 3 * _nNumShellsFieldYR * _nNumBinsFieldYR;  // *3: number of sections: 0: all, 1: upper, 2: lower section

#ifndef NDEBUG
	cout << "_numComponents    = " << _numComponents    << endl;
	cout << "_nNumShellsFieldYR = " << _nNumShellsFieldYR << endl;
	cout << "_nNumBinsFieldYR   = " << _nNumBinsFieldYR   << endl;
	cout << "_nNumValsFieldYR   = " << _nNumValsFieldYR   << endl;
#endif

	// Offsets
	_nOffsetFieldYR = new uint64_t***[3];
	for(uint8_t dim=0; dim<3; ++dim)
	{
		_nOffsetFieldYR[dim] = new uint64_t**[3];
		for(uint8_t sec=0; sec<3; ++sec)
		{
			_nOffsetFieldYR[dim][sec] = new uint64_t*[_numComponents];
			for(uint32_t cid = 0; cid<_numComponents; ++cid)
				_nOffsetFieldYR[dim][sec][cid] = new uint64_t[_nNumShellsFieldYR];
		}
	}

	// Init offsets
	{
		uint64_t nOffset = 0;
		for(uint8_t dim = 0; dim<3; ++dim){
			for(uint8_t sec = 0; sec<3; ++sec){
				for(uint32_t cid = 0; cid<_numComponents; ++cid){
					for(uint32_t si = 0; si<_nNumShellsFieldYR; ++si){
						_nOffsetFieldYR[dim][sec][cid][si] = nOffset;
						nOffset += _nNumBinsFieldYR;
					}
				}
			}
		}
	}

	// Scalar quantities
	// [direction all|+|-][component][position]
	_nNumMoleculesFieldYRLocal  = new unsigned long[_nNumValsFieldYR];
	_nNumMoleculesFieldYRGlobal = new unsigned long[_nNumValsFieldYR];

	// output profiles
	_dDensityFieldYR = new double[_nNumValsFieldYR];
	_dInvShellVolumesFieldYR = new double[_nNumValsFieldYR];

	// section shell volume factor
	 float faSecFac[3] = {1., 2., 2.};

	for(uint32_t cid=0; cid<_numComponents; ++cid){
		for(uint32_t si = 0; si<_nNumShellsFieldYR; ++si)
		{
			double dShellVolume = _dShellVolumesFieldYR[si];
			double dInvShellVolume = 1. / dShellVolume;

			for(uint8_t sec=0; sec<3; ++sec)
			{
				float fSecFac = faSecFac[sec];
				uint64_t nOffset = _nOffsetFieldYR[0][sec][cid][si];
				for(uint32_t bi = 0; bi<_nNumBinsFieldYR; ++bi)
				{
					mardyn_assert( (nOffset+bi) < _nNumValsFieldYR );
					_dInvShellVolumesFieldYR[nOffset+bi] = dInvShellVolume * fSecFac;
				}
			}
		}
	}

	// init sampling data structures
	this->ResetLocalValuesFieldYR();

	// discretisation
	this->DoDiscretisationFieldYR(RS_DIMENSION_Y);
}


void SampleRegion::DoDiscretisationProfiles(int nDimension)
{
	if(false == _SamplingEnabledProfiles)
		return;

	if(_bDiscretisationDoneProfiles == true)  // if allready done -> return
		return;

	double* dLowerCorner = this->GetLowerCorner();

	// calc Bin midpoints
	for(unsigned int s = 0; s < _nNumBinsProfiles; s++)
	{
		_dBinMidpointsProfiles[s] = (s + 0.5) * _dBinWidthProfiles + dLowerCorner[nDimension];
	}

	_bDiscretisationDoneProfiles = true;
}


void SampleRegion::DoDiscretisationVDF(int nDimension)
{
	if(false == _SamplingEnabledVDF)
		return;

	if(_bDiscretisationDoneVDF == true)  // if allready done -> return
		return;

	// calc discrete velocity values
	for(auto&& csp : _vecComponentSpecificParamsVDF)
	{
		for(unsigned int vi=0; vi< csp.numVelocityClasses; ++vi)
			csp.dDiscreteVelocityValues.at(vi) = csp.dVelocityClassWidth * (vi + 0.5);
	}

	// calc Bin midpoints
	double* dLowerCorner = this->GetLowerCorner();
	for(uint32_t bi=0; bi<_numBinsVDF; ++bi)
		_dBinMidpointsVDF.at(bi) = (bi + 0.5) * _dBinWidthVDF + dLowerCorner[nDimension];

	_bDiscretisationDoneVDF = true;
}

void SampleRegion::DoDiscretisationFieldYR(int nDimension)
{
	if(false == _SamplingEnabledFieldYR)
		return;

	if(_bDiscretisationDoneFieldYR == true)  // if allready done -> return
		return;

	double* dLowerCorner = this->GetLowerCorner();

	// calc Bin midpoints
	for(unsigned int bi = 0; bi < _nNumBinsFieldYR; bi++)
	{
		_dBinMidpointsFieldYR[bi] = (bi + 0.5) * _dBinWidthFieldYR + dLowerCorner[nDimension];
	}
	// calc Shell midpoints
	for(unsigned int si = 0; si < _nNumShellsFieldYR; si++)
	{
		_dShellMidpointsFieldYR[si] = (si + 0.5) * _dShellWidthFieldYR + dLowerCorner[nDimension];
	}

	_bDiscretisationDoneFieldYR = true;
}

void SampleRegion::SampleProfiles(Molecule* molecule, int nDimension)
{
	if(false == _SamplingEnabledProfiles)
		return;

	unsigned int nPosIndex;
	unsigned int nIndexMax = _nNumBinsProfiles - 1;

	// do not reset profile matrices here!!!
	// BUT: reset profile before calling this function!!!

	// calc position index
	double* dLowerCorner = this->GetLowerCorner();
	double dPosRelative = molecule->r(nDimension) - dLowerCorner[nDimension];

	nPosIndex = (unsigned int) floor(dPosRelative / _dBinWidthProfiles);

	// ignore outer (halo) molecules
	if(nPosIndex > nIndexMax)  // negative values will be ignored to: cast to unsigned int --> high value
		return;

	unsigned int cid = molecule->componentid() + 1;  // id starts internally with 0
	unsigned int nRotDOF = molecule->component()->getRotationalDegreesOfFreedom();
	double d2EkinTrans = molecule->U2_trans();
	double d2EkinRot   = molecule->U2_rot();
	double v[3];
	v[0] = molecule->v(0);
	v[1] = molecule->v(1);
	v[2] = molecule->v(2);
	double F[3];
	F[0] = molecule->F(0);
	F[1] = molecule->F(1);
	F[2] = molecule->F(2);
	double v2[3];
	v2[0] = v[0]*v[0];
	v2[1] = v[1]*v[1];
	v2[2] = v[2]*v[2];

	// Loop over directions: all (+/-) | only (+) | only (-)
	for(unsigned int dir = 0; dir<3; ++dir)
	{
		// only (+)
		if(1==dir && v[1] < 0.)
			continue;
		// only (-)
		if(2==dir && v[1] > 0.)
			continue;

		mardyn_assert(_nOffsetScalar[dir][0  ] + nPosIndex < _nNumValsScalar);
		mardyn_assert(_nOffsetScalar[dir][cid] + nPosIndex < _nNumValsScalar);

		// Scalar quantities
		_nNumMoleculesLocal[ _nOffsetScalar[dir][0  ] + nPosIndex ] ++;  // all components
		_nNumMoleculesLocal[ _nOffsetScalar[dir][cid] + nPosIndex ] ++;  // specific component
		_nRotDOFLocal      [ _nOffsetScalar[dir][0  ] + nPosIndex ] += nRotDOF;
		_nRotDOFLocal      [ _nOffsetScalar[dir][cid] + nPosIndex ] += nRotDOF;

		_d2EkinRotLocal    [ _nOffsetScalar[dir][0  ] + nPosIndex ] += d2EkinRot;
		_d2EkinRotLocal    [ _nOffsetScalar[dir][cid] + nPosIndex ] += d2EkinRot;

		// Vector quantities
		// Loop over dimensions  x, y, z (vector components)
		for(unsigned int dim = 0; dim<3; ++dim)
		{
			mardyn_assert(_nOffsetVector[dim][dir][0  ] + nPosIndex < _nNumValsVector);
			mardyn_assert(_nOffsetVector[dim][dir][cid] + nPosIndex < _nNumValsVector);

			_dVelocityLocal       [ _nOffsetVector[dim][dir][0  ] + nPosIndex ] += v[dim];
			_dVelocityLocal       [ _nOffsetVector[dim][dir][cid] + nPosIndex ] += v[dim];
			_dSquaredVelocityLocal[ _nOffsetVector[dim][dir][0  ] + nPosIndex ] += v2[dim];
			_dSquaredVelocityLocal[ _nOffsetVector[dim][dir][cid] + nPosIndex ] += v2[dim];
			_dForceLocal          [ _nOffsetVector[dim][dir][0  ] + nPosIndex ] += F[dim];
			_dForceLocal          [ _nOffsetVector[dim][dir][cid] + nPosIndex ] += F[dim];
		}
	}
}


void SampleRegion::SampleVDF(Molecule* molecule, int nDimension)
{

	if(false == _SamplingEnabledVDF)
		return;

	// return if discretisation is not done yet
	// reason: v_max has to be determined first (when not set manually)
	if(false == _bDiscretisationDoneVDF)
		return;

	uint32_t cid = molecule->componentid()+1;  // 0: all components
	const ComponentSpecificParamsVDF& csp = _vecComponentSpecificParamsVDF.at(cid);
	if(false == csp.bSamplingEnabled)
		return;
	uint32_t nComponentOffset = csp.nOffsetDataStructure;

	// calc bin index / offset
	double* dLowerCorner = this->GetLowerCorner();
	double dPosRelative = molecule->r(nDimension) - dLowerCorner[nDimension];
	uint32_t nBinIndex = (uint32_t) floor(dPosRelative * _dInvBinWidthVDF);
	uint32_t numVelocityClasses   = csp.numVelocityClasses;
	uint32_t nBinOffset = nBinIndex * numVelocityClasses;
	uint32_t nOffset = nComponentOffset + nBinOffset;

	// ignore outer (halo) molecules
	uint32_t nIndexMax = _numBinsVDF - 1;
	if(nBinIndex > nIndexMax)
		return;

	double absVal;
	_f2ptr(absVal, molecule);
	absVal = sqrt(absVal);
	double dInvVelocityClassWidth = csp.dInvVelocityClassWidth;
	uint32_t nVelocityClassIndex = (uint32_t)(floor(absVal * dInvVelocityClassWidth) );

	// calculate velocity vector indices for velocity components
	double v[3];
	uint32_t naVelocityClassIndex[3];

	_fptr(v, molecule);  // sample either velocity or force vector
	for(unsigned int d=0; d<3; ++d) {
		//v[d] = molecule->v(d);
		naVelocityClassIndex[d] = (uint32_t)(floor( fabs( v[d] ) * dInvVelocityClassWidth) );
	}

	// respect finite resolution of velocity
	uint32_t nIndexMaxVelo = numVelocityClasses - 1;
	for(unsigned int d=0; d<3; ++d) {
		if(naVelocityClassIndex[d] > nIndexMaxVelo)
			return;
	}

	// sampling
	bool bSampleComponentSum = _vecComponentSpecificParamsVDF.at(0).bSamplingEnabled;
	// velocity components
	for(unsigned int d=0; d<3; ++d)
	{
		uint8_t ptrIndex = 2*(v[1]>0.)+(v[d]>0.);
		_dataPtrs.at(d).at(ptrIndex)[ nOffset + naVelocityClassIndex[d] ]++;
		if(true == bSampleComponentSum)
			_dataPtrs.at(d).at(ptrIndex)[ nBinOffset + naVelocityClassIndex[d] ]++;
	}

	// absolute velocity
	if(nVelocityClassIndex > nIndexMaxVelo)  // v_abs > v_max? (respect finite resolution of velocity)
		return;

	if(v[1] > 0.)  // particle flux in positive y-direction
	{
		if(true == bSampleComponentSum)
			_VDF_pjy_abs_local[nBinOffset + nVelocityClassIndex]++;
		_VDF_pjy_abs_local[nOffset + nVelocityClassIndex]++;
	}
	else  // particle flux in negative y-direction
	{
		if(true == bSampleComponentSum)
			_VDF_njy_abs_local[nBinOffset + nVelocityClassIndex]++;
		_VDF_njy_abs_local[nOffset + nVelocityClassIndex]++;
	}
}

void SampleRegion::SampleFieldYR(Molecule* molecule)
{
	if(false == _SamplingEnabledFieldYR)
		return;

	uint32_t nPosIndexY;
	uint32_t nIndexMaxY = _nNumBinsFieldYR - 1;
	uint32_t nPosIndexR;
	uint32_t nIndexMaxR = _nNumShellsFieldYR - 1;

	// do not reset profile matrices here!!!
	// BUT: reset profile before calling this function!!!

	// calc position index
	double* dLowerCorner = this->GetLowerCorner();
	double dPosRelativeX = molecule->r(0) - (dLowerCorner[0] + this->GetWidth(0)*0.5);
	double dPosRelativeY = molecule->r(1) -  dLowerCorner[1];
	double dPosRelativeZ = molecule->r(2) - (dLowerCorner[2] + this->GetWidth(2)*0.5);

	nPosIndexY = (unsigned int) floor(dPosRelativeY / _dBinWidthFieldYR);
	double dRadiusSquared = (dPosRelativeX*dPosRelativeX + dPosRelativeZ*dPosRelativeZ);
//	double dRadius = sqrt(dPosRelativeX*dPosRelativeX + dPosRelativeZ*dPosRelativeZ);
	nPosIndexR = (unsigned int) floor(dRadiusSquared / _dShellWidthSquaredFieldYR);

	// ignore outer (halo) molecules
	if(nPosIndexY > nIndexMaxY)  // negative values will be ignored too: cast to unsigned int --> high value
		return;

	// ignore outer (halo) molecules
	if(nPosIndexR > nIndexMaxR)  // negative values will be ignored too: cast to unsigned int --> high value
		return;

	unsigned int cid = molecule->componentid() + 1;  // id starts internally with 0

	for(uint8_t sec=0; sec<3; ++sec)
	{
		mardyn_assert(_nOffsetFieldYR[0][sec][0  ][nPosIndexR] + nPosIndexY < _nNumValsFieldYR);
		mardyn_assert(_nOffsetFieldYR[0][sec][cid][nPosIndexR] + nPosIndexY < _nNumValsFieldYR);
	}

	// Scalar quantities
	_nNumMoleculesFieldYRLocal[ _nOffsetFieldYR[0][0][0  ][nPosIndexR] + nPosIndexY ] ++;  // all components
	_nNumMoleculesFieldYRLocal[ _nOffsetFieldYR[0][0][cid][nPosIndexR] + nPosIndexY ] ++;  // specific component

	if(dPosRelativeX >= 0.)
	{
		_nNumMoleculesFieldYRLocal[ _nOffsetFieldYR[0][1][0  ][nPosIndexR] + nPosIndexY ] ++;  // all components
		_nNumMoleculesFieldYRLocal[ _nOffsetFieldYR[0][1][cid][nPosIndexR] + nPosIndexY ] ++;  // specific component
	}
	else
	{
		_nNumMoleculesFieldYRLocal[ _nOffsetFieldYR[0][2][0  ][nPosIndexR] + nPosIndexY ] ++;  // all components
		_nNumMoleculesFieldYRLocal[ _nOffsetFieldYR[0][2][cid][nPosIndexR] + nPosIndexY ] ++;  // specific component
	}
}

void SampleRegion::CalcGlobalValuesProfiles(DomainDecompBase* domainDecomp, Domain* domain)
{
	if(false == _SamplingEnabledProfiles)
		return;

	// perform reduce operation, process further calculations
#ifdef ENABLE_MPI

	// Scalar quantities
	// [direction all|+|-][component][position]
	MPI_Reduce( _nNumMoleculesLocal, _nNumMoleculesGlobal, _nNumValsScalar, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce( _nRotDOFLocal,       _nRotDOFGlobal,       _nNumValsScalar, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce( _d2EkinRotLocal,     _d2EkinRotGlobal,     _nNumValsScalar, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	// Vector quantities
	// [dimension x|y|z][direction all|+|-][component][position]
	MPI_Reduce( _dVelocityLocal,        _dVelocityGlobal,        _nNumValsVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce( _dSquaredVelocityLocal, _dSquaredVelocityGlobal, _nNumValsVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce( _dForceLocal,           _dForceGlobal,           _nNumValsVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

#else
	// Scalar quantities
	for(unsigned int i = 0; i < _nNumValsScalar; ++i)
	{
		_nNumMoleculesGlobal[i] = _nNumMoleculesLocal[i];
		_nRotDOFGlobal[i]       = _nRotDOFLocal[i];
		_d2EkinRotGlobal[i]     = _d2EkinRotLocal[i];
	}

	// Vector quantities
	for(unsigned int i = 0; i < _nNumValsVector; ++i)
	{
		_dVelocityGlobal[i]        = _dVelocityLocal[i];
		_dSquaredVelocityGlobal[i] = _dSquaredVelocityLocal[i];
		_dForceGlobal[i]           = _dForceLocal[i];
	}
#endif

	int rank = domainDecomp->getRank();
	//  int numprocs = domainDecomp->getNumProcs();

	// only root process writes out data
	if(rank != 0)
		 return;

	// reset data structures of kinetic energy, before cumsum is for all components is calculated
	this->ResetOutputDataProfiles();

	double dInvertBinVolume = 1. / _dBinVolumeProfiles;

	double dInvertNumSamples = (double) (_writeFrequencyProfiles);  // TODO: perhabs in future different from writeFrequency
	dInvertNumSamples = 1. / dInvertNumSamples;

	// sum for cid == 0 (all components)
	for(uint8_t dir=0; dir<3; ++dir)
	{
		for(uint8_t dim=0; dim<3; ++dim)
		{
			for(uint8_t cid=1; cid<_numComponents; ++cid)
			{
				unsigned long offsetN    =      _nOffsetScalar[dir][cid];
				unsigned long offsetComp = _nOffsetVector[dim][dir][cid];
				unsigned long offsetSum  = _nOffsetVector[dim][dir][  0];

				for(uint32_t s = 0; s < _nNumBinsProfiles; ++s)
				{
					mardyn_assert( (offsetSum+s)  < _nNumValsVector );
					mardyn_assert( (offsetComp+s) < _nNumValsVector );

					// Ekin trans.
					double v2 = _dSquaredVelocityGlobal[offsetComp+s];
					double d2EkinTrans = v2 * _vecMass.at(cid);
					_d2EkinTransComp[offsetComp+s] = d2EkinTrans;
					_d2EkinTransComp[offsetSum+s] += d2EkinTrans;
					// Ekin drift
					unsigned long N = _nNumMoleculesGlobal[offsetN+s];
					double dInvN = (N > 0) ? (1. / ( (double)(N) ) ) : (0.);
					double vd = _dVelocityGlobal[offsetComp+s];
					double d2EkinDrift = vd*vd * _vecMass.at(cid) * dInvN;
					_d2EkinDriftComp[offsetComp+s] = d2EkinDrift;
					_d2EkinDriftComp[offsetSum+s] += d2EkinDrift;
				}
			}
		}
	}

	// sum x, y, z --> scalar
	for(uint8_t dir=0; dir<3; ++dir)
	{
		for(uint8_t dim=0; dim<3; ++dim)
		{
			for(uint8_t cid=0; cid<_numComponents; ++cid)
			{
				unsigned long offsetDim = _nOffsetVector[dim][dir][cid];
				unsigned long offsetSum = _nOffsetVector[  0][dir][cid];  // sum x, y, z --> scalar

				for(uint32_t s = 0; s < _nNumBinsProfiles; ++s)
				{
					mardyn_assert( (offsetDim+s) < _nNumValsVector );
					mardyn_assert( (offsetSum+s) < _nNumValsScalar );

					_d2EkinTrans[offsetSum+s] += _d2EkinTransComp[offsetDim+s];
					_d2EkinDrift[offsetSum+s] += _d2EkinDriftComp[offsetDim+s];
				}
			}
		}
	}

	// Scalar quantities
	for(unsigned int i = 0; i < _nNumValsScalar; ++i)
	{
		unsigned long nNumMolecules = _nNumMoleculesGlobal[i];
		unsigned long nDOF_trans = nNumMolecules * 3;
		unsigned long nDOF_rot   = _nRotDOFGlobal[i];
		unsigned long nDOF_total = nDOF_trans + nDOF_rot;
		double dInvDOF_trans = (nDOF_trans > 0) ? (1. / ( (double)(nDOF_trans) ) ) : (0.);
		double dInvDOF_rot   = (nDOF_rot   > 0) ? (1. / ( (double)(nDOF_rot  ) ) ) : (0.);
		double dInvDOF_total = (nDOF_total > 0) ? (1. / ( (double)(nDOF_total) ) ) : (0.);

		// density profile
		_dDensity [i] = nNumMolecules * _dInvertBinVolSamplesProfiles;

		double d2EkinTrans = _d2EkinTrans[i];
		double d2EkinDrift = _d2EkinDrift[i];
		double d2EkinRot   = _d2EkinRotGlobal[i];
		double d2EkinTotal = d2EkinTrans + d2EkinRot;
		double d2EkinT     = d2EkinTotal - d2EkinDrift;
		_d2EkinRot[i]   = d2EkinRot;
		_d2EkinTotal[i] = d2EkinTotal;
		_d2EkinT[i]     = d2EkinT;
		_dTemperature[i]      = d2EkinT     * dInvDOF_total;
		_dTemperatureTrans[i] = d2EkinTrans * dInvDOF_trans;
		_dTemperatureRot[i]   = d2EkinRot   * dInvDOF_rot;
	}

	// Vector quantities
	for(unsigned int dim = 0; dim<3; ++dim)
	{
		unsigned long nDimOffset = _nNumValsScalar*dim;

		for(unsigned int i = 0; i < _nNumValsScalar; ++i)
		{
			unsigned long nNumMolecules = _nNumMoleculesGlobal[i];
			double dInvertNumMolecules = 1. / ( (double)(nNumMolecules) );

			_dDriftVelocity[nDimOffset+i] = _dVelocityGlobal[nDimOffset+i] * dInvertNumMolecules;
			_dForce        [nDimOffset+i] = _dForceGlobal   [nDimOffset+i] * dInvertNumMolecules;
			double d2EkinTrans = _d2EkinTransComp[nDimOffset+i];
			double d2EkinDrift = _d2EkinDriftComp[nDimOffset+i];
//			cout << "nDimOffset+i = " << nDimOffset+i << endl;
//			cout << "dEkinTrans = " << dEkinTrans << endl;
//			cout << "dEkinDrift = " << dEkinDrift << endl;
			_dTemperatureComp[nDimOffset+i] = (d2EkinTrans - d2EkinDrift) * dInvertNumMolecules;
		}
	}
}


void SampleRegion::CalcGlobalValuesVDF()
{
	if(false == _SamplingEnabledVDF)
		return;

	#ifdef ENABLE_MPI

		// positive y-direction
		MPI_Reduce( _VDF_pjy_abs_local, _VDF_pjy_abs_global, _numValsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

		MPI_Reduce( _VDF_pjy_pvx_local, _VDF_pjy_pvx_global, _numValsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce( _VDF_pjy_pvy_local, _VDF_pjy_pvy_global, _numValsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce( _VDF_pjy_pvz_local, _VDF_pjy_pvz_global, _numValsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

		MPI_Reduce( _VDF_pjy_nvx_local, _VDF_pjy_nvx_global, _numValsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce( _VDF_pjy_nvz_local, _VDF_pjy_nvz_global, _numValsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);


		// negative y-direction
		MPI_Reduce( _VDF_njy_abs_local, _VDF_njy_abs_global, _numValsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

		MPI_Reduce( _VDF_njy_pvx_local, _VDF_njy_pvx_global, _numValsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce( _VDF_njy_pvz_local, _VDF_njy_pvz_global, _numValsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

		MPI_Reduce( _VDF_njy_nvx_local, _VDF_njy_nvx_global, _numValsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce( _VDF_njy_nvy_local, _VDF_njy_nvy_global, _numValsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce( _VDF_njy_nvz_local, _VDF_njy_nvz_global, _numValsVDF, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

	#else
		for(uint32_t vi=0; vi<_numValsVDF; ++vi)
		{
			_VDF_pjy_abs_global[vi] = _VDF_pjy_abs_local[vi];
			_VDF_pjy_pvx_global[vi] = _VDF_pjy_pvx_local[vi];
			_VDF_pjy_pvy_global[vi] = _VDF_pjy_pvy_local[vi];
			_VDF_pjy_pvz_global[vi] = _VDF_pjy_pvz_local[vi];
			_VDF_pjy_nvx_global[vi] = _VDF_pjy_nvx_local[vi];
			_VDF_pjy_nvz_global[vi] = _VDF_pjy_nvz_local[vi];

			_VDF_njy_abs_global[vi] = _VDF_njy_abs_local[vi];
			_VDF_njy_pvx_global[vi] = _VDF_njy_pvx_local[vi];
			_VDF_njy_pvz_global[vi] = _VDF_njy_pvz_local[vi];
			_VDF_njy_nvx_global[vi] = _VDF_njy_nvx_local[vi];
			_VDF_njy_nvy_global[vi] = _VDF_njy_nvy_local[vi];
			_VDF_njy_nvz_global[vi] = _VDF_njy_nvz_local[vi];
		}
	#endif
}

void SampleRegion::CalcGlobalValuesFieldYR(DomainDecompBase* domainDecomp, Domain* domain)
{
	if(false == _SamplingEnabledFieldYR)
		return;

	// perform reduce operation, process further calculations
#ifdef ENABLE_MPI

	// Scalar quantities
	// [dimension x|y|z][component][positionR][positionY]
	MPI_Reduce( _nNumMoleculesFieldYRLocal, _nNumMoleculesFieldYRGlobal, _nNumValsFieldYR, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

#else
	// Scalar quantities
	for(unsigned int i = 0; i < _nNumValsFieldYR; ++i)
	{
		_nNumMoleculesFieldYRGlobal[i] = _nNumMoleculesFieldYRLocal[i];
	}
#endif

	int rank = domainDecomp->getRank();
	//  int numprocs = domainDecomp->getNumProcs();

	// only root process writes out data
	if(rank != 0)
		 return;

	double dInvertNumSamples = (double) (_writeFrequencyFieldYR);  // TODO: perhabs in future different from writeFrequency
	dInvertNumSamples = 1. / dInvertNumSamples;

	// Scalar quantities
	for(unsigned int i=0; i<_nNumValsFieldYR; ++i)
	{
		uint64_t nNumMolecules = _nNumMoleculesFieldYRGlobal[i];

		// density profile
//		_dDensityFieldYR[i] = nNumMolecules * _dInvShellVolumesFieldYR[i] * dInvertNumSamples;  // equidistant
		_dDensityFieldYR[i] = nNumMolecules * _dInvShellVolumeFieldYR * dInvertNumSamples;
	}
}


void SampleRegion::WriteDataProfiles(DomainDecompBase* domainDecomp, unsigned long simstep, Domain* domain)
{
	if(false == _SamplingEnabledProfiles)
		return;

	// sampling starts after initial timestep (_initSamplingVDF) and with respect to write frequency (_writeFrequencyVDF)
	if( simstep <= _initSamplingProfiles )
		return;

	if ( (simstep - _initSamplingProfiles) % _writeFrequencyProfiles != 0 )
		return;

	// calc global values
	this->CalcGlobalValuesProfiles(domainDecomp, domain);

	// reset local values
	this->ResetLocalValuesProfiles();

	// writing .dat-files
	std::stringstream filenamestream_scal[3];
	std::stringstream filenamestream_vect[3];
	std::string dir_prefix[3] = {"all", "pos", "neg"};
	for(uint8_t dir=0; dir<3; dir++) {
		filenamestream_scal[dir] << "scalquant_" << dir_prefix[dir] << "_reg" << this->GetID() << "_TS" << fill_width('0', 9) << simstep << ".dat";
		filenamestream_vect[dir] << "vectquant_" << dir_prefix[dir] << "_reg" << this->GetID() << "_TS" << fill_width('0', 9) << simstep << ".dat";
	}

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	// int numprocs = domainDecomp->getNumProcs();
	if (rank!= 0)
		return;
#endif

	for(uint8_t dir=0; dir<3; ++dir)
	{
		std::stringstream outputstream_scal;
		std::stringstream outputstream_vect;

		// header
		outputstream_scal << "                     pos";
		outputstream_vect << "                     pos";
		for(uint32_t cid = 0; cid < _numComponents; ++cid)
		{
			// scalar
			outputstream_scal << "            DOF_total[" << cid << "]";
			outputstream_scal << "            DOF_trans[" << cid << "]";
			outputstream_scal << "              DOF_rot[" << cid << "]";
			outputstream_scal << "                  rho[" << cid << "]";
			outputstream_scal << "           Ekin_total[" << cid << "]";
			outputstream_scal << "           Ekin_trans[" << cid << "]";
			outputstream_scal << "           Ekin_drift[" << cid << "]";
			outputstream_scal << "             Ekin_rot[" << cid << "]";
			outputstream_scal << "               Ekin_T[" << cid << "]";
			outputstream_scal << "                 Epot[" << cid << "]";
			outputstream_scal << "                    T[" << cid << "]";
			outputstream_scal << "              T_trans[" << cid << "]";
			outputstream_scal << "                T_rot[" << cid << "]";

			// vector
			outputstream_vect << "                   Fx[" << cid << "]";
			outputstream_vect << "                   Fy[" << cid << "]";
			outputstream_vect << "                   Fz[" << cid << "]";
			outputstream_vect << "                   vx[" << cid << "]";
			outputstream_vect << "                   vy[" << cid << "]";
			outputstream_vect << "                   vz[" << cid << "]";
			outputstream_vect << "         Ekin_trans,x[" << cid << "]";
			outputstream_vect << "         Ekin_trans,y[" << cid << "]";
			outputstream_vect << "         Ekin_trans,z[" << cid << "]";
			outputstream_vect << "         Ekin_drift,x[" << cid << "]";
			outputstream_vect << "         Ekin_drift,y[" << cid << "]";
			outputstream_vect << "         Ekin_drift,z[" << cid << "]";
			outputstream_vect << "                   Tx[" << cid << "]";
			outputstream_vect << "                   Ty[" << cid << "]";
			outputstream_vect << "                   Tz[" << cid << "]";
		}
		outputstream_scal << endl;
		outputstream_vect << endl;

		// data
		for(uint32_t s=0; s<_nNumBinsProfiles; ++s)
		{
			outputstream_scal << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dBinMidpointsProfiles[s];
			outputstream_vect << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dBinMidpointsProfiles[s];

			for(uint8_t cid=0; cid<_numComponents; ++cid)
			{
				// scalar
				mardyn_assert( _nOffsetScalar[dir][cid]+s < _nNumValsScalar );

				unsigned long offset = _nOffsetScalar[dir][cid]+s;
				unsigned long DOF_trans = _nNumMoleculesGlobal[offset] * 3;
				unsigned long DOF_rot   = _nRotDOFGlobal[offset];
				unsigned long DOF_total = DOF_trans + DOF_rot;
				double rho = _dDensity[offset];
				double Ekin_total = _d2EkinTotal[offset] * 0.5;
				double Ekin_trans = _d2EkinTrans[offset] * 0.5;
				double Ekin_drift = _d2EkinDrift[offset] * 0.5;
				double Ekin_rot   = _d2EkinRot  [offset] * 0.5;
				double Ekin_T     = _d2EkinT    [offset] * 0.5;
				double Epot = 0.;
				double T = _dTemperature[offset];
				double T_trans = _dTemperatureTrans[offset];
				double T_rot   = _dTemperatureRot[offset];
				outputstream_scal << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << DOF_total;
				outputstream_scal << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << DOF_trans;
				outputstream_scal << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << DOF_rot;
				outputstream_scal << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << rho;
				outputstream_scal << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << Ekin_total * _dInvertNumSamplesProfiles;;
				outputstream_scal << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << Ekin_trans * _dInvertNumSamplesProfiles;;
				outputstream_scal << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << Ekin_drift * _dInvertNumSamplesProfiles;;
				outputstream_scal << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << Ekin_rot   * _dInvertNumSamplesProfiles;;
				outputstream_scal << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << Ekin_T     * _dInvertNumSamplesProfiles;;
				outputstream_scal << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << Epot       * _dInvertNumSamplesProfiles;;
				outputstream_scal << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << T;
				outputstream_scal << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << T_trans;
				outputstream_scal << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << T_rot;

				// vector
				mardyn_assert( _nOffsetVector[0][dir][cid]+s < _nNumValsVector );
				mardyn_assert( _nOffsetVector[1][dir][cid]+s < _nNumValsVector );
				mardyn_assert( _nOffsetVector[2][dir][cid]+s < _nNumValsVector );

				unsigned long offset_x = _nOffsetVector[0][dir][cid]+s;
				unsigned long offset_y = _nOffsetVector[1][dir][cid]+s;
				unsigned long offset_z = _nOffsetVector[2][dir][cid]+s;

				outputstream_vect << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dForce[offset_x];
				outputstream_vect << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dForce[offset_y];
				outputstream_vect << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dForce[offset_z];
				outputstream_vect << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dDriftVelocity[offset_x];
				outputstream_vect << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dDriftVelocity[offset_y];
				outputstream_vect << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dDriftVelocity[offset_z];
				outputstream_vect << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _d2EkinTransComp[offset_x] * _dInvertNumSamplesProfiles * 0.5;
				outputstream_vect << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _d2EkinTransComp[offset_y] * _dInvertNumSamplesProfiles * 0.5;
				outputstream_vect << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _d2EkinTransComp[offset_z] * _dInvertNumSamplesProfiles * 0.5;
				outputstream_vect << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _d2EkinDriftComp[offset_x] * _dInvertNumSamplesProfiles * 0.5;
				outputstream_vect << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _d2EkinDriftComp[offset_y] * _dInvertNumSamplesProfiles * 0.5;
				outputstream_vect << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _d2EkinDriftComp[offset_z] * _dInvertNumSamplesProfiles * 0.5;
				outputstream_vect << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dTemperatureComp[offset_x];
				outputstream_vect << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dTemperatureComp[offset_y];
				outputstream_vect << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << _dTemperatureComp[offset_z];
			} // loop: cid
			outputstream_scal << endl;
			outputstream_vect << endl;
		} // loop: pos

		// open file for writing
		// scalar
		ofstream fileout_scal(filenamestream_scal[dir].str().c_str(), ios::out);
		fileout_scal << outputstream_scal.str();
		fileout_scal.close();
		// vector
		ofstream fileout_vect(filenamestream_vect[dir].str().c_str(), ios::out);
		fileout_vect << outputstream_vect.str();
		fileout_vect.close();
	} // loop: dir
}


void SampleRegion::WriteDataVDF(DomainDecompBase* domainDecomp, unsigned long simstep)
{
	if(false == _SamplingEnabledVDF)
		return;

	// sampling starts after initial timestep (_initSamplingVDF) and with respect to write frequency (_writeFrequencyVDF)
	if( simstep <= _initSamplingVDF )
		return;

	if ( (simstep - _initSamplingVDF) % _writeFrequencyVDF != 0 )
		return;

	// calc global values
	this->CalcGlobalValuesVDF();  // calculate global velocity distribution sums

	// reset local values
	this->ResetLocalValuesVDF();

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	// int numprocs = domainDecomp->getNumProcs();
	if (rank != 0)
		return;
#endif

	const uint8_t numDataStructuresVDF = 12;
	const uint8_t numFiles = numDataStructuresVDF+1;
	const std::array<std::string,numFiles> strSubPrefixes = {
			"_pjy_abs",
			"_pjy_pvx",
			"_pjy_pvy",
			"_pjy_pvz",
			"_pjy_nvx",
			"_pjy_nvz",
			"_njy_abs",
			"_njy_pvx",
			"_njy_pvz",
			"_njy_nvx",
			"_njy_nvy",
			"_njy_nvz",
			"_classes"};

	std::array<uint64_t*,12> dataPtrs;
	dataPtrs.at(0)  = _VDF_pjy_abs_global;
	dataPtrs.at(1)  = _VDF_pjy_pvx_global;
	dataPtrs.at(2)  = _VDF_pjy_pvy_global;
	dataPtrs.at(3)  = _VDF_pjy_pvz_global;
	dataPtrs.at(4)  = _VDF_pjy_nvx_global;
	dataPtrs.at(5)  = _VDF_pjy_nvz_global;
	dataPtrs.at(6)  = _VDF_njy_abs_global;
	dataPtrs.at(7)  = _VDF_njy_pvx_global;
	dataPtrs.at(8)  = _VDF_njy_pvz_global;
	dataPtrs.at(9)  = _VDF_njy_nvx_global;
	dataPtrs.at(10) = _VDF_njy_nvy_global;
	dataPtrs.at(11) = _VDF_njy_nvz_global;

	for(uint64_t cid=0; cid<_numComponents; ++cid)
	{
		const ComponentSpecificParamsVDF& csp = _vecComponentSpecificParamsVDF.at(cid);
		if(false == csp.bSamplingEnabled)
			continue;

		std::stringstream sstrOutput[numFiles];
		std::stringstream sstrFilename[numFiles];

		// concat filename string
		for(uint32_t fi=0; fi<numFiles; ++fi)
			sstrFilename[fi] << _fnamePrefixVDF << "_reg" << this->GetID() << "_cid" << cid << strSubPrefixes.at(fi) << "_TS" << fill_width('0', 9) << simstep << ".dat";

		// write to string streams
		sstrOutput[numFiles-1] << "            classes_cid" << cid << endl;
		for(auto&& dvv : csp.dDiscreteVelocityValues)
		{
			sstrOutput[numFiles-1] << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << dvv;
			sstrOutput[numFiles-1] << endl;
		}

		// loop over vdf data structures
		for(uint8_t di=0; di<numDataStructuresVDF; ++di)
		{
			// meta data
			sstrOutput[di].write(reinterpret_cast<const char*>(&csp.numVelocityClasses), 4);
			sstrOutput[di].write(reinterpret_cast<const char*>(&_numBinsVDF), 4);

			// vdf data
			uint64_t* data = dataPtrs.at(di);
			for(uint64_t vi=csp.nOffsetDataStructure; vi<csp.nOffsetDataStructure+csp.numVals; ++vi)
			{
				mardyn_assert(vi < _numValsVDF);
				sstrOutput[di].write(reinterpret_cast<const char*>(&data[vi]), 8);
//				sstrOutput[di].write(reinterpret_cast<const char*>(&vi), 8);
			}
		}

		// write to file streams
		for(uint8_t fi=0; fi<numFiles; ++fi)
		{
			ofstream ofs(sstrFilename[fi].str().c_str(), ios::out);
			ofs << sstrOutput[fi].str();
			ofs.close();
		}
	}

	// bin midpoint positions
	std::stringstream sstrOutput;
	std::stringstream sstrFilename;
	sstrFilename << "VDF_reg" << this->GetID() << "_bin_coords" << "_TS" << fill_width('0', 9) << simstep << ".dat";

	// write to string stream
	sstrOutput << "                  coords" << endl;
	for(auto&& bmp : _dBinMidpointsVDF)
	{
		sstrOutput << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << bmp;
		sstrOutput << endl;
	}
	// write to file stream
	ofstream ofs(sstrFilename.str().c_str(), ios::out);
	ofs << sstrOutput.str();
	ofs.close();
}


void SampleRegion::WriteDataFieldYR(DomainDecompBase* domainDecomp, unsigned long simstep, Domain* domain)
{
	if(false == _SamplingEnabledFieldYR)
		return;

	// sampling starts after initial timestep (_initSamplingVDF) and with respect to write frequency (_writeFrequencyVDF)
	if( simstep <= _initSamplingFieldYR )
		return;

	if ( (simstep - _initSamplingFieldYR) % _writeFrequencyFieldYR != 0 )
		return;

	// calc global values
	this->CalcGlobalValuesFieldYR(domainDecomp, domain);

	// reset local values
	this->ResetLocalValuesFieldYR();

#ifdef ENABLE_MPI
	int rank = domainDecomp->getRank();
	// int numprocs = domainDecomp->getNumProcs();
	if (rank!= 0)
		return;
#endif

	for(uint8_t sec=0; sec<3; ++sec)
	{
		// writing .dat-files
		std::stringstream filenamestream;
		filenamestream << _strFilePrefixFieldYR << "_sec" << (uint32_t)sec << "_reg" << this->GetID() << "_TS" << fill_width('0', 9) << simstep << ".dat";

		std::stringstream outputstream;
		uint64_t nOffset = 0;

		if(RSFT_ASCII == _nFileTypeFieldYR)
		{
			for(uint32_t si=0; si<_nNumShellsFieldYR; ++si)
			{
				nOffset = _nOffsetFieldYR[0][sec][0][si];
				for(uint32_t bi=0; bi<_nNumBinsFieldYR; ++bi)
				{
					mardyn_assert( (nOffset+bi) < _nNumValsFieldYR );
					outputstream << FORMAT_SCI_MAX_DIGITS << _dDensityFieldYR[nOffset+bi];
				}
				outputstream << endl;
			}
			ofstream fileout(filenamestream.str().c_str(), std::ios::out);
			fileout << outputstream.str();
			fileout.close();
		}
		else
		{
			outputstream.write(reinterpret_cast<const char*>(&_nNumBinsFieldYR), 4);
			outputstream.write(reinterpret_cast<const char*>(&_nNumShellsFieldYR), 4);
			for(uint32_t si=0; si<_nNumShellsFieldYR; ++si)
			{
				nOffset = _nOffsetFieldYR[0][sec][0][si];
				for(uint32_t bi=0; bi<_nNumBinsFieldYR; ++bi)
				{
					mardyn_assert( (nOffset+bi) < _nNumValsFieldYR );
					double dVal = _dDensityFieldYR[nOffset+bi];
					outputstream.write(reinterpret_cast<const char*>(&dVal), 8);
				}
			}
			ofstream fileout(filenamestream.str().c_str(), std::ios::out | std::ios::binary);
			fileout << outputstream.str();
			fileout.close();
		}
	}
}


// private methods

void SampleRegion::ResetLocalValuesVDF()
{
	if(false == _SamplingEnabledVDF)
		return;

	for(uint32_t vi=0; vi<_numValsVDF; ++vi)
	{
		_VDF_pjy_abs_local[vi] = 0;
		_VDF_pjy_pvx_local[vi] = 0;
		_VDF_pjy_pvy_local[vi] = 0;
		_VDF_pjy_pvz_local[vi] = 0;
		_VDF_pjy_nvx_local[vi] = 0;
		_VDF_pjy_nvz_local[vi] = 0;

		_VDF_njy_abs_local[vi] = 0;
		_VDF_njy_pvx_local[vi] = 0;
		_VDF_njy_pvz_local[vi] = 0;
		_VDF_njy_nvx_local[vi] = 0;
		_VDF_njy_nvy_local[vi] = 0;
		_VDF_njy_nvz_local[vi] = 0;
	}
}


void SampleRegion::ResetLocalValuesProfiles()
{
	if(false == _SamplingEnabledProfiles)
		return;

	// Scalar quantities
	for(unsigned int i = 0; i < _nNumValsScalar; ++i)
	{
		_nNumMoleculesLocal[i] = 0;
		_nRotDOFLocal[i] = 0;
		_d2EkinRotLocal[i]   = 0.;
	}

	// Vector quantities
	for(unsigned int i = 0; i < _nNumValsVector; ++i)
	{
		_dVelocityLocal[i] = 0.;
		_dSquaredVelocityLocal[i] = 0.;
		_dForceLocal[i] = 0.;
	}
}

void SampleRegion::ResetOutputDataProfiles()
{
	if(false == _SamplingEnabledProfiles)
		return;

	// Scalar quantities
	for(unsigned int i = 0; i < _nNumValsScalar; ++i)
	{
		_d2EkinTrans[i] = 0.;
		_d2EkinDrift[i] = 0.;
	}

	// Vector quantities
	for(unsigned int i = 0; i < _nNumValsVector; ++i)
	{
		_d2EkinTransComp[i] = 0.;
		_d2EkinDriftComp[i] = 0.;
	}
}

void SampleRegion::ResetLocalValuesFieldYR()
{
	if(false == _SamplingEnabledFieldYR)
		return;

	// Scalar quantities
	for(unsigned int i = 0; i < _nNumValsFieldYR; ++i)
	{
		_nNumMoleculesFieldYRLocal[i] = 0;
	}
}

void SampleRegion::UpdateSlabParameters()
{
	mardyn_assert(0 > 1);
	return;  // Do not update these parameters by now. TODO: Get rid of this??

	double dWidth = this->GetWidth(1);
	double* dLowerCorner = this->GetLowerCorner();

	// profiles
	if(true == _SamplingEnabledProfiles)
	{
		_nNumBinsProfiles = round(dWidth / _dBinWidthProfilesInit);
		_dBinWidthProfiles = dWidth / ( (double)(_nNumBinsProfiles) );

		// recalculate Bin midpoint positions
		for(unsigned int s = 0; s < _nNumBinsProfiles; s++)
			_dBinMidpointsProfiles[s] = (s + 0.5) * _dBinWidthProfiles + dLowerCorner[1];
	}

	// VDF
	if(true == _SamplingEnabledVDF)
	{
		_numBinsVDF = round(dWidth / _dBinWidthVDFInit);
		_dBinWidthVDF = dWidth / ( (double)(_numBinsVDF) );

		// recalculate Bin midpoint positions
		for(unsigned int s = 0; s < _numBinsVDF; s++)
			_dBinMidpointsVDF[s] = (s + 0.5) * _dBinWidthVDF + dLowerCorner[1];
	}
}

// class RegionSampling

RegionSampling::RegionSampling(Domain* domain, DomainDecompBase* domainDecomp)
: ControlInstance(domain, domainDecomp)
{
}


RegionSampling::~RegionSampling()
{
	// free memory
	for(auto&& rptr : _vecSampleRegions)
	{
		delete rptr;
		rptr = nullptr;
	}
}

void RegionSampling::readXML(XMLfileUnits& xmlconfig)
{
	// add regions
	uint32_t numRegions = 0;
	uint32_t nRegID = 0;
	XMLfile::Query query = xmlconfig.query("region");
	numRegions = query.card();
	global_log->info() << "RegionSampling: Number of sampling regions: " << numRegions << endl;
	if(numRegions < 1) {
		global_log->warning() << "RegionSampling: No region parameters specified. Program exit ..." << endl;
		Simulation::exit(-1);
	}
	string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query::const_iterator outputRegionIter;
	for( outputRegionIter = query.begin(); outputRegionIter; outputRegionIter++ )
	{
		xmlconfig.changecurrentnode( outputRegionIter );
		double lc[3];
		double uc[3];
		std::string strVal[3];
		std::string strControlType;

		// coordinates
		xmlconfig.getNodeValue("coords/lcx", lc[0]);
		xmlconfig.getNodeValue("coords/lcy", lc[1]);
		xmlconfig.getNodeValue("coords/lcz", lc[2]);
		xmlconfig.getNodeValue("coords/ucx", strVal[0]);
		xmlconfig.getNodeValue("coords/ucy", strVal[1]);
		xmlconfig.getNodeValue("coords/ucz", strVal[2]);
		// read upper corner
		for(uint8_t d=0; d<3; ++d)
			uc[d] = (strVal[d] == "box") ? GetDomain()->getGlobalLength(d) : atof(strVal[d].c_str() );

		global_log->info() << "RegionSampling->region["<<nRegID<<"]: lower corner: " << lc[0] << ", " << lc[1] << ", " << lc[2] << endl;
		global_log->info() << "RegionSampling->region["<<nRegID<<"]: upper corner: " << uc[0] << ", " << uc[1] << ", " << uc[2] << endl;

		// add regions
		SampleRegion* region = new SampleRegion(this, lc, uc);
		this->AddRegion(region);

		// observer mechanism
		uint32_t refCoordsID[6] = {0, 0, 0, 0, 0, 0};
		xmlconfig.getNodeValue("coords/lcx@refcoordsID", refCoordsID[0]);
		xmlconfig.getNodeValue("coords/lcy@refcoordsID", refCoordsID[1]);
		xmlconfig.getNodeValue("coords/lcz@refcoordsID", refCoordsID[2]);
		xmlconfig.getNodeValue("coords/ucx@refcoordsID", refCoordsID[3]);
		xmlconfig.getNodeValue("coords/ucy@refcoordsID", refCoordsID[4]);
		xmlconfig.getNodeValue("coords/ucz@refcoordsID", refCoordsID[5]);

		bool bIsObserver = (refCoordsID[0]+refCoordsID[1]+refCoordsID[2]+refCoordsID[3]+refCoordsID[4]+refCoordsID[5]) > 0;

		if(true == bIsObserver)
		{
			region->PrepareAsObserver(refCoordsID);

			if(global_simulation->GetDistControl() != NULL)
				global_simulation->GetDistControl()->registerObserver(region);
			else
			{
				global_log->error() << "RegionSampling->region["<<region->GetID()<<"]: Initialization of feature DistControl is needed before! Program exit..." << endl;
				exit(-1);
			}
		}

		region->readXML(xmlconfig);
		nRegID++;

	}  // for( outputRegionIter = query.begin(); outputRegionIter; outputRegionIter++ )
}

void RegionSampling::AddRegion(SampleRegion* region)
{
	_vecSampleRegions.push_back(region);
}

void RegionSampling::Init()
{
	// init data structures
	std::vector<SampleRegion*>::iterator it;

	for(it=_vecSampleRegions.begin(); it!=_vecSampleRegions.end(); ++it)
	{
		(*it)->InitSamplingProfiles(RS_DIMENSION_Y);
		(*it)->InitSamplingVDF(RS_DIMENSION_Y);
		(*it)->InitSamplingFieldYR(RS_DIMENSION_Y);
	}
}

void RegionSampling::DoSampling(Molecule* mol, DomainDecompBase* domainDecomp, unsigned long simstep)
{
	// sample profiles and vdf
	std::vector<SampleRegion*>::iterator it;

	for(it=_vecSampleRegions.begin(); it!=_vecSampleRegions.end(); ++it)
	{
		(*it)->SampleProfiles(mol, RS_DIMENSION_Y);
		(*it)->SampleVDF(mol, RS_DIMENSION_Y);
		(*it)->SampleFieldYR(mol);
	}
}

void RegionSampling::WriteData(DomainDecompBase* domainDecomp, unsigned long simstep, Domain* domain)
{
	// write out profiles and vdf
	std::vector<SampleRegion*>::iterator it;

	for(it=_vecSampleRegions.begin(); it!=_vecSampleRegions.end(); ++it)
	{
		(*it)->WriteDataProfiles(domainDecomp, simstep, domain);
		(*it)->WriteDataVDF(domainDecomp, simstep);
		(*it)->WriteDataFieldYR(domainDecomp, simstep, domain);
	}
}

void RegionSampling::PrepareRegionSubdivisions()
{
	std::vector<SampleRegion*>::iterator it;

	for(it=_vecSampleRegions.begin(); it!=_vecSampleRegions.end(); ++it)
	{
		(*it)->PrepareSubdivisionProfiles();
		(*it)->PrepareSubdivisionVDF();
		(*it)->PrepareSubdivisionFieldYR();
	}
}

