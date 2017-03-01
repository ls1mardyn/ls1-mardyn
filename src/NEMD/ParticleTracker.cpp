/*
 * ParticleTracker.cpp
 *
 *  Created on: 27.01.2017
 *      Author: mheinen
 */

#include "NEMD/ParticleTracker.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>


// class TrackerState
void TrackerState::ChangeState(ParticleTracker* tracker, TrackerState* state)
{
	tracker->ChangeState(state);
}


// class TStateSleep
TrackerState* TStateSleep::_exemplar = NULL;

TrackerState* TStateSleep::Exemplar()
{
	if(_exemplar == NULL)
		_exemplar = new TStateSleep;

	return _exemplar;
}

void TStateSleep::PreLoopAction(ParticleTracker* context, unsigned long simstep)
{
}
void TStateSleep::LoopAction(ParticleTracker* context, Molecule* mol)
{
}
void TStateSleep::PostLoopAction(ParticleTracker* context)
{
}

// class TStateCollectMolIDs
TrackerState* TStateCollectMolIDs::_exemplar = NULL;

TrackerState* TStateCollectMolIDs::Exemplar()
{
	if(_exemplar == NULL)
		_exemplar = new TStateCollectMolIDs;

	return _exemplar;
}

void TStateCollectMolIDs::PreLoopAction(ParticleTracker* context, unsigned long simstep)
{
	context->SetSimstep(simstep);

	if(false == context->SimstepInsideRange() )
	{
		context->WriteToFileIDs();
		context->ChangeState(TStateSleep::Exemplar() );
	}
}
void TStateCollectMolIDs::LoopAction(ParticleTracker* context, Molecule* mol)
{
	context->CollectIDs(mol);
}
void TStateCollectMolIDs::PostLoopAction(ParticleTracker* context)
{
}

// class TStateCollectData
TrackerState* TStateCollectData::_exemplar = NULL;

TrackerState* TStateCollectData::Exemplar()
{
	if(_exemplar == NULL)
		_exemplar = new TStateCollectData;

	return _exemplar;
}

void TStateCollectData::PreLoopAction(ParticleTracker* context, unsigned long simstep)
{
	context->SetSimstep(simstep);

	// reset sample map
	context->ResetSampleMap();

	if(false == context->SimstepInsideRange() )
		context->ChangeState(TStateSleep::Exemplar() );
}
void TStateCollectData::LoopAction(ParticleTracker* context, Molecule* mol)
{
	context->SampleData(mol);
}
void TStateCollectData::PostLoopAction(ParticleTracker* context)
{
	context->AppendDataMap();

	if(context->GetSimstep() % context->GetWriteFreqData() != 0)
		return;

	context->CalcGlobalValuesData();
	context->WriteToFileData();
	context->ClearSampleDataVector();
}


ParticleTracker::ParticleTracker( DomainDecompBase* domainDecomp, unsigned long simstepStart,
		unsigned long simstepStop) : _domainDecomp(domainDecomp), _simstep(0), _writeFreqData(1000),
		_simstepStart(simstepStart), _simstepStop(simstepStop)
{
	_trackerState = TStateSleep::Exemplar();
}

ParticleTracker::~ParticleTracker()
{
}

void ParticleTracker::Prepare()
{
	unsigned long numIDs = this->ReadMoleculeIDs();

	if(numIDs > 0)
	{
		this->CreateDataFiles();
		this->InitDatastructures();
		this->ChangeState(TStateCollectData::Exemplar() );
	}
	else
	{
		this->CreateMoleculeIDsFile();
		this->ChangeState(TStateCollectMolIDs::Exemplar() );
	}

}

void ParticleTracker::PreLoopAction(unsigned long simstep)
{
	_trackerState->PreLoopAction(this, simstep);
}

void ParticleTracker::LoopAction(Molecule* mol)
{
	_trackerState->LoopAction(this, mol);
}

void ParticleTracker::PostLoopAction()
{
	_trackerState->PostLoopAction(this);
}

unsigned int ParticleTracker::ReadMoleculeIDs()
{
	_vecMoleculeIDs.clear();

	std::stringstream sstrFilename;
	sstrFilename << "moleculeIDs.dat";
	std::ifstream ifs(sstrFilename.str().c_str(), std::ios::out);

	if( ifs.is_open() )
	{
		while( !ifs.eof() )
		{
			std::string strToken;
			getline(ifs, strToken);
			unsigned long mid = atoi(strToken.c_str() );
			if(mid > 0)
				_vecMoleculeIDs.push_back(mid);
		}
		ifs.close();
	}
	return _vecMoleculeIDs.size();
}

void ParticleTracker::ShowSampleMap()
{
	if(_domainDecomp->getRank() != 0)
		return;

	std::cout << std::endl;
	std::cout << "map content:" << std::endl;
	for(auto& mit : _mapMolIDsDoubleVals)
	{
		std::vector<double>::iterator it;
		std::vector<double>& vec = mit.second;
		for(it=vec.begin(); it!=vec.end(); ++it)
		{
			std::cout << std::setw(8) << std::setprecision(3) << *it << " | ";
		}
		std::cout << std::endl;
	}
}

void ParticleTracker::ResetSampleMap()
{
	for(auto& mit : _mapMolIDsDoubleVals)
	{
	std::vector<double>::iterator it;
		std::vector<double>& vec = mit.second;
		for(it=vec.begin(); it!=vec.end(); ++it)
		{
			*it = 0.0;
		}
	}
}

void ParticleTracker::CreateDataFiles()
{
	// create sampling files
	_vecQuantyNamesDouble = {
		"U_kin", "U_trans", "U_rot",
		"F[0]", "F[1]", "F[2]", "Fabs",
		"M[0]", "M[1]", "M[2]", "Mabs",
		"D[0]", "D[1]", "D[2]", "Dabs",
	};

	// files only created once by root process
	if(_domainDecomp->getRank() != 0)
		return;

	// create filename vector
	std::vector<std::string> vecFileNames;
	vecFileNames.push_back("cid");

	std::vector<std::string>::iterator it;
	for(auto& it : vecFileNames)
	{
		std::string strQuantyNameDouble = it;
		vecFileNames.push_back(strQuantyNameDouble);
	}

	for(auto& it : vecFileNames)
	{
		std::stringstream sstrFilename;
		sstrFilename << it << ".dat";
		std::ofstream ofs( sstrFilename.str().c_str(), std::ios::out);
		std::stringstream outputstream;
		for(auto id : _vecMoleculeIDs)
		{
			std::stringstream sstrColHeader;
			sstrColHeader << "mol" << id;
			outputstream << std::setw(24) << sstrColHeader.str();
		}
		outputstream << std::endl;
		ofs << outputstream.str();
		ofs.close();
	}
}

void ParticleTracker::CreateMoleculeIDsFile()
{
	if(_domainDecomp->getRank() != 0) // create file
		return;

	std::ofstream ofs( "moleculeIDs.dat", std::ios::out);
	ofs << "             moleculeIDs" << std::endl;
	ofs.close();
}

void ParticleTracker::InitDatastructures()
{
	for(auto id : _vecMoleculeIDs)
	{
		std::vector<double> vec(_vecQuantyNamesDouble.size(), 0.0);
		_mapMolIDsDoubleVals.insert(std::pair<unsigned long, std::vector<double> >(id, vec) );
#ifndef NDEBUG
		std::cout << "id = " << id << std::endl;
#endif
	}
}

void ParticleTracker::CalcGlobalValuesData()
{
	// calc global values
	size_t numDoubleVals = _vecSampleDataDouble.begin()->begin()->second.size();
	size_t numMolIDs = _vecSampleDataDouble.begin()->size();
	size_t numSimsteps = _vecSampleDataDouble.size();

	_domainDecomp->collCommInit(numDoubleVals*numMolIDs*numSimsteps);

#ifndef NDEBUG
	std::cout << "[" << _domainDecomp->getRank() << "]: numDoubleVals = " << numDoubleVals << std::endl;
	std::cout << "[" << _domainDecomp->getRank() << "]: numMolIDs = " << numMolIDs << std::endl;
	std::cout << "[" << _domainDecomp->getRank() << "]: numSimsteps = " << numSimsteps << std::endl;
#endif

	for(auto vit=_vecSampleDataDouble.begin(); vit!=_vecSampleDataDouble.end(); ++vit)
	{
		for(auto mit=vit->begin(); mit!=vit->end(); ++mit)
		{
			for(auto vit2=mit->second.begin(); vit2!=mit->second.end(); ++vit2)
			{
				_domainDecomp->collCommAppendDouble(*vit2);
			}
		}
	}
	// reduce
	_domainDecomp->collCommAllreduceSum();
	for(auto vit=_vecSampleDataDouble.begin(); vit!=_vecSampleDataDouble.end(); ++vit)
	{
		for(auto mit=vit->begin(); mit!=vit->end(); ++mit)
		{
			for(auto vit2=mit->second.begin(); vit2!=mit->second.end(); ++vit2)
			{
				(*vit2) = _domainDecomp->collCommGetDouble();
			}
		}
	}
	_domainDecomp->collCommFinalize();
}

void ParticleTracker::WriteToFileData()
{
	// write to files
	if(_domainDecomp->getRank() == 0)
	{
		std::cout << "write to files" << std::endl;

		size_t numDoubleVals = _vecSampleDataDouble.begin()->begin()->second.size();

		for(size_t qi=0; qi<numDoubleVals; qi++)
		{
			std::stringstream sstrFilename;
			sstrFilename << _vecQuantyNamesDouble.at(qi) << ".dat";
			std::stringstream outputstream;

			for(auto vit : _vecSampleDataDouble)
			{
				for(auto mit : vit)
				{
					outputstream << std::setw(24) << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << mit.second.at(qi);
				}
				outputstream << std::endl;
			}

			std::ofstream fileout(sstrFilename.str().c_str(), std::ios::app);
			fileout << outputstream.str();
			fileout.close();
		}

	} // if rank == 0
}

void ParticleTracker::WriteToFileIDs()
{
	if(_domainDecomp->getRank() != 0)
		return;

	std::ofstream ofs("moleculeIDs.dat", std::ios::app);
	std::stringstream outputstream;

	for(auto id : _mapCollectIDs)
	{
		outputstream << std::setw(24) << id.first << std::endl;
	}
	ofs << outputstream.str();
	ofs.close();
}

void ParticleTracker::CollectIDs(Molecule* mol)
{
	unsigned char cid = mol->componentid();

	if(cid == 2)
	{
		unsigned long id = mol->id();

		if( _mapCollectIDs.count(id) == 0)
			_mapCollectIDs.insert(std::pair<unsigned long, char>(id, cid) );
	}
}

void ParticleTracker::SampleData(Molecule* mol)
{
	unsigned long id = mol->id();

	std::map<unsigned long, std::vector<double> >::iterator mit;
	mit = _mapMolIDsDoubleVals.find(id);

	if(mit != _mapMolIDsDoubleVals.end() )
	{
		double F[3];
		double M[3];
		double D[3];

		for(auto d=0; d<3; ++d)
		{
			F[d] = mol->F(d);
			M[d] = mol->M(d);
			D[d] = mol->D(d);
		}

		std::vector<double>& vec = mit->second;
		vec.at(0) = mol->U_kin();
		vec.at(1) = mol->U_trans();
		vec.at(2) = mol->U_rot();

		vec.at(3) = F[0];
		vec.at(4) = F[0];
		vec.at(5) = F[0];
		double Fabs = F[0]*F[0] + F[1]*F[1] + F[2]*F[2];
		vec.at(6) = Fabs;

		vec.at(7) = M[0];
		vec.at(8) = M[0];
		vec.at(9) = M[0];
		double Mabs = M[0]*M[0] + M[1]*M[1] + M[2]*M[2];
		vec.at(10) = Mabs;

		vec.at(11) = D[0];
		vec.at(12) = D[0];
		vec.at(13) = D[0];
		double Dabs = D[0]*D[0] + D[1]*D[1] + D[2]*D[2];
		vec.at(14) = Dabs;
	}
}

