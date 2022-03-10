#include "DriftCtrl.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"
#include "utils/FileUtils.h"

#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;
using Log::global_log;

DriftCtrl::DriftCtrl()
{
}

DriftCtrl::~DriftCtrl()
{
}

void DriftCtrl::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain)
{
	global_log->debug() << "[DriftCtrl] Enabled in range yl,yr=" << _range.yl << "," << _range.yr << std::endl;
	
	// number of components
	uint32_t numComponents = domain->getNumberOfComponents() + 1;  // + 1 because component 0 stands for all components
	_sampling.resize(numComponents);
	for(uint32_t cid = 0; cid < numComponents; ++cid) {
		// local
		_sampling.at(cid).numParticles.local.resize(_range.subdivision.numBins);
		_sampling.at(cid).momentum.at(0).local.resize(_range.subdivision.numBins);
		_sampling.at(cid).momentum.at(1).local.resize(_range.subdivision.numBins);
		_sampling.at(cid).momentum.at(2).local.resize(_range.subdivision.numBins);
		// global
		_sampling.at(cid).numParticles.global.resize(_range.subdivision.numBins);
		_sampling.at(cid).momentum.at(0).global.resize(_range.subdivision.numBins);
		_sampling.at(cid).momentum.at(1).global.resize(_range.subdivision.numBins);
		_sampling.at(cid).momentum.at(2).global.resize(_range.subdivision.numBins);
		// corr
		_sampling.at(cid).mom_corr.at(0).resize(_range.subdivision.numBins);
		_sampling.at(cid).mom_corr.at(1).resize(_range.subdivision.numBins);
		_sampling.at(cid).mom_corr.at(2).resize(_range.subdivision.numBins);
	}
	// reset local values
	for(uint32_t cid = 0; cid < numComponents; ++cid)
	{
		for(uint32_t yPosID = 0; yPosID < _range.subdivision.numBins; ++yPosID)
		{
			_sampling.at(cid).numParticles.local.at(yPosID) = 0;
			_sampling.at(cid).momentum.at(0).local.at(yPosID) = 0.;
			_sampling.at(cid).momentum.at(1).local.at(yPosID) = 0.;
			_sampling.at(cid).momentum.at(2).local.at(yPosID) = 0.;
		}
	}
	global_log->info() << "[DriftCtrl] Init data structures for " << numComponents << " components." << std::endl;
	
	// init files
	uint64_t simstep = global_simulation->getSimulationStep();
	{
		const std::string fname = "DriftCtrl_drift.dat";
		std::ofstream ofs;
		ofs.open(fname, std::ios::out);
		ofs << "     simstep";
		for(uint32_t yPosID = 0; yPosID < _range.subdivision.numBins; ++yPosID)
			ofs  << "                 bin" << fill_width('0', 4) << yPosID;
		ofs << std::endl;
		ofs.close();
	}
	{
		const std::string fname = "DriftCtrl_numParticles.dat";
		std::ofstream ofs;
		ofs.open(fname, std::ios::out);
		ofs << "     simstep";
		for(uint32_t yPosID = 0; yPosID < _range.subdivision.numBins; ++yPosID)
			ofs  << "     bin" << fill_width('0', 4) << yPosID;
		ofs << std::endl;
		ofs.close();
	}
}

void DriftCtrl::readXML(XMLfileUnits& xmlconfig)
{
	// control
	_control.freq.sample = 10;
	_control.freq.control = 100;
	_control.freq.write = 10000;
	xmlconfig.getNodeValue("control/freq/sample", _control.freq.sample);
	xmlconfig.getNodeValue("control/freq/control", _control.freq.control);
	xmlconfig.getNodeValue("control/freq/write", _control.freq.write);
	
	// range
	_range.yl = 0.;
	_range.yr = 100.;
	xmlconfig.getNodeValue("range/yl", _range.yl);
	xmlconfig.getNodeValue("range/yr", _range.yr);
	global_log->info() << "[DriftCtrl] Enabled in range yl,yr=" << _range.yl << "," << _range.yr << std::endl;
	_range.width = _range.yr - _range.yl;
	// subdivision
	_range.subdivision.binWidth.init = 10.;
	_range.subdivision.binWidth.actual = 10.;
	_range.subdivision.numBins = 10;
	xmlconfig.getNodeValue("range/subdivision/binwidth", _range.subdivision.binWidth.init);
	uint32_t numBins = 0;
	numBins = static_cast<uint32_t>(_range.width / _range.subdivision.binWidth.init);
	_range.subdivision.binWidth.actual = _range.width / static_cast<double>(numBins);
	_range.subdivision.numBins = numBins;

	// target drift
	_target.drift.at(0) = 0.;
	_target.drift.at(1) = 0.;
	_target.drift.at(2) = 0.;
	_target.cid = 1;
	xmlconfig.getNodeValue("target/cid", _target.cid);
	xmlconfig.getNodeValue("target/drift/vx", _target.drift.at(0));
	xmlconfig.getNodeValue("target/drift/vy", _target.drift.at(1));
	xmlconfig.getNodeValue("target/drift/vz", _target.drift.at(2));
	global_log->info() << "[DriftCtrl] Target drift vx,vy,vz="
		<< _target.drift.at(0) << "," << _target.drift.at(1) << "," << _target.drift.at(2) << ", cid=" << _target.cid<< endl;
}

void DriftCtrl::beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep)
{
	int nRank = domainDecomp->getRank();
	
	// sample
	if(simstep % _control.freq.sample == 0)
	{
		for(auto it = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it) {
			// check if inside range
			double yPos = it->r(1);
			if(yPos <= _range.yl || yPos > _range.yr)
				continue;
			
			const uint32_t cid_zb = it->componentid();
			const uint32_t cid_ub = cid_zb+1;
			const uint32_t yPosID = floor( (yPos-_range.yl) / _range.subdivision.binWidth.actual);
			const double mass = it->mass();
			std::array<double,3> p;  // momentum
			p.at(0) = it->v(0) * mass;
			p.at(1) = it->v(1) * mass;
			p.at(2) = it->v(2) * mass;

			_sampling.at(cid_ub).numParticles.local.at(yPosID)++;
			_sampling.at(cid_ub).momentum.at(0).local.at(yPosID) += p.at(0);
			_sampling.at(cid_ub).momentum.at(1).local.at(yPosID) += p.at(1);
			_sampling.at(cid_ub).momentum.at(2).local.at(yPosID) += p.at(2);
		}
	}
	
	// control
	if(simstep % _control.freq.control == 0)
	{
		uint32_t numComponents = _sampling.size();
		uint64_t numVals = static_cast<uint64_t>(numComponents) * 4 * _range.subdivision.numBins;
		uint64_t numValsCheck = 0;
		domainDecomp->collCommInit(numVals);
		// append local values
		for(uint32_t cid = 0; cid < numComponents; ++cid)
		{
			for(uint32_t yPosID = 0; yPosID < _range.subdivision.numBins; ++yPosID)
			{
				domainDecomp->collCommAppendUnsLong(_sampling.at(cid).numParticles.local.at(yPosID));
				domainDecomp->collCommAppendDouble(_sampling.at(cid).momentum.at(0).local.at(yPosID));
				domainDecomp->collCommAppendDouble(_sampling.at(cid).momentum.at(1).local.at(yPosID));
				domainDecomp->collCommAppendDouble(_sampling.at(cid).momentum.at(2).local.at(yPosID));
				numValsCheck += 4;
			}
		}
		//~ cout << "numVals ?= numValsCheck: " << numVals << "?=" << numValsCheck << endl;
		// reduce
		domainDecomp->collCommAllreduceSum();
		// collect global values
		for(uint32_t cid = 0; cid < numComponents; ++cid)
		{
			for(uint32_t yPosID = 0; yPosID < _range.subdivision.numBins; ++yPosID)
			{
				uint64_t numParticles = domainDecomp->collCommGetUnsLong();
				if(numParticles < 1)
					numParticles = 1;
				double invNumParticles = 1./static_cast<double>(numParticles);
				_sampling.at(cid).numParticles.global.at(yPosID) = numParticles;
				//~ cout << "[" << nRank << "]: cid=" << cid << ",yPosID=" << yPosID << ",numParticles=" << numParticles << endl;
				_sampling.at(cid).momentum.at(0).global.at(yPosID) = domainDecomp->collCommGetDouble() * invNumParticles;
				_sampling.at(cid).momentum.at(1).global.at(yPosID) = domainDecomp->collCommGetDouble() * invNumParticles;
				_sampling.at(cid).momentum.at(2).global.at(yPosID) = domainDecomp->collCommGetDouble() * invNumParticles;
				
				// reset local values
				_sampling.at(cid).numParticles.local.at(yPosID) = 0;
				_sampling.at(cid).momentum.at(0).local.at(yPosID) = 0.;
				_sampling.at(cid).momentum.at(1).local.at(yPosID) = 0.;
				_sampling.at(cid).momentum.at(2).local.at(yPosID) = 0.;
			}
		}
		// finalize
		domainDecomp->collCommFinalize();
		
		// calc correction
		for(uint32_t cid = 0; cid < numComponents; ++cid)
		{
			for(uint32_t yPosID = 0; yPosID < _range.subdivision.numBins; ++yPosID)
			{
				_sampling.at(cid).mom_corr.at(0).at(yPosID) = _target.drift.at(0) - _sampling.at(cid).momentum.at(0).global.at(yPosID);
				_sampling.at(cid).mom_corr.at(1).at(yPosID) = _target.drift.at(1) - _sampling.at(cid).momentum.at(1).global.at(yPosID);
				_sampling.at(cid).mom_corr.at(2).at(yPosID) = _target.drift.at(2) - _sampling.at(cid).momentum.at(2).global.at(yPosID);
			}
		}
		
		// do correction
		for(auto it = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); it.isValid(); ++it) {
			// check if inside range
			double yPos = it->r(1);
			if(yPos <= _range.yl || yPos > _range.yr)
				continue;
			
			// check if target component
			const uint32_t cid_zb = it->componentid();
			const uint32_t cid_ub = cid_zb+1;
			if ((cid_ub != _target.cid) and (_target.cid != 0)) {
				continue;
			}
			
			uint32_t yPosID = floor( (yPos-_range.yl) / _range.subdivision.binWidth.actual);
			const double mass = it->mass();
			it->setv(0, it->v(0) + _sampling.at(cid_ub).mom_corr.at(0).at(yPosID)/mass ); 
			it->setv(1, it->v(1) + _sampling.at(cid_ub).mom_corr.at(1).at(yPosID)/mass );
			it->setv(2, it->v(2) + _sampling.at(cid_ub).mom_corr.at(2).at(yPosID)/mass );
		}
	}
	
	// write to file
	if(simstep % _control.freq.write == 0 && nRank == 0)
	{
		{
			const std::string fname = "DriftCtrl_drift.dat";
			std::ofstream ofs;
			ofs.open(fname, std::ios::app);
			ofs << setw(12) << simstep;
			for(uint32_t yPosID = 0; yPosID < _range.subdivision.numBins; ++yPosID) {
				ofs << FORMAT_SCI_MAX_DIGITS << _sampling.at(_target.cid).momentum.at(1).global.at(yPosID);
			}
			ofs << std::endl;
			ofs.close();
		}
		{
			const std::string fname = "DriftCtrl_numParticles.dat";
			std::ofstream ofs;
			ofs.open(fname, std::ios::app);
			ofs << setw(12) << simstep;
			for(uint32_t yPosID = 0; yPosID < _range.subdivision.numBins; ++yPosID) {
				ofs << setw(12) << _sampling.at(_target.cid).numParticles.global.at(yPosID);
			}
			ofs << std::endl;
			ofs.close();
		}
	}
}
