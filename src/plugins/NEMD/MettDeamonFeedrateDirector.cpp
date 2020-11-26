#include "MettDeamonFeedrateDirector.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"
#include "utils/FileUtils.h"
#include "plugins/Mirror.h"
#include "plugins/NEMD/MettDeamon.h"

#include <cstdlib>
#include <cstdint>
#include <list>
#include <numeric>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;
using Log::global_log;

MettDeamonFeedrateDirector::MettDeamonFeedrateDirector()
	:
	_mirror_id(100)
{
	uint32_t numComponents = global_simulation->getEnsemble()->getComponents()->size();
	_particleManipCount.deleted.local.resize(numComponents+1);
	_particleManipCount.deleted.global.resize(numComponents+1);
	_particleManipCount.reflected.local.resize(numComponents+1);
	_particleManipCount.reflected.global.resize(numComponents+1);

	// reset local values
	this->resetLocalValues();
}

MettDeamonFeedrateDirector::~MettDeamonFeedrateDirector()
{
}

void MettDeamonFeedrateDirector::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain)
{
	global_log->debug() << "MettDeamonFeedrateDirector enabled." << std::endl;

	// set actual feedrate
	MettDeamon* mettDeamon = nullptr;
	std::list<PluginBase*>& plugins = *(global_simulation->getPluginList() );
	for (auto&& pit:plugins) {
		std::string name = pit->getPluginName();
		if(name == "MettDeamon")
			mettDeamon = dynamic_cast<MettDeamon*>(pit);
	}
	if(nullptr != mettDeamon) {
		// init _feedrate.sum
		_feedrate.sum = 0;
		for (std::list<double>::iterator it=_feedrate.list.begin(); it != _feedrate.list.end(); ++it) {
			_feedrate.sum += *it;
		}
		_feedrate.avg = _feedrate.sum * 1./(double)(_feedrate.list.size());
		mettDeamon->setActualFeedrate(_feedrate.avg);
	}
}

void MettDeamonFeedrateDirector::readXML(XMLfileUnits& xmlconfig)
{
	// Mirror pluginID
	_mirror_id = 100;
	xmlconfig.getNodeValue("mirror/pluginID", _mirror_id);

	// Update control
	_updateControl.sampledTimestepCount = 0;
	_updateControl.updateFreq = 1000;
	xmlconfig.getNodeValue("mirror/control/update_freq", _updateControl.updateFreq);
	global_log->info() << "MettDeamonFeedrateDirector: update frequency of Mirror = " << _updateControl.updateFreq << endl;

	// feedrate
	_feedrate.init = 0.;
	_feedrate.actual = 0.;
	_feedrate.sum = 0.;
	_feedrate.list.clear();
	_feedrate.numvals = 1;
	_feedrate.avg = 0.;
	xmlconfig.getNodeValue("mirror/feedrate/numvals", _feedrate.numvals);
	_feedrate.list.resize(_feedrate.numvals);
	bool bRet = false;
	bRet = xmlconfig.getNodeValue("mirror/feedrate/init", _feedrate.init);
	if(bRet) {
		std::fill(_feedrate.list.begin(), _feedrate.list.end(), _feedrate.init);
		_feedrate.actual = _feedrate.init;
	}
	// after restart
	std::string strCSV = "0";
	bRet = false;
	bRet = xmlconfig.getNodeValue("mirror/feedrate/list", strCSV);
	if(bRet) {
		_feedrate.list.clear();
		this->csv_str2list(strCSV, _feedrate.list);
	}
#ifndef NDEBUG
	cout << "feedrate.list:" << endl;
	for (std::list<double>::iterator it=_feedrate.list.begin(); it != _feedrate.list.end(); ++it)
		std::cout << ' ' << *it;
	cout << endl;
#endif
	// init actual feed rate
	_feedrate.actual = _feedrate.list.back();

//	_forceConstant = 100.;
//	xmlconfig.getNodeValue("forceConstant", _forceConstant);
//	global_log->info() << "MettDeamonFeedrateDirector: force constant = " << _forceConstant << endl;

	// restart information
	_restart.writefreq = 1000;
	xmlconfig.getNodeValue("restart/writefreq", _restart.writefreq);
}

void MettDeamonFeedrateDirector::beforeForces(
		ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		unsigned long simstep)
{
	Mirror* mirror = nullptr;
	MettDeamon* mettDeamon = nullptr;
	std::list<PluginBase*>& plugins = *(global_simulation->getPluginList() );
	for (auto&& pit:plugins) {
		std::string name = pit->getPluginName();
		if(name == "Mirror" && mirror == nullptr) {
			mirror = dynamic_cast<Mirror*>(pit);
			if(mirror->getPluginID() != _mirror_id)
				mirror = nullptr;
		}
		if(name == "MettDeamon")
			mettDeamon = dynamic_cast<MettDeamon*>(pit);
	}

	// Check if other plugins were found
	if(nullptr == mirror) {
		global_log->error() << "No Mirror plugin found in plugin list. Program exit ..." << std::endl;
		Simulation::exit(-2004);
	}
	if(nullptr == mettDeamon) {
		global_log->error() << "No MettDeamon plugin found in plugin list. Program exit ..." << std::endl;
		Simulation::exit(-2004);
	}

	// Get number of deleted/reflected particles from Mirror plugin
	uint32_t cid = 0;
	_particleManipCount.deleted.local.at(cid) += mirror->getDeletedParticlesCountLocal(cid);
	_particleManipCount.reflected.local.at(cid) += mirror->getReflectedParticlesCountLocal(cid);
	_updateControl.sampledTimestepCount++;

	// Calc and update new feedrate for MettDeamon plugin
	if(_updateControl.sampledTimestepCount == _updateControl.updateFreq)
	{
		_updateControl.sampledTimestepCount = 0;  // reset sampling control
		this->calcFeedrate(mettDeamon);
		mettDeamon->setActualFeedrate(_feedrate.avg);
	}

	// Write out restart information
	this->writeRestartfile();
}

void MettDeamonFeedrateDirector::afterForces(
        ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
        unsigned long simstep
)
{
}

void MettDeamonFeedrateDirector::calcFeedrate(MettDeamon* mettDeamon)
{
	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();
	uint32_t cid = 0;
	domainDecomp.collCommInit(1);
	domainDecomp.collCommAppendUnsLong(_particleManipCount.deleted.local.at(cid) );
	domainDecomp.collCommAllreduceSum();
	_particleManipCount.deleted.global.at(cid) = domainDecomp.collCommGetUnsLong();
	domainDecomp.collCommFinalize();

	// reset local values
	this->resetLocalValues();

	double dInvSampledTimestepCount = 1. / (double)(_updateControl.updateFreq);
	double deletedParticlesPerTimestep = _particleManipCount.deleted.global.at(cid) * dInvSampledTimestepCount;
	_feedrate.actual = deletedParticlesPerTimestep * mettDeamon->getInvDensityArea();

	// calc avg over prior values in _feedrate.list
	_feedrate.sum += _feedrate.actual;
	_feedrate.sum -= _feedrate.list.front();
	_feedrate.list.push_back(_feedrate.actual);
	if(_feedrate.list.size() > _feedrate.numvals)
		_feedrate.list.pop_front();
	else
		_feedrate.sum = std::accumulate(_feedrate.list.begin(), _feedrate.list.end(), 0.0);
	double dInvNumvals = 1./(double)(_feedrate.list.size());
	_feedrate.avg = _feedrate.sum * dInvNumvals;

#ifndef NDEBUG
	cout << "feedrate.list:" << endl;
	for (std::list<double>::iterator it=_feedrate.list.begin(); it != _feedrate.list.end(); ++it)
		std::cout << ' ' << *it;
	cout << endl;
	cout << "_particleManipCount.deleted.local.at(cid)=" << _particleManipCount.deleted.local.at(cid) << endl;
	cout << "_particleManipCount.deleted.global.at(cid)=" << _particleManipCount.deleted.global.at(cid) << endl;
	cout << "deletedParticlesPerTimestep=" << deletedParticlesPerTimestep << endl;
	cout << "_feedrate.actual=" << _feedrate.actual << endl;
	cout << "_feedrate.sum=" << _feedrate.sum << endl;
	cout << "_feedrate.avg=" << _feedrate.avg << endl;
#endif
}

void MettDeamonFeedrateDirector::resetLocalValues()
{
	for(auto& it:_particleManipCount.deleted.local)
		it = 0;

	for(auto& it:_particleManipCount.reflected.local)
		it = 0;
}

void MettDeamonFeedrateDirector::csv_str2list(const std::string& strCSV, std::list<double>& list)
{
    std::stringstream ss(strCSV);

    for (double i; ss >> i;) {
        list.push_back(i);
        if (ss.peek() == ',')
            ss.ignore();
    }
}

void MettDeamonFeedrateDirector::writeRestartfile()
{
	uint64_t simstep = global_simulation->getSimulationStep();
	if(0 != simstep % _restart.writefreq)
		return;

	DomainDecompBase& domainDecomp = global_simulation->domainDecomposition();

	if(domainDecomp.getRank() != 0)
		return;

	// write restart info in XML format
	{
		std::stringstream fnamestream;
		fnamestream << "MettDeamonFeedrateDirectorRestart" << "_TS" << fill_width('0', 9) << simstep << ".xml";
		std::ofstream ofs(fnamestream.str().c_str(), std::ios::out);
		std::stringstream outputstream;
		ofs << "<?xml version='1.0' encoding='UTF-8'?>" << endl;
		ofs << "<feedrate>" << endl;
		ofs << "\t<numvals>" << _feedrate.numvals << "</numvals>" << endl;
		ios::fmtflags f( ofs.flags() );
		ofs << "\t<list>" << FORMAT_SCI_MAX_DIGITS_WIDTH_21 << _feedrate.list.front();
		std::list<double>::iterator it=_feedrate.list.begin();
		for(std::advance(it, 1); it!=_feedrate.list.end(); ++it)
			ofs << "," << FORMAT_SCI_MAX_DIGITS_WIDTH_21 << *it;
		ofs << "</list>" << endl;
		ofs.flags(f);  // restore default format flags
		ofs << "</feedrate>" << endl;

		ofs << outputstream.str();
		ofs.close();
	}
}
