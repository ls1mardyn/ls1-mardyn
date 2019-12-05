#include "MettDeamonFeedrateDirector.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"
#include "plugins/Mirror.h"
#include "plugins/NEMD/MettDeamon.h"

#include <cstdlib>
#include <cstdint>
#include <list>

using namespace std;
using Log::global_log;

MettDeamonFeedrateDirector::MettDeamonFeedrateDirector()
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

//	_forceConstant = 100.;
//	xmlconfig.getNodeValue("forceConstant", _forceConstant);
//	global_log->info() << "MettDeamonFeedrateDirector: force constant = " << _forceConstant << endl;
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
		double feedrate_old = _feedrate;
		this->calcFeedrate(mettDeamon);
		if(_feedrate != feedrate_old)
			mettDeamon->setActualFeedrate(_feedrate);
	}
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
	_feedrate = deletedParticlesPerTimestep * mettDeamon->getInvDensityArea();
	cout << "MDFRD: No. deleted particles = " << _particleManipCount.deleted.global.at(cid) << endl;
	cout << "dInvSampledTimestepCount = " << dInvSampledTimestepCount << endl;
	cout << "mettDeamon->getInvDensityArea() = " << mettDeamon->getInvDensityArea() << endl;
}

void MettDeamonFeedrateDirector::resetLocalValues()
{
	for(auto& it:_particleManipCount.deleted.local)
		it = 0;

	for(auto& it:_particleManipCount.reflected.local)
		it = 0;
}
