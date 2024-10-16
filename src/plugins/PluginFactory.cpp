//
// Created by Moritz Kruegener on 4/17/2018.
//

#include "PluginFactory.h"
#include "Domain.h"
#include "utils/mardyn_assert.h"

#include <map>
#include <string>
#include <vector>

#include "utils/Logger.h"
#include "utils/String_utils.h"

// Output plugins
#include "io/CavityWriter.h"
#include "io/CheckpointWriter.h"
#include "io/CommunicationPartnerWriter.h"
#include "io/DecompWriter.h"
#include "io/EnergyLogWriter.h"
#include "io/FlopRateWriter.h"
#include "io/GammaWriter.h"
#include "io/HaloParticleWriter.h"
#include "io/KDTreePrinter.h"
#include "io/LoadBalanceWriter.h"
#include "io/MPICheckpointWriter.h"
#include "io/MaxWriter.h"
#include "io/MmpldWriter.h"
#include "io/MmspdBinWriter.h"
#include "io/MmspdWriter.h"
#include "io/ODF.h"
#include "io/PovWriter.h"
#include "io/RDF.h"
#include "io/ResultWriter.h"
#include "io/SysMonOutput.h"
#include "io/TimerWriter.h"
#include "io/XyzWriter.h"

// General plugins
#include "plugins/COMaligner.h"
#include "plugins/DirectedPM.h"
#include "plugins/Dropaccelerator.h"
#include "plugins/Dropaligner.h"
#include "plugins/ExamplePlugin.h"
#include "plugins/FixRegion.h"
#include "plugins/InMemoryCheckpointing.h"
#include "plugins/LoadImbalanceThroughSleepPlugin.h"
#include "plugins/MaxCheck.h"
#include "plugins/Mirror.h"
#include "plugins/MirrorSystem.h"
#include "plugins/NEMD/DensityControl.h"
#include "plugins/NEMD/DistControl.h"
#include "plugins/NEMD/DriftCtrl.h"
#include "plugins/NEMD/ExtractPhase.h"
#include "plugins/NEMD/MettDeamon.h"
#include "plugins/NEMD/MettDeamonFeedrateDirector.h"
#include "plugins/NEMD/PosNegComp.h"
#include "plugins/NEMD/RegionSampling.h"
#include "plugins/NEMD/VelocityExchange.h"
#include "plugins/Permittivity.h"
#include "plugins/SpatialProfile.h"
#include "plugins/TestPlugin.h"
#include "plugins/VectorizationTuner.h"
#include "plugins/WallPotential.h"
#include "plugins/EnergyRAPL.h"
#ifdef ENABLE_ADIOS2
#include "io/Adios2Writer.h"
#endif

#ifdef VTK
#include "io/vtk/VTKGridWriter.h"
#include "io/vtk/VTKMoleculeWriter.h"
#endif

#ifdef MAMICO_COUPLING
#include "plugins/MamicoCoupling.h"
#endif

/** @brief Register all default plugins with base PluginBase
 *
 * @param createInstance  pointer to a function returning an instance of the plugin object.
 */
template <>
void PluginFactory<PluginBase>::registerDefaultPlugins() {
	Log::global_log->debug() << "REGISTERING PLUGINS" << std::endl;

#ifdef ENABLE_ADIOS2
	REGISTER_PLUGIN(Adios2Writer);
#endif
	REGISTER_PLUGIN(COMaligner);
	REGISTER_PLUGIN(CavityWriter);
	REGISTER_PLUGIN(CheckpointWriter);
	REGISTER_PLUGIN(CommunicationPartnerWriter);
	REGISTER_PLUGIN(DecompWriter);
	REGISTER_PLUGIN(DirectedPM);
	REGISTER_PLUGIN(Dropaccelerator);
	REGISTER_PLUGIN(Dropaligner);
	REGISTER_PLUGIN(EnergyLogWriter);
	REGISTER_PLUGIN(ExamplePlugin);
	REGISTER_PLUGIN(FixRegion);
	REGISTER_PLUGIN(FlopRateWriter);
	REGISTER_PLUGIN(GammaWriter);
	REGISTER_PLUGIN(HaloParticleWriter);
	REGISTER_PLUGIN(InMemoryCheckpointing);
	REGISTER_PLUGIN(SpatialProfile);
	REGISTER_PLUGIN(KDTreePrinter);
	REGISTER_PLUGIN(LoadbalanceWriter);
	REGISTER_PLUGIN(LoadImbalanceThroughSleepPlugin);
	REGISTER_PLUGIN(MPICheckpointWriter);
	REGISTER_PLUGIN(MaxCheck);
	REGISTER_PLUGIN(MaxWriter);
	REGISTER_PLUGIN(Mirror);
	REGISTER_PLUGIN(MirrorSystem);
	REGISTER_PLUGIN(MmpldWriter);
	REGISTER_PLUGIN(MmspdBinWriter);
	REGISTER_PLUGIN(MmspdWriter);
	REGISTER_PLUGIN(ODF);
	REGISTER_PLUGIN(Permittivity);
	REGISTER_PLUGIN(PovWriter);
	REGISTER_PLUGIN(RDF);
	REGISTER_PLUGIN(RegionSampling);
	REGISTER_PLUGIN(VelocityExchange);
	REGISTER_PLUGIN(MettDeamon);
	REGISTER_PLUGIN(MettDeamonFeedrateDirector);
	REGISTER_PLUGIN(PosNegComp);
	REGISTER_PLUGIN(DensityControl);
	REGISTER_PLUGIN(DistControl);
	REGISTER_PLUGIN(DriftCtrl);
	REGISTER_PLUGIN(ExtractPhase);
	REGISTER_PLUGIN(ResultWriter);
	REGISTER_PLUGIN(SysMonOutput);
	REGISTER_PLUGIN(TestPlugin);
	REGISTER_PLUGIN(TimerWriter);
	REGISTER_PLUGIN(VectorizationTuner);
	REGISTER_PLUGIN(WallPotential);
	REGISTER_PLUGIN(XyzWriter);
	REGISTER_PLUGIN(EnergyRAPL);
#ifdef VTK
	REGISTER_PLUGIN(VTKMoleculeWriter);
#ifndef MARDYN_AUTOPAS
	REGISTER_PLUGIN(VTKGridWriter);
#endif
#endif

#ifdef MAMICO_COUPLING
	REGISTER_PLUGIN(MamicoCoupling);
#endif
}

/** @brief Enable selected plugins */
template <>
long PluginFactory<PluginBase>::enablePlugins(std::list<PluginBase*>& _plugins, XMLfileUnits& xmlconfig,
											  std::string category, Domain* _domain) {
	std::string oldpath = xmlconfig.getcurrentnodepath();

	// plugins
	long numPlugins = 0;
	XMLfile::Query query = xmlconfig.query(category);
	numPlugins = query.card();
	Log::global_log->info() << "Number of plugins with tag " << category << ": " << numPlugins << std::endl;
	if (numPlugins < 1) {
		Log::global_log->warning() << "No plugins specified for tag " << category << "." << std::endl;
	}

	for (auto pluginIter = query.begin(); pluginIter; ++pluginIter) {
		xmlconfig.changecurrentnode(pluginIter);
		std::string pluginname("");
		xmlconfig.getNodeValue("@name", pluginname);
		bool enabled = true;
		xmlconfig.getNodeValue("@enabled", enabled);
		if (not enabled) {
			Log::global_log->debug() << "skipping disabled plugin: " << pluginname << std::endl;
			continue;
		}
		Log::global_log->info() << "Enabling plugin: " << pluginname << std::endl;

		PluginBase* plugin = this->create(pluginname);
		if (plugin == nullptr) {
			Log::global_log->warning() << "Could not create plugin using factory: " << pluginname << std::endl;
		}

		//@TODO: add plugin specific functions

		if (pluginname == "MmpldWriter") {
			// @todo this should be handled in the MMPLD Writer readXML()
			std::string sphere_representation = "simple";
			xmlconfig.getNodeValue("@type", sphere_representation);
			delete plugin;
			if ("simple" == sphere_representation) {
				plugin = new MmpldWriterSimpleSphere();
			} else if ("multi" == sphere_representation) {
				plugin = new MmpldWriterMultiSphere();
			} else {
				std::ostringstream error_message;
				error_message << "[MMPLD Writer] Unknown sphere representation type: " << sphere_representation << std::endl;
				MARDYN_EXIT(error_message.str());
			}
		} else if (pluginname == "DomainProfiles") {
			plugin = this->create("DensityProfileWriter");
			// TODO: add _domain access (via Simularion)
			_domain->readXML(xmlconfig);
		}

		if (nullptr != plugin) {
			plugin->readXML(xmlconfig);
			_plugins.push_back(plugin);
		} else {
			Log::global_log->warning() << "Unknown plugin " << pluginname << std::endl;
		}
	}

	xmlconfig.changecurrentnode(oldpath);

	return numPlugins;
}
