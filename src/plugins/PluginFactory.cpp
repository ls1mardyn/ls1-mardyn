//
// Created by Moritz Kruegener on 4/17/2018.
//

#include "PluginFactory.h"
#include "Domain.h"
#include "Simulation.h"

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
#include "io/VISWriter.h"
#include "io/XyzWriter.h"

// General plugins
#include "plugins/COMaligner.h"
#include "plugins/DirectedPM.h"
#include "plugins/Dropaccelerator.h"
#include "plugins/Dropaligner.h"
#include "plugins/ExamplePlugin.h"
#include "plugins/FixRegion.h"
#include "plugins/InMemoryCheckpointing.h"
#include "plugins/MaxCheck.h"
#include "plugins/Mirror.h"
#include "plugins/MirrorSystem.h"
#include "plugins/NEMD/RegionSampling.h"
#include "plugins/Permittivity.h"
#include "plugins/NEMD/MettDeamon.h"
#include "plugins/NEMD/MettDeamonFeedrateDirector.h"
#include "plugins/NEMD/PosNegComp.h"
#include "plugins/NEMD/DistControl.h"
#include "plugins/NEMD/DriftCtrl.h"
#include "plugins/NEMD/ExtractPhase.h"
#include "plugins/SpatialProfile.h"
#include "plugins/TestPlugin.h"
#include "plugins/VectorizationTuner.h"
#include "plugins/WallPotential.h"

#ifdef VTK
#include "io/vtk/VTKGridWriter.h"
#include "io/vtk/VTKMoleculeWriter.h"
#endif

/** @brief Register all default plugins with base PluginBase
 *
 * @param createInstance  pointer to a function returning an instance of the plugin object.
 */
template <>
void PluginFactory<PluginBase>::registerDefaultPlugins() {
	global_log->debug() << "REGISTERING PLUGINS" << endl;

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
	REGISTER_PLUGIN(MettDeamon);
	REGISTER_PLUGIN(MettDeamonFeedrateDirector);
	REGISTER_PLUGIN(PosNegComp);
	REGISTER_PLUGIN(DistControl);
	REGISTER_PLUGIN(DriftCtrl);
	REGISTER_PLUGIN(ExtractPhase);
	REGISTER_PLUGIN(ResultWriter);
	REGISTER_PLUGIN(SysMonOutput);
	REGISTER_PLUGIN(TestPlugin);
	REGISTER_PLUGIN(TimerWriter);
	REGISTER_PLUGIN(VectorizationTuner);
	REGISTER_PLUGIN(VISWriter);
	REGISTER_PLUGIN(WallPotential);
	REGISTER_PLUGIN(XyzWriter);
#ifdef VTK
	REGISTER_PLUGIN(VTKMoleculeWriter);
#ifndef MARDYN_AUTOPAS
	REGISTER_PLUGIN(VTKGridWriter);
#endif
#endif
}

/** @brief Enable selected plugins */
template <>
long PluginFactory<PluginBase>::enablePlugins(std::list<PluginBase*>& _plugins, XMLfileUnits& xmlconfig,
											  std::string category, Domain* _domain) {
	string oldpath = xmlconfig.getcurrentnodepath();

	// plugins
	long numPlugins = 0;
	XMLfile::Query query = xmlconfig.query(category);
	numPlugins = query.card();
	global_log->info() << "Number of plugins with tag " << category << ": " << numPlugins << endl;
	if (numPlugins < 1) {
		global_log->warning() << "No plugins specified for tag " << category << "." << endl;
	}

	for (auto pluginIter = query.begin(); pluginIter; ++pluginIter) {
		xmlconfig.changecurrentnode(pluginIter);
		string pluginname("");
		xmlconfig.getNodeValue("@name", pluginname);
		bool enabled = true;
		xmlconfig.getNodeValue("@enabled", enabled);
		if (not enabled) {
			global_log->debug() << "skipping disabled plugin: " << pluginname << endl;
			continue;
		}
		global_log->info() << "Enabling plugin: " << pluginname << endl;

		PluginBase* plugin = this->create(pluginname);
		if (plugin == nullptr) {
			global_log->warning() << "Could not create plugin using factory: " << pluginname << endl;
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
				global_log->error() << "[MMPLD Writer] Unknown sphere representation type: " << sphere_representation
									<< endl;
				Simulation::exit(-1);
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
			global_log->warning() << "Unknown plugin " << pluginname << endl;
		}
	}

	xmlconfig.changecurrentnode(oldpath);

	return numPlugins;
}
