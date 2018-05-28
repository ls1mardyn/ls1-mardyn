//
// Created by Moritz Kruegener on 4/17/2018.
//

#include "PluginFactory.h"
#include "Simulation.h"
#include "Domain.h"

#include <map>
#include <string>
#include <vector>

#include "utils/Logger.h"
#include "utils/String_utils.h"

// Output plugins
#include "io/CavityWriter.h"
#include "io/CheckpointWriter.h"
#include "io/DecompWriter.h"
#include "io/DensityProfileWriter.h"
#include "io/EnergyLogWriter.h"
#include "io/FlopRateWriter.h"
#include "io/GammaWriter.h"
#include "io/LoadBalanceWriter.h"
#include "io/MPICheckpointWriter.h"
#include "io/MmpldWriter.h"
#include "io/MmspdBinWriter.h"
#include "io/MmspdWriter.h"
#include "io/PovWriter.h"
#include "io/RDF.h"
#include "io/ResultWriter.h"
#include "io/SysMonOutput.h"
#include "io/VISWriter.h"
#include "io/XyzWriter.h"
#include "io/MaxWriter.h"

// General plugins
#include "plugins/COMaligner.h"
#include "plugins/Mirror.h"

#ifdef VTK
#include "io/vtk/VTKMoleculeWriter.h"
#include "io/vtk/VTKGridWriter.h"
#endif

#include "utils/TestPlugin.h"

/** @brief Register all default plugins with base PluginBase
	 *
	 * @param createInstance  pointer to a function returning an instance of the plugin object.
	 */
template<>
void PluginFactory<PluginBase>::registerDefaultPlugins(){
    global_log -> debug() << "REGISTERING PLUGINS" << endl;
    REGISTER_PLUGIN(TestPlugin);
    REGISTER_PLUGIN(CavityWriter);
    REGISTER_PLUGIN(CheckpointWriter);
    REGISTER_PLUGIN(DecompWriter);
    REGISTER_PLUGIN(DensityProfileWriter);
    REGISTER_PLUGIN(EnergyLogWriter);
    REGISTER_PLUGIN(FlopRateWriter);
    REGISTER_PLUGIN(GammaWriter);
    REGISTER_PLUGIN(LoadbalanceWriter);
    REGISTER_PLUGIN(MPICheckpointWriter);
    REGISTER_PLUGIN(MmpldWriter);
    REGISTER_PLUGIN(MmspdBinWriter);
    REGISTER_PLUGIN(MmspdWriter);
    REGISTER_PLUGIN(PovWriter);
    REGISTER_PLUGIN(RDF);
    REGISTER_PLUGIN(ResultWriter);
    REGISTER_PLUGIN(SysMonOutput);
    REGISTER_PLUGIN(VISWriter);
    REGISTER_PLUGIN(XyzWriter);
    REGISTER_PLUGIN(MaxWriter);
    REGISTER_PLUGIN(COMaligner);
    REGISTER_PLUGIN(Mirror);

#ifdef VTK
    REGISTER_PLUGIN(VTKMoleculeWriter);
    REGISTER_PLUGIN(VTKGridWriter);
#endif
}

/** @brief Enable selected plugins */
template<>
long PluginFactory<PluginBase>::enablePlugins(std::list<PluginBase*>& _plugins, XMLfileUnits& xmlconfig, std::string category, Domain* _domain) {
    string oldpath = xmlconfig.getcurrentnodepath();

    // plugins
    long numPlugins = 0;
    XMLfile::Query query = xmlconfig.query(category);
    numPlugins = query.card();
    global_log->info() << "Number of plugins with tag " << category << ": " << numPlugins << endl;
    if(numPlugins < 1) {
        global_log->warning() << "No plugins specified for tag" << category << "." << endl;
    }

    for (auto pluginIter = query.begin(); pluginIter; ++pluginIter) {
        xmlconfig.changecurrentnode( pluginIter );
        string pluginname("");
        xmlconfig.getNodeValue("@name", pluginname);
        bool enabled = true;
        xmlconfig.getNodeValue("@enabled", enabled);
        if(not enabled) {
            global_log->debug() << "skipping disabled plugin: " << pluginname << endl;
            continue;
        }
        global_log->info() << "Enabling plugin: " << pluginname << endl;


        PluginBase* plugin = this->create(pluginname);
        if(plugin == nullptr) {
            global_log->warning() << "Could not create plugin using factory: " << pluginname << endl;
        }

        //@TODO: add plugin specific functions

        if(pluginname == "MmpldWriter") {
            // @todo this should be handled in the MMPLD Writer readXML()
            std::string sphere_representation = "simple";
            xmlconfig.getNodeValue("@type", sphere_representation);
            delete plugin;
            if("simple" == sphere_representation) {
                plugin = new MmpldWriterSimpleSphere();
            } else if("multi" == sphere_representation) {
                plugin = new MmpldWriterMultiSphere ();
            } else {
                global_log->error() << "[MMPLD Writer] Unknown sphere representation type: " << sphere_representation << endl;
                Simulation::exit(-1);
            }
        }
            // temporary
        else if(pluginname == "VectorizationTuner") {
            //plugin = new VectorizationTuner(_cutoffRadius, _LJCutoffRadius, &_cellProcessor);
        }
        else if(pluginname == "DomainProfiles") {
            plugin = this->create("DensityProfileWriter");
            // TODO: add _domain access (via Simularion)
            _domain->readXML(xmlconfig);
        }

        if(nullptr != plugin) {
            plugin->readXML(xmlconfig);
            _plugins.push_back(plugin);
        } else {
            global_log->warning() << "Unknown plugin " << pluginname << endl;
        }
    }

    xmlconfig.changecurrentnode(oldpath);

    return numPlugins;
}
