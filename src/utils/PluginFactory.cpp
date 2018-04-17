//
// Created by Moritz Kruegener on 4/17/2018.
//

#include "PluginFactory.h"
#include "Simulation.h"
#include "Domain.h"

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