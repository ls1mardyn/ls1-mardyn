#ifndef SRC_UTILS_PLUGINFACTORY_H_
#define SRC_UTILS_PLUGINFACTORY_H_

#include <map>
#include <string>
#include <vector>
#include <Simulation.h>

#include "utils/Logger.h"
#include "utils/String_utils.h"

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
//#include "io/ResultWriter.h"
#include "io/SysMonOutput.h"
//#include "io/VISWriter.h"
#include "io/XyzWriter.h"
#include "io/MaxWriter.h"

#ifdef VTK
#include "io/vtk/VTKMoleculeWriter.h"
#include "io/vtk/VTKGridWriter.h"
#endif

#include "utils/testPlugin.h"

class Simulation;



/** @brief macro used to simplify the registration of plugins in the constructor.
 * @param NAME class name of the plugin
 */
#define REGISTER_PLUGIN(NAME) registerPlugin(&(NAME::createInstance));


/** @brief Plugin Factory
 *
 * Provides a common interface to register plugins based on the template given as parameter.
 * The interface must include the following method:
 * \code{.cpp}
 *   virtual std::string getPluginName() = 0; // return the name of the plugin
 * \endcode
 * Also each plugin has to implement the following static method:
 * \code{.cpp}
 *   static T* createInstance() { return new MyPlugin(); }  // return an instance object
 * \endcode
 * Plugins are registered in the PluginFactory constructor using the REGISTER_PLUGIN(NAME) macro.
 */
template <typename T>
class PluginFactory {

typedef T* createInstanceFunc();

private:
	std::map<std::string, createInstanceFunc*> _pluginFactoryMap;

public:
	PluginFactory() {}
	~PluginFactory() {}


	/** @brief Register an output plugin
	 *
	 * @param createInstance  pointer to a function returning an instance of the plugin object.
	 */
	void registerPlugin(createInstanceFunc* createInstance) {
		T *pluginInstance = createInstance();
		std::string pluginname = pluginInstance->getPluginName();
		Log::global_log->debug() << "Registering plugin with name " << pluginname << std::endl;
		delete pluginInstance;
		if( _pluginFactoryMap.count(pluginname) > 0 ) {
			Log::global_log->warning() << "Skipping already registered plugin with name " << pluginname << std::endl;
			Log::global_log->debug() << "Registered plugins: " << string_utils::join(getPluginNames(), ", ") << std::endl;
			return;
		}
		_pluginFactoryMap[pluginname] = createInstance;
	}

	/** @brief Register an output plugin
	 *
	 * @param createInstance  pointer to a function returning an instance of the plugin object.
	 */
	void registerDefaultPlugins() {
		global_log -> debug() << "REGISTERING PLUGINS" << endl;
        REGISTER_PLUGIN(testPlugin);
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
		//REGISTER_PLUGIN(ResultWriter);
		REGISTER_PLUGIN(SysMonOutput);
		//REGISTER_PLUGIN(VISWriter);
		REGISTER_PLUGIN(XyzWriter);
		REGISTER_PLUGIN(MaxWriter);
#ifdef VTK
		REGISTER_PLUGIN(VTKMoleculeWriter);
		REGISTER_PLUGIN(VTKGridWriter);
#endif
	}

	/** @brief Get all names of registered plugins */
	std::vector<std::string> getPluginNames() {
		std::vector<std::string> pluginNames;
		for(auto const &plugin : _pluginFactoryMap) {
			pluginNames.push_back(plugin.first);
		}
		return pluginNames;
	}

	/** @brief Create a new instance of plugin */
	T* create(std::string pluginname) {
		auto existing = _pluginFactoryMap.find(pluginname);
		if( existing != _pluginFactoryMap.end() ) {
			return existing->second(); /* call createInstance for plugin */
		}
		Log::global_log->warning() << "Plugin not found: " << pluginname << std::endl;
		return nullptr;
	}

	long enablePlugins(std::list<PluginBase*>& _plugins, XMLfileUnits& xmlconfig, std::string category, Domain* _domain);

	};
#endif  // SRC_UTILS_PLUGINFACTORY_H_