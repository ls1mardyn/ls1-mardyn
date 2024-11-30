#ifndef SRC_UTILS_PLUGINFACTORY_H_
#define SRC_UTILS_PLUGINFACTORY_H_

#include <map>
#include <string>
#include <vector>
#include "../Simulation.h"

#include "../utils/Logger.h"
#include "../utils/String_utils.h"

#include "PluginBase.h"



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

	/** @brief Register a plugin
	 *
	 * @param createInstance  pointer to a function returning an instance of the plugin object.
	 */
	void registerDefaultPlugins();

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
		MARDYN_EXIT("Requested plugin not found: " + pluginname + ". Is it enabled in CMake?");
		return nullptr;
	}

	long enablePlugins(std::list<PluginBase*>& _plugins, XMLfileUnits& xmlconfig, std::string category, Domain* _domain);
};
#endif  // SRC_UTILS_PLUGINFACTORY_H_
