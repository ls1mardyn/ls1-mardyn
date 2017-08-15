#ifndef OUTPUT_PLUGIN_FACTORY
#define OUTPUT_PLUGIN_FACTORY

#include <map>
#include <string>
#include <vector>

#include "io/OutputBase.h"


typedef OutputBase* createInstanceFunc();

/** @brief Output Plugin Factory
 *
 * Provide a comman interface to register output plugins based on the OutputBase interface.
 * Therefore each plugin has to implement the two static methods:
 * \code{.cpp}
 *   static OutputBase* createInstance(); // returning an instance object
 *   static std::string getPluginName();  // returning the name of the plugin as used in the XML file
 * \endcode
 * Plugins are registered in the OutputPluginFactory constructor using the REGISTER_PLUGIN(NAME) macro.
 */
class OutputPluginFactory {
public:
	OutputPluginFactory();
	~OutputPluginFactory(){}

	/** @brief Register an output plugin
	 *
	 * @par pluginname      name of the plugin
	 * @par createInstance  pointer to a function returning an instance of the plugin object.
	 */
	void registerPlugin(std::string pluginname, createInstanceFunc* createInstance);

	/** @brief Get all names of registered plugins */
	std::vector<string> getPluginNames();

	/** @brief Create a new instance of plugin */
	OutputBase* create(std::string pluginname);

private:
	std::map<std::string, createInstanceFunc*> _outputPluginFactoryMap;
};

#endif  // OUTPUT_PLUGIN_FACTORY
