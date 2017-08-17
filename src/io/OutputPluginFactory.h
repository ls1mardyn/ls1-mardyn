#ifndef OUTPUT_PLUGIN_FACTORY
#define OUTPUT_PLUGIN_FACTORY

#include <map>
#include <string>
#include <vector>

class OutputBase;


typedef OutputBase* createInstanceFunc();

/** @brief Output Plugin Factory
 *
 * Provide a comman interface to register output plugins based on the OutputBase interface.
 * Therefore each plugin has to implement the following static method:
 * \code{.cpp}
 *   static OutputBase* createInstance() { return new MyPlugin(); }  // returning an instance object
 * \endcode
 * Plugins are registered in the OutputPluginFactory constructor using the REGISTER_PLUGIN(NAME) macro.
 */
class OutputPluginFactory {
public:
	OutputPluginFactory();
	~OutputPluginFactory(){}

	/** @brief Register an output plugin
	 *
	 * @param createInstance  pointer to a function returning an instance of the plugin object.
	 */
	void registerPlugin(createInstanceFunc* createInstance);

	/** @brief Get all names of registered plugins */
	std::vector<std::string> getPluginNames();

	/** @brief Create a new instance of plugin */
	OutputBase* create(std::string pluginname);

private:
	std::map<std::string, createInstanceFunc*> _outputPluginFactoryMap;
};

#endif  // OUTPUT_PLUGIN_FACTORY
