#ifndef OUTPUT_PLUGIN_FACTORY
#define OUTPUT_PLUGIN_FACTORY

#include "utils/PluginFactory.h"

class OutputBase;


/** @brief Output Plugin Factory
 *
 * Provide a common interface to register output plugins based on the OutputBase interface.
 * Therefore each plugin has to implement the following static method:
 * \code{.cpp}
 *   static OutputBase* createInstance() { return new MyPlugin(); }  // returning an instance object
 * \endcode
 * Plugins are registered in the OutputPluginFactory constructor using the REGISTER_PLUGIN(NAME) macro.
 */
class OutputPluginFactory : public PluginFactory<OutputBase>{
public:
	OutputPluginFactory();
};

#endif  // OUTPUT_PLUGIN_FACTORY
