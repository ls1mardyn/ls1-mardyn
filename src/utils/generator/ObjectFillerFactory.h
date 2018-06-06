#ifndef OBJECTFILLER_FACTORY
#define OBJECTFILLER_FACTORY

#include "plugins/PluginFactory.h"

class ObjectFillerBase;


/** @brief Filler Factory
 *
 * Provide a common interface to register Filler based on the FillerBase interface.
 * Therefore each plugin has to implement the following static method:
 * \code{.cpp}
 *   static OutputBase* createInstance() { return new MyPlugin(); }  // returning an instance object
 * \endcode
 * Plugins are registered in the OutputPluginFactory constructor using the REGISTER_PLUGIN(NAME) macro.
 */
class ObjectFillerFactory : public PluginFactory<ObjectFillerBase>{
public:
	   ObjectFillerFactory();
};

#endif  // OBJECTFILLER_FACTORY
