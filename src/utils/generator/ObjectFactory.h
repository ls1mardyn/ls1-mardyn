#ifndef SRC_UTILS_GENERATOR_OBJECTFACTORY_H_
#define SRC_UTILS_GENERATOR_OBJECTFACTORY_H_

#include "utils/PluginFactory.h"
#include "utils/generator/Objects.h"

class ObjectFactory : public PluginFactory<Object> {
public:
	ObjectFactory();
};

#endif  // SRC_UTILS_GENERATOR_OBJECTFACTORY_H_
