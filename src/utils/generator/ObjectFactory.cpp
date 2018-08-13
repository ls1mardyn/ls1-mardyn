#include "utils/generator/ObjectFactory.h"

#include "utils/generator/Objects.h"

ObjectFactory::ObjectFactory() {
	REGISTER_PLUGIN(Cuboid);
	REGISTER_PLUGIN(Sphere);
	REGISTER_PLUGIN(Cylinder);
	REGISTER_PLUGIN(ObjectUnification);
	REGISTER_PLUGIN(ObjectSubtractor);
	REGISTER_PLUGIN(ObjectIntersection);
	REGISTER_PLUGIN(ObjectShifter);
}
