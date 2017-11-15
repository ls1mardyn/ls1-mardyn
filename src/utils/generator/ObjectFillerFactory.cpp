#include "utils/generator/ObjectFillerFactory.h"

#include "utils/generator/GridFiller.h"
#include "utils/generator/ReplicaFiller.h"

ObjectFillerFactory::ObjectFillerFactory(){
	REGISTER_PLUGIN(GridFiller);
	REGISTER_PLUGIN(ReplicaFiller);
}
