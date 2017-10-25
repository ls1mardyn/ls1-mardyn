#include "utils/generator/ObjectFillerFactory.h"

#include "utils/generator/GridFiller.h"

ObjectFillerFactory::ObjectFillerFactory(){
	REGISTER_PLUGIN(GridFiller);
}
