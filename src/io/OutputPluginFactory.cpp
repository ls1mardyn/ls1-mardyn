#include "io/OutputPluginFactory.h"

#include "io/OutputBase.h"
#include "utils/Logger.h"

#include "io/CavityWriter.h"
#include "io/CheckpointWriter.h"
#include "io/DecompWriter.h"
/** @todo fix Interface missmatch */
// #include "io/FlopRateWriter.h"
#include "io/GammaWriter.h"
#include "io/MPICheckpointWriter.h"
#include "io/MmpldWriter.h"
#include "io/MmspdBinWriter.h"
#include "io/MmspdWriter.h"
#include "io/PovWriter.h"
#include "io/RDF.h"
#include "io/ResultWriter.h"
#include "io/SysMonOutput.h"
#include "io/VISWriter.h"
#include "io/XyzWriter.h"
#include "io/MaxWriter.h"

#ifdef VTK
#include "io/vtk/VTKMoleculeWriter.h"
#include "io/vtk/VTKGridWriter.h"
#endif

#define REGISTER_PLUGIN(NAME) registerPlugin(&(NAME::createInstance));

OutputPluginFactory::OutputPluginFactory():
_outputPluginFactoryMap(){
	REGISTER_PLUGIN(CavityWriter);
	REGISTER_PLUGIN(CheckpointWriter);
	REGISTER_PLUGIN(DecompWriter);
/** @todo fix Interface missmatch */
// 	REGISTER_PLUGIN(FlopRateWriter);
	REGISTER_PLUGIN(GammaWriter);
	REGISTER_PLUGIN(MPICheckpointWriter);
	REGISTER_PLUGIN(MmpldWriter);
	REGISTER_PLUGIN(MmspdBinWriter);
	REGISTER_PLUGIN(MmspdWriter);
	REGISTER_PLUGIN(PovWriter);
//	REGISTER_PLUGIN(RDF);
	REGISTER_PLUGIN(ResultWriter);
	REGISTER_PLUGIN(SysMonOutput);
	REGISTER_PLUGIN(VISWriter);
	REGISTER_PLUGIN(XyzWriter);
	REGISTER_PLUGIN(MaxWriter);
#ifdef VTK
	REGISTER_PLUGIN(VTKMoleculeWriter);
	REGISTER_PLUGIN(VTKGridWriter);
#endif /* VTK */
}

#include "utils/String_utils.h"

void OutputPluginFactory::registerPlugin(createInstanceFunc* createInstance) {
	OutputBase *pluginInstance = createInstance();
	std::string pluginname = pluginInstance->getPluginName();
	Log::global_log->info() << "Registering plugin with name " << pluginname << std::endl;
	delete pluginInstance;
	if( _outputPluginFactoryMap.count(pluginname) > 0 ) {
		Log::global_log->warning() << "Skipping already registered plugin with name " << pluginname << std::endl;
		Log::global_log->debug() << "Registered plugins: " << string_utils::join(getPluginNames(),", ") << std::endl;
		return;
	}
	_outputPluginFactoryMap[pluginname] = createInstance;
}

std::vector<string> OutputPluginFactory::getPluginNames() {
	std::vector<string> pluginNames;
	for(auto pluginIter : _outputPluginFactoryMap) {
		pluginNames.push_back(pluginIter.first);
	}
	return pluginNames;
}

OutputBase * OutputPluginFactory::create(std::string pluginname) {
	auto existing = _outputPluginFactoryMap.find(pluginname);
	if( existing != _outputPluginFactoryMap.end() ) {
		return existing->second(); /* call createInstance for plugin */
	}
	Log::global_log->warning() << "Plugin not found: " << pluginname << std::endl;
	return nullptr;
}
