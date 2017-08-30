#include "io/OutputPluginFactory.h"

#include "io/OutputBase.h"
#include "utils/Logger.h"

#include "io/CavityWriter.h"
#include "io/CheckpointWriter.h"
#include "io/DecompWriter.h"
#include "io/EnergyLogWriter.h"
/** @todo fix Interface missmatch */
// #include "io/FlopRateWriter.h"
#include "io/GammaWriter.h"
#include "io/LoadBalanceWriter.h"
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


OutputPluginFactory::OutputPluginFactory() {
	REGISTER_PLUGIN(CavityWriter);
	REGISTER_PLUGIN(CheckpointWriter);
	REGISTER_PLUGIN(DecompWriter);
	REGISTER_PLUGIN(EnergyLogWriter);
/** @todo fix Interface missmatch */
// 	REGISTER_PLUGIN(FlopRateWriter);
	REGISTER_PLUGIN(GammaWriter);
	REGISTER_PLUGIN(LoadbalanceWriter);
	REGISTER_PLUGIN(MPICheckpointWriter);
	REGISTER_PLUGIN(MmpldWriter);
	REGISTER_PLUGIN(MmspdBinWriter);
	REGISTER_PLUGIN(MmspdWriter);
	REGISTER_PLUGIN(PovWriter);
	REGISTER_PLUGIN(RDF);
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
