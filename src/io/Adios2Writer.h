#ifndef ADIOS2_WRITER_H_
#define ADIOS2_WRITER_H_
/*
 * \file Adios2Writer.h
 *
 * Allows to write ADIOS2 phase space series files for visualization with Megamol.
 *
 */
#ifdef ENABLE_ADIOS2

#include "molecules/Component.h"
#include "molecules/MoleculeForwardDeclaration.h"
#include "plugins/PluginBase.h"

#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

#include <adios2.h>
#ifdef ENABLE_MPI
#include <mpi.h>
#endif

class Adios2Writer : public PluginBase {
#ifdef MARDYN_DPDP
	using PRECISION = double;
#else
	using PRECISION = float;
#endif

public:
	Adios2Writer() = default;
	~Adios2Writer() override = default;

	void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;

	/** @brief Read in XML configuration for Adios2Writer.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <outputplugin name="Adios2Writer">
		 <outputfile>STRING</outputfile>
		 <adios2enginetype><!-- For possible engines see the ADIOS2 doc (default: BP4) --></adios2enginetype>
		 <writefrequency>INTEGER</writefrequency>
		 <compression><!-- Enables compression. Supported compression libs: SZ and ZFP. ADIOS2 should be compiled with these libraries. (default: none) --> </compression>
		 <compressionaccuracy><!-- Parameter for the SZ compression lib (default: 0.00001) --></compressionaccuracy>
		 <compressionrate><!-- Parameter for the ZFP compression lib (default: 8) --></compressionrate>
		 <appendmode><!-- Enables the append mode to append data to existing files/checkpoints (default: OFF) --></appendmode>
		 <numfiles><!-- Number of Aggregators in ADIOS2 IO object (default: -1) --></numfiles>
	   </outputplugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	void testInit(std::vector<Component>& comps, const std::string outfile = "mardyn.bp",
				  const std::string adios2enginetype = "BP4", const unsigned long writefrequency = 50000,
				  const std::string compression = "none", const std::string compression_accuracy = "0.00001",
				  const std::string compression_rate = "8");

	void beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
								unsigned long simstep) override;

	void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;

	void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) override;

	void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
				 unsigned long simstep) override;

	void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) override;

	std::string getPluginName() override { return std::string("Adios2Writer"); }

	static PluginBase* createInstance() { return new Adios2Writer(); }

protected:
	//
private:
	void resetContainers();
	void clearContainers();
	void defineVariables(const uint64_t global, const uint64_t offset, const uint64_t local, const int numProcs,
						 const int rank);
	void initAdios2();
	std::vector<Component> _comps;
	// output filename, from XML
	std::string _outputfile;
	std::string _adios2enginetype;
	uint32_t _writefrequency;
	std::stringstream _xmlstream;
	std::string _compression;
	std::string _compression_accuracy;
	std::string _compression_rate;
	std::string _append_mode;
	int32_t _num_files;
	// variables to write, see documentation
	std::map<std::string, std::variant<std::vector<double>, std::vector<float>, std::vector<uint64_t>>> _vars;
	// main instance
	std::shared_ptr<adios2::ADIOS> _inst;
	std::shared_ptr<adios2::Engine> _engine;
	std::shared_ptr<adios2::IO> _io;
	adios2::Operator _compressionOperator;
	bool _singleCenter = true;
};
#endif
#endif /* ADIOS2_WRITER_H_*/
