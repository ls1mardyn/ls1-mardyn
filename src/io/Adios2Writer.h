#ifndef ADIOS2_WRITER_H_
#define ADIOS2_WRITER_H_
/*
 * \file Adios2Writer.h
 *
 * Allows to write ADIOS2 phase space series files for visualization with Megamol.
 *
 */
#ifdef ENABLE_ADIOS2
#include "molecules/MoleculeForwardDeclaration.h"
#include "plugins/PluginBase.h"
#include "molecules/Component.h"

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
public:
	Adios2Writer(){};
	virtual ~Adios2Writer(){};

	void init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);

	/** @brief Read in XML configuration for Adios2Writer.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
	   <outputplugin name="Adios2Writer">
		 <outputfile>STRING</outputfile>
		 <adios2enginetype><!-- For possible engines see the ADIOS2 doc (default: BP4) --></adios2enginetype>
		 <writefrequency>INTEGER</writefrequency>
	   </outputplugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

	void testInit(std::vector<Component>& comps, const std::string outfile = "mardyn.bp", const std::string adios2enginetype = "BP4", const unsigned long writefrequency = 50000,
				  const std::string compression = "none", const std::string compression_accuracy = "0.00001",
				  const std::string compression_rate = "8");

	void beforeEventNewTimestep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
								unsigned long simstep);

	void beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep);

	void afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep);

	void endStep(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
				 unsigned long simstep);

	void finish(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);

	std::string getPluginName() { return std::string("Adios2Writer"); }

	static PluginBase* createInstance() { return new Adios2Writer(); }

protected:
	//
private:
	void resetContainers();
	void clearContainers();
	void defineVariables(uint64_t global, uint64_t offset, uint64_t local, int numProcs, int rank);
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
	// variables to write, see documentation
	std::map<std::string, std::variant<std::vector<double>, std::vector<uint64_t>>> _vars;
	// std::map<std::string, std::vector<double>> vars;
	// main instance
	std::shared_ptr<adios2::ADIOS> _inst;
	std::shared_ptr<adios2::Engine> _engine;
	std::shared_ptr<adios2::IO> _io;
	adios2::Operator _compressionOperator;
	bool _singleCenter = true;
};
#endif
#endif /* ADIOS2_WRITER_H_*/