#ifndef ADIOS2_WRITER_H_
#define ADIOS2_WRITER_H_
/*
 * \file Adios2Writer.h
 *
 * Allows to write ADIOS2 phase space series files for visualization with Megamol.
 */

#include "molecules/MoleculeForwardDeclaration.h"
#include "plugins/PluginBase.h"

#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

#include <adios2.h>
#include <mpi.h>


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
		 <adios2enginetype><!-- For possible engines see the ADIOS2 doc --></adios2enginetype>
		 <writefrequency>INTEGER</writefrequency>
	   </outputplugin>
	   \endcode
	 */
	void readXML(XMLfileUnits& xmlconfig) override;

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
	void initAdios2();
	// output filename, from XML
	std::string _outputfile;
	std::string _adios2enginetype;
	uint32_t _writefrequency;
	double current_time;
	std::stringstream _xmlstream;
	// variables to write, see documentation
	std::map<std::string, std::variant<std::vector<double>, std::vector<uint64_t>>> vars;
	// std::map<std::string, std::vector<double>> vars;
	// main instance
	std::shared_ptr<adios2::ADIOS> inst;
	std::shared_ptr<adios2::Engine> engine;
	std::shared_ptr<adios2::IO> io;
};

#endif /* ADIOS2_WRITER_H_*/