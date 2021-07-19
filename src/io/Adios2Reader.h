#ifndef ADIOS2_READER_H_
#define ADIOS2_READER_H_
/*
 * \file Adios2Reader.h
 *
 * Allows to read ADIOS2 phase space series files as checkpoints.
 *
 */
#ifdef ENABLE_ADIOS2

#include "io/InputBase.h"
#include "molecules/MoleculeForwardDeclaration.h"
#include "utils/Logger.h"

#include <chrono>
#include <set>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <memory>
#include <map>
#include <errno.h>

#include "adios2.h"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif

class Adios2Reader : public InputBase {
public:
    Adios2Reader() {}; 
    virtual ~Adios2Reader() {};

    void init(ParticleContainer* particleContainer,
            DomainDecompBase* domainDecomp, Domain* domain
    );

	/** @brief Read in XML configuration for Adios2Reader.
	 *
	 * The following xml object structure is handled by this method:
	 * \code{.xml}
		<phasespacepoint>
			<file type="adios2">
				<filename>STRING</filename>
				<adios2enginetype>STRING<!-- For possible engines see the ADIOS2 doc (default: BP4) --></adios2enginetype>
				<adios2step>INTEGER<!-- Step in the ADIOS2 file --></adios2step>
				<mode>STRING<!-- Possible options are rootOnly and equalRanks(Not implemented yet) --></mode>
			</file>
		</phasespacepoint>
	   \endcode
	*/
	void readXML(XMLfileUnits& xmlconfig);

	void readPhaseSpaceHeader(Domain* domain, double timestep);

    std::string getPluginName() {
        return std::string("Adios2Reader");
    }

   unsigned long readPhaseSpace(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp);


protected:
    // 
private:
    void initAdios2();
    void rootOnlyRead(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp);
    void equalRanksRead(ParticleContainer* particleContainer, Domain* domain, DomainDecompBase* domainDecomp) {};
    // output filename, from XML
    std::string _inputfile;
    std::string _adios2enginetype;
    std::string _mode;
    int _step;
    uint64_t particle_count;
    // main instance
    std::shared_ptr<adios2::ADIOS> inst;
    std::shared_ptr<adios2::Engine> engine;  
    std::shared_ptr<adios2::IO> io;

    template <typename T> 
    void doTheRead(std::string var_name, std::vector<T>& container) {
      uint64_t num = 1;
      auto advar = io->InquireVariable<T>(var_name);
      advar.SetStepSelection({_step,1});
      auto info = engine->BlocksInfo(advar, 0);
      auto shape = info[0].Count;
      std::for_each(shape.begin(), shape.end(), [&](decltype(num) n) { num *= n; });
      container.resize(num);
      engine->Get<T>(advar, container);
    }

    template <typename T> 
    void doTheRead(std::string var_name, std::vector<T>& container, uint64_t buffer, uint64_t offset) {

      adios2::Dims readsize({buffer});
      adios2::Dims offs({offset});

      Log::global_log->info() << "[Adios2Reader]: Var name " << var_name << std::endl;
      auto advar = io->InquireVariable<T>(var_name);
      advar.SetStepSelection({_step,1});
      Log::global_log->info() << "[Adios2Reader]: buffer " << buffer << " offset " << offset << std::endl;
      for (auto entry: advar.Shape())
        Log::global_log->info() << "[Adios2Reader]: shape " << entry << std::endl;
      advar.SetSelection(adios2::Box<adios2::Dims>(offs, readsize));
      container.resize(buffer);
      engine->Get<T>(advar, container);
    }

    template <typename T> 
    void doTheRead(std::string var_name, T& value) {
      auto advar = io->InquireVariable<T>(var_name);
      advar.SetStepSelection({_step,1});
      engine->Get<T>(advar, value);
    }

};
#endif // ENABLE_ADIOS2
#endif /* ADIOS2_READER_H_*/