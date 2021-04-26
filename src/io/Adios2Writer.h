/*
 * AdiosWriter.h
 *
 *  Created on: 11 Jul 2018
 *      Author: Oliver Fernandes
 */

///
/// \file AdiosWriter.h
/// Insitu Megamol Plugin Header. See the Adios2Writer class description for a manual on how to use the plugin
///

#pragma once

#include "plugins/PluginBase.h"
#include "molecules/MoleculeForwardDeclaration.h"

#include <chrono>
#include <set>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <memory>
#include <map>
#include <errno.h>

#include <variant>

#ifdef ENABLE_ADIOS2

#include <adios2.h>
#include <mpi.h>

class Adios2Writer : public PluginBase {
public:
    Adios2Writer() {};
    virtual ~Adios2Writer() {};

    void init(ParticleContainer* particleContainer,
            DomainDecompBase* domainDecomp, Domain* domain
    );

    void readXML(XMLfileUnits& xmlconfig);

    void beforeEventNewTimestep(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    );

    void beforeForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    );

    void afterForces(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            unsigned long simstep
    );

    void endStep(
            ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
            Domain* domain, unsigned long simstep
    );
    
    void finish(ParticleContainer* particleContainer,
            DomainDecompBase* domainDecomp, Domain* domain);

    std::string getPluginName() {
        return std::string("Adios2Writer");
    }

    static PluginBase* createInstance() { return new Adios2Writer(); }

protected:
    // 
private:
    void initAdios2();
    // output filename, from XML
    std::string fname;
    uint32_t _writefrequency;
    double current_time;
    // variables to write, see documentation
    std::map<std::string, std::variant<std::vector<double>, std::vector<uint64_t>>> vars;
    //std::map<std::string, std::vector<double>> vars;
    // main instance
    std::shared_ptr<adios2::ADIOS> inst;
    std::shared_ptr<adios2::Engine> engine;  
    std::shared_ptr<adios2::IO> io;
};
#endif // ENABLE_ADIOS2
