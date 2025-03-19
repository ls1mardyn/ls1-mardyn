#include "io/CavityWriter.h"

#include <fstream>
#include <sstream>

#include "Common.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "utils/mardyn_assert.h"
#include "Simulation.h"
#include "ensemble/CavityEnsemble.h"

void CavityWriter::readXML(XMLfileUnits &xmlconfig) {
    _writeFrequency = 1;
    xmlconfig.getNodeValue("writefrequency", _writeFrequency);
    Log::global_log->info() << "[CavityWriter] Write frequency: " << _writeFrequency << std::endl;

    _outputPrefix = "mardyn";
    xmlconfig.getNodeValue("outputprefix", _outputPrefix);
    Log::global_log->info() << "[CavityWriter] Output prefix: " << _outputPrefix << std::endl;

    xmlconfig.getNodeValue("incremental", _incremental);
    Log::global_log->info() << "[CavityWriter] Incremental numbers: " << _incremental << std::endl;

    xmlconfig.getNodeValue("appendTimestamp", _appendTimestamp);
    Log::global_log->info() << "[CavityWriter] Append timestamp: " << _appendTimestamp << std::endl;


    xmlconfig.getNodeValue("maxNeighbours", _maxNeighbors);
    if (_maxNeighbors <= 0) {
        std::ostringstream error_message;
        error_message << "[CavityWriter] Invalid number of maxNeighbors: " << _maxNeighbors << std::endl;
        MARDYN_EXIT(error_message.str());
    }

    xmlconfig.getNodeValue("radius", _radius);
    if (_radius <= 0.0f) {
        std::ostringstream error_message;
        error_message << "[CavityWriter] Invalid size of radius: " << _radius << std::endl;
        MARDYN_EXIT(error_message.str());
    }
    xmlconfig.getNodeValue("Nx", _Nx);
    if (_Nx <= 0) {
        std::ostringstream error_message;
        error_message << "[CavityWriter] Invalid number of cells Nx: " << _Nx << std::endl;
        MARDYN_EXIT(error_message.str());
    }
    xmlconfig.getNodeValue("Ny", _Ny);
    if (_Ny <= 0) {
        std::ostringstream error_message;
        error_message << "[CavityWriter] Invalid number of cells Ny: " << _Ny << std::endl;
        MARDYN_EXIT(error_message.str());
    }
    xmlconfig.getNodeValue("Nz", _Nz);
    if (_Nz <= 0) {
        std::ostringstream error_message;
        error_message << "[CavityWriter] Invalid number of cells Nz: " << _Nz << std::endl;
        MARDYN_EXIT(error_message.str());
    }

    // Default Control Volume is entire Domain
    for (int d = 0; d < 3; d++) {
        _controlVolume[d * 2] = 0;
        _controlVolume[d * 2 + 1] = global_simulation->getDomain()->getGlobalLength(d);
    }
    xmlconfig.getNodeValue("ControlVolume/x0", _controlVolume[0]);
    xmlconfig.getNodeValue("ControlVolume/x1", _controlVolume[1]);
    xmlconfig.getNodeValue("ControlVolume/y0", _controlVolume[2]);
    xmlconfig.getNodeValue("ControlVolume/y1", _controlVolume[3]);
    xmlconfig.getNodeValue("ControlVolume/z0", _controlVolume[4]);
    xmlconfig.getNodeValue("ControlVolume/z1", _controlVolume[5]);
    for (int d = 0; d < 3; d++) {
        if (_controlVolume[d * 2] > _controlVolume[d * 2 + 1]) {
            std::ostringstream error_message;
            error_message << "[CavityWriter] Lower Bound of Control Volume may not be larger than upper bound." << std::endl;
            MARDYN_EXIT(error_message.str());
        }
        if (_controlVolume[d * 2] < 0 ||
            _controlVolume[d * 2 + 1] > global_simulation->getDomain()->getGlobalLength(d)) {
            std::ostringstream error_message;
            error_message << "[CavityWriter] Control volume bounds may not be outside of domain boundaries." << std::endl;
            MARDYN_EXIT(error_message.str());
        }
    }

    // Get components to check
    // get root
    std::string oldpath = xmlconfig.getcurrentnodepath();

    this->_mcav = std::map<unsigned, CavityEnsemble *>();

    // iterate over all components
    XMLfile::Query query = xmlconfig.query("componentid");
    for (auto pluginIter = query.begin(); pluginIter; ++pluginIter) {
        xmlconfig.changecurrentnode(pluginIter);
        int componentID = -1;
        xmlconfig.getNodeValue(xmlconfig.getcurrentnodepath(), componentID);
        Log::global_log->info() << "[CavityWriter] Component: " << componentID << std::endl;
        CavityEnsemble *cav = new CavityEnsemble();
        _mcav[componentID] = cav;
    }

    // go back to root
    xmlconfig.changecurrentnode(oldpath);
}

void CavityWriter::init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
                        Domain *domain) {

    // set initial values for CavityEnsemble

    double Tcur = domain->getGlobalCurrentTemperature();

    /* FIXME: target temperature from thermostat ID 0 or 1? */
    double Ttar = domain->severalThermostats() ? domain->getTargetTemperature(1)
                                               : domain->getTargetTemperature(0);
    if ((Tcur < 0.85 * Ttar) || (Tcur > 1.15 * Ttar))
        Tcur = Ttar;

    std::map<unsigned, CavityEnsemble *>::iterator ceit;
    for (ceit = _mcav.begin(); ceit != _mcav.end(); ceit++) {

        // setup
        ceit->second->setSystem(domain->getGlobalLength(0), domain->getGlobalLength(1), domain->getGlobalLength(2),
                                _maxNeighbors, _radius);
        int rank = domainDecomp->getRank();
        double min[3], max[3];
        domainDecomp->getBoundingBoxMinMax(domain, min, max);
        ceit->second->setControlVolume(_controlVolume[0], _controlVolume[2], _controlVolume[4], _controlVolume[1],
                                       _controlVolume[3], _controlVolume[5]);
        ceit->second->setSubdomain(rank, min[0], max[0], min[1], max[1], min[2], max[2]);
        ceit->second->submitTemperature(Tcur);
        int cID = ceit->first;
        Component *c = global_simulation->getEnsemble()->getComponent(cID);
        Log::global_log->info() << "[Cavity Writer] init: " << cID << std::endl;
        ceit->second->init(c, _Nx, _Ny, _Nz);
        Log::global_log->info() << "[Cavity Writer] init done: " << cID << std::endl;
    }

}

void CavityWriter::beforeEventNewTimestep(
        ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
        unsigned long simstep
) {
    if (simstep >= global_simulation->getInitStatistics() && simstep % _writeFrequency == 0) {
        std::map<unsigned, CavityEnsemble *>::iterator ceit;
        for (ceit = this->_mcav.begin(); ceit != this->_mcav.end(); ceit++) {
            ceit->second->preprocessStep();
        }
    }
}

void CavityWriter::afterForces(
        ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
        unsigned long simstep
) {

    if (simstep >= global_simulation->getInitStatistics() && simstep % _writeFrequency == 0) {
        std::map<unsigned, CavityEnsemble *>::iterator ceit;
        for (ceit = this->_mcav.begin(); ceit != this->_mcav.end(); ceit++) {

            ceit->second->cavityStep(particleContainer);
            ceit->second->communicateNumCavities(domainDecomp);

        }
    }
}

void CavityWriter::endStep(ParticleContainer * /*particleContainer*/, DomainDecompBase *domainDecomp,
                           Domain * /*domain*/, unsigned long simstep) {

    if (simstep % _writeFrequency == 0) {
        std::map<unsigned, CavityEnsemble *>::iterator ceit;

        std::map<unsigned, std::stringstream *> cav_filenamestream;
        for (ceit = _mcav.begin(); ceit != _mcav.end(); ceit++) {
            cav_filenamestream[ceit->first] = new std::stringstream;
            *cav_filenamestream[ceit->first] << _outputPrefix << "-c" << ceit->first;
        }

        if (_incremental) {
            unsigned long numTimesteps = _simulation.getNumTimesteps();
            int num_digits = (int) ceil(log(double(numTimesteps / _writeFrequency)) / log(10.));
            for (ceit = _mcav.begin(); ceit != _mcav.end(); ceit++) {
                *cav_filenamestream[ceit->first] << "-" << aligned_number(simstep / _writeFrequency, num_digits, '0');
            }
        }

        for (ceit = _mcav.begin(); ceit != _mcav.end(); ceit++) {
            *cav_filenamestream[ceit->first] << ".cav.xyz";
            Log::global_log->info() << "[CavityWriter] outputName: " << cav_filenamestream[ceit->first]->str() << std::endl;
        }

        int ownRank = domainDecomp->getRank();
        if (ownRank == 0) {
            for (ceit = _mcav.begin(); ceit != _mcav.end(); ceit++) {
                std::ofstream cavfilestream(cav_filenamestream[ceit->first]->str().c_str());
                cavfilestream << ceit->second->numCavities() << std::endl;
                cavfilestream << "comment line" << std::endl;
                cavfilestream.close();
            }
        }

        for (int process = 0; process < domainDecomp->getNumProcs(); process++) {
            domainDecomp->barrier();
            if (ownRank == process) {
                for (ceit = _mcav.begin(); ceit != _mcav.end(); ceit++) {

                    std::ofstream cavfilestream(cav_filenamestream[ceit->first]->str().c_str(), std::ios::app);

                    std::map<unsigned long, Molecule *> tcav = ceit->second->activeParticleContainer();
                    std::map<unsigned long, Molecule *>::iterator tcit;
                    for (tcit = tcav.begin(); tcit != tcav.end(); tcit++) {
                        //global_log->info() << "[CavityWriter] output6" << std::endl;

                        if (ceit->first == 0) { cavfilestream << "C "; }
                        else if (ceit->first == 1) { cavfilestream << "N "; }
                        else if (ceit->first == 2) { cavfilestream << "O "; }
                        else if (ceit->first == 3) { cavfilestream << "F "; }
                        else { cavfilestream << "Ne "; }
                        cavfilestream << tcit->second->r(0) << "\t" << tcit->second->r(1) << "\t" << tcit->second->r(2)
                                      << "\n";
                    }

                    cavfilestream.close();
                }
            }
        }
    }
}

void CavityWriter::finish(ParticleContainer * /*particleContainer*/, DomainDecompBase * /*domainDecomp*/,
                          Domain * /*domain*/ ) {}
