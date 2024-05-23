/*
 * CapiSampling.cpp
 *
 *  Created on: May 2024
 *      Author: jniemann
 */

#include "CapiSampling.h"

#include <math.h>
#include <limits>
#include <algorithm>
#include <map>

#include "particleContainer/ParticleContainer.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "utils/Random.h"
#include "utils/FileUtils.h"
#include "Simulation.h"
#include "longRange/LongRangeCorrection.h"


CapiSampling::CapiSampling() {}

void CapiSampling::init(ParticleContainer* /* particleContainer */, DomainDecompBase* domainDecomp, Domain* domain) {

    _globalBoxLength[0] = domain->getGlobalLength(0);
    _globalBoxLength[1] = domain->getGlobalLength(1);
    _globalBoxLength[2] = domain->getGlobalLength(2);

    _yLiqLower = _globalBoxLength[1]/2 - _minimumFilmWidth/2;
    _yLiqUpper = _globalBoxLength[1]/2 + _minimumFilmWidth/2;
    _yGasLeft = _minimumGasPhase;
    _yGasRight = _globalBoxLength[1] - _minimumGasPhase;

    //check for consistency of bulk phase coordinates:
    if ((_yGasLeft >= 0) or (_yGasLeft >=_yLiqLower) or (_yLiqLower>= _yLiqUpper) or (_yLiqUpper >= _yGasRight))
    {
        Log::global_log->error() << "["<< getPluginName()<<"] Bulk Phase Coordinates are invalid! " << std::endl;
        Log::global_log->error() << "["<< getPluginName()<<"] Coords are: _yGasLeft:_yLiqLower:_yLiqUpper:_yGasRight==" 
                                 <<  _yGasLeft<<":"<<_yLiqLower<<":"<<_yLiqUpper<<":"<<_yGasRight<< std::endl;
        Simulation::exit(-1);
    }
    

    // check if numBins is valid
    if ((_numBinsX < 1) or (_numBinsY <1) or (_numBinsZ <1)) {
        Log::global_log->error() << "["<< getPluginName()<<"] Please select valid vlaues for number of Bins in X and Z direction" << std::endl;
        Simulation::exit(-1);
    }
    _binwidthX = (_globalBoxLength[0]/_numBinsX);
    _binwidthY = (_globalBoxLength[1]/_numBinsY);
    _binwidthZ = (_globalBoxLength[2]/_numBinsZ);
    _cellVolume = _binwidthX * _binwidthZ;

    Log::global_log->info() << "["<< getPluginName()<<"] _binwidthX:_binwidthZ : " <<  _binwidthX<<":"<<_binwidthZ << std::endl;

    // Entry per bin; all components sampled as one
    _lenVector = _numBinsX * _numBinsY * _numBinsZ;

    resizeVectors();
    resetVectors();
}

void CapiSampling::readXML(XMLfileUnits& xmlconfig) {

    xmlconfig.getNodeValue("numBinsX", _numBinsX);
    xmlconfig.getNodeValue("numBinsZ", _numBinsZ);
    xmlconfig.getNodeValue("writefrequency", _writeFrequency);
    xmlconfig.getNodeValue("start", _startSampling);
    xmlconfig.getNodeValue("stop", _stopSampling);
    xmlconfig.getNodeValue("minimumFilmWidth", _minimumFilmWidth);
    xmlconfig.getNodeValue("minimumGasPhase", _minimumGasPhase);

    Log::global_log->info() << "["<< getPluginName()<<"] Start:WriteFreq:Stop: " << _startSampling << " : " << _writeFrequency << " : " << _stopSampling << std::endl;
    Log::global_log->info() << "["<< getPluginName()<<"] All components treated as single one" << std::endl;
}

void CapiSampling::afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) {

    // Sampling starts after _startSampling and is conducted up to _stopSampling
    if ((simstep < _startSampling) or (simstep > _stopSampling)) {
        return;
    }

    CommVar<std::vector<unsigned long>> numMolecules_step;

    numMolecules_step.local.resize(_lenVector);
    numMolecules_step.global.resize(_lenVector);

    std::fill(numMolecules_step.local.begin(), numMolecules_step.local.end(), 0ul);
    std::fill(numMolecules_step.global.begin(), numMolecules_step.global.end(), 0ul);

    for (auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit) {
        const double rx = pit->r(0);
        const double ry = pit->r(1);
        const double rz = pit->r(2);
        const unsigned int indexX = std::min(_numBinsX, static_cast<unsigned int>(rx/_binwidthX));  // Index of bin of height
        const unsigned int indexY = std::min(_numBinsY, static_cast<unsigned int>(ry/_binwidthY));  // Index of bin of height
        const unsigned int indexZ = std::min(_numBinsZ, static_cast<unsigned int>(rz/_binwidthZ));  // Index of bin of radius
        const unsigned int index =  _numBinsY*_numBinsX*indexY  + _numBinsX*indexX + indexZ;

        numMolecules_step.local[index] ++;

        // counting for rhol & rhov:
        unsigned long numMoleculesLiq_step = 0ul;
        unsigned long numMoleculesGas_step = 0ul;
        if (ry < _yGasLeft){ numMoleculesGas_step++;}
        else if (ry < _yLiqLower){ continue;}
        else if (ry < _yLiqUpper){ numMoleculesLiq_step++;}
        else if (ry < _yGasRight){ continue;}
        else { numMoleculesGas_step++;}
    }
    // #TODO: continue with runtime-calculation of rho_liq / rho_gas

// Gather quantities. Note: MPI_Reduce instead of MPI_Allreduce! Therefore, only root has correct values
#ifdef ENABLE_MPI
    MPI_Reduce(numMolecules_step.local.data(), numMolecules_step.global.data(), _lenVector, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    for (unsigned long i = 0; i < _lenVector; i++) {
        numMolecules_step.global[i] = numMolecules_step.local[i];
    }
#endif

    // Only root knows real quantities (MPI_Reduce instead of MPI_Allreduce)
    // Accumulate data
    if (domainDecomp->getRank() == 0) {
        for (unsigned long i = 0; i < _lenVector; i++) {
            _numMolecules_accum[i]  += numMolecules_step.global[i];
            _countSamples[i]++;
        }
    }

    // Write out data every _writeFrequency step
    if ( (simstep - _startSampling) % _writeFrequency == 0 ) {

        if (domainDecomp->getRank() == 0) {
            // Write output file
            std::stringstream ss;
            ss << std::setw(9) << std::setfill('0') << simstep;
            const std::string fname = "CapiSampling_TS"+ss.str()+".dat";
            std::ofstream ofs;
            ofs.open(fname, std::ios::out);

            // head
            ofs << std::setw(24) << "yIndex";   // Bin Index (y) is written in colum 
            ofs << std::setw(24) << "x";        // Bin position (x) is written in colum
            for(int iz = 0; iz < _numBinsZ; iz++ ){
                double rz = (.5 + iz) * _binwidthZ;
                ofs << std::setw(24) << rz;        // Bin position (z) is written in row
            }

            //content
            for(int ix = 0; ix<_numBinsX; ix++){
                double rx = (.5 + ix) * _binwidthX;
                for(int iy = 0; iy<_numBinsY; iy++){
                    ofs << std::setw(24) << iy;        // Bin Index (y)
                    ofs << std::setw(24) << rx;        // Bin position (x)
                    for (int iz = 0; iz < _numBinsZ; iz++)
                    {
                        int index = _numBinsY*_numBinsX*iy  + _numBinsX*ix + iz;

                        unsigned long numSamples {0ul};
                        double numMolsPerStep {std::nan("0")}; // Not an int as particles change bin during simulation
                        double rho {0.0};
                        if (_countSamples[index] > 0ul) {
                            const double numMols_accum = static_cast<double>(_numMolecules_accum[index]);
                            numSamples = _countSamples[index];

                            numMolsPerStep = numMols_accum/numSamples;
                            rho         = numMolsPerStep / _cellVolume;
                        }
                        ofs << std::setw(24) << rho;        
                    }
                    ofs << std::endl;
                }
            }
            ofs.close();

        // Reset vectors to zero
        resetVectors();
        }
    }
}

// Resize vectors
void CapiSampling::resizeVectors() {

    _numMolecules_accum.resize(_lenVector);
    _countSamples.resize(_lenVector);
}

// Fill vectors with zeros
void CapiSampling::resetVectors() {
    std::fill(_numMolecules_accum.begin(), _numMolecules_accum.end(), 0ul);
    _numMolecules_liq_accum = 0;  
    _numMolecules_vap_accum = 0;
    
    std::fill(_countSamples.begin(), _countSamples.end(), 0ul);
    _countSamples_liq = 0;
    _countSamples_vap = 0;
}
