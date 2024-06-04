/*
 * CoalescenceSampling.cpp
 *
 *  Created on: May 2024
 *      Author: niemann
 */

#include "CoalescenceSampling.h"

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



CoalescenceSampling::CoalescenceSampling() {}

void CoalescenceSampling::init(ParticleContainer* /* particleContainer */, DomainDecompBase* domainDecomp, Domain* domain) {
    Log::global_log->info() << "["<< getPluginName()<<"] ----------------INIT CALLED: "  << std::endl;
    _rhoOutputFilename = _outputPrefix +  "rho_r"+ std::to_string(_radius)+".csv";

    _globalBoxLength[0] = domain->getGlobalLength(0);
    _globalBoxLength[1] = domain->getGlobalLength(1);
    _globalBoxLength[2] = domain->getGlobalLength(2);

    _yLower = std::max(_yLower, 0.f);
    _yUpper = std::min(_yUpper, (float) _globalBoxLength[1]);
    _radius = std::min(_radius, 0.5f *  (float) std::min(_globalBoxLength[0],_globalBoxLength[2]));

    _numBinsGlobal = static_cast<unsigned int>((_yUpper - _yLower)/_binwidth);

    if ((_yUpper - _yLower)/_binwidth != static_cast<float>(_numBinsGlobal)) {
        Log::global_log->info() << "["<< getPluginName()<<"] Can not divide domain without remainder in y-direction! " << std::endl;
        
        float yCenter = (_yUpper + _yLower)/2;
        float adjustedWidth = ((float) _numBinsGlobal) * _binwidth;
        _yLower = yCenter - adjustedWidth/2.f;
        _yUpper = yCenter + adjustedWidth/2.f;
                
        Log::global_log->info() << "["<< getPluginName()<<"] Changing (_yLower, _yUpper) to ("<<_yLower<<", "<<_yUpper<<")." << std::endl;
    }


    // prepare rho Outputfile 
    if (domainDecomp->getRank() == 0) {
		std::ofstream ofs(_rhoOutputFilename, std::ios::out);
		ofs << std::setw(24) << "simstep";
        for (unsigned long i = 0; i < _numBinsGlobal; i++) {
            ofs << std::setw(24) << _yLower+(i+0.5)*_binwidth;   
        }
		ofs << std::endl;
		ofs.close();
	}

    resizeVectors();
    resetVectors();
}



void CoalescenceSampling::readXML(XMLfileUnits& xmlconfig) {

    Log::global_log->info() << "["<< getPluginName()<<"] ----------------READXML CALLED: "  << std::endl;
    xmlconfig.getNodeValue("binwidth", _binwidth);
    xmlconfig.getNodeValue("start", _startSampling);
    xmlconfig.getNodeValue("writefrequency", _writeFrequency);
    xmlconfig.getNodeValue("stop", _stopSampling);
    xmlconfig.getNodeValue("yLower", _yLower);
    xmlconfig.getNodeValue("yUpper", _yUpper);
    xmlconfig.getNodeValue("radius", _radius);
    xmlconfig.getNodeValue("outputPrefix", _outputPrefix);


    Log::global_log->info() << "["<< getPluginName()<<"] Start:WriteFreq:Stop: " << _startSampling << " : " << _writeFrequency << " : " << _stopSampling << std::endl;
    Log::global_log->info() << "["<< getPluginName()<<"] Binwidth: " << _binwidth << std::endl;
    Log::global_log->info() << "["<< getPluginName()<<"] All components treated as single one" << std::endl;
}        



void CoalescenceSampling::afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) {

    // Sampling starts after _startSampling and is conducted up to _stopSampling
    if ((simstep < _startSampling) or (simstep > _stopSampling)) {
        return;
    }

    // Do not write or sample data directly after (re)start in the first step
    // if ( simstep == _simulation.getNumInitTimesteps() ) {
    //     return;
    // }

    CommVar<std::vector<unsigned long>> numMolecules_step;
    CommVar<std::vector<double>> mass_step;
    CommVar<std::vector<double>> ekin_step;                        // Including drift energy
    std::array<CommVar<std::vector<double>>, 3> ekinVect_step;
    std::array<CommVar<std::vector<double>>, 3> velocityVect_step;
    std::array<CommVar<std::vector<double>>, 3> virialVect_step;

    numMolecules_step.local.resize(_numBinsGlobal);
    mass_step.local.resize(_numBinsGlobal);

    numMolecules_step.global.resize(_numBinsGlobal);
    mass_step.global.resize(_numBinsGlobal);

    std::fill(numMolecules_step.local.begin(), numMolecules_step.local.end(), 0ul);

    std::fill(numMolecules_step.global.begin(), numMolecules_step.global.end(), 0ul);

    for (auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit) {
        // Do not consider particles outside of [_yLower,_yUpper]
        const double ry = pit->r(1);
        if (ry<_yLower or ry>_yUpper) { continue; }
        // Do not consider particles outside of most outer radius
        const double distCenter_x = pit->r(0)-0.5*_globalBoxLength[0];
        const double distCenter_z = pit->r(2)-0.5*_globalBoxLength[2];
        const double distCenter = std::sqrt(std::pow(distCenter_x,2) + std::pow(distCenter_z,2));
        if (distCenter >= _radius) { continue; }
        const unsigned int index = std::min(_numBinsGlobal, static_cast<unsigned int>((ry-_yLower)/_binwidth));  // Index of bin of height

        numMolecules_step.local[index] ++;
    }

// Gather quantities. Note: MPI_Reduce instead of MPI_Allreduce! Therefore, only root has correct values
#ifdef ENABLE_MPI
    MPI_Reduce(numMolecules_step.local.data(), numMolecules_step.global.data(), _numBinsGlobal, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    for (unsigned long i = 0; i < _numBinsGlobal; i++) {
        numMolecules_step.global[i] = numMolecules_step.local[i];
    }
#endif

    // Only root knows real quantities (MPI_Reduce instead of MPI_Allreduce)
    // Accumulate data
    if (domainDecomp->getRank() == 0) {
        for (unsigned long i = 0; i < _numBinsGlobal; i++) {
            const unsigned long numMols = numMolecules_step.global[i];
            _numMolecules_accum[i] += numMols;
            _countSamples[i]++;
        }
    }

    // Write out data every _writeFrequency step
    if ( (simstep - _startSampling) % _writeFrequency == 0 ) {
        // Log::global_log->info() << "------------------------------------------A_-------------------";
        // write rho Output 
        if (domainDecomp->getRank() == 0) {
            std::ofstream ofs(_rhoOutputFilename, std::ios::app);
            ofs << std::setw(24) << simstep;

            for (unsigned long i = 0; i < _numBinsGlobal; i++) {
                unsigned long numSamples {0ul};
                double numMolsPerStep {std::nan("0")}; // Not an int as particles change bin during simulation
                double rho {0.0};
                if (_countSamples[i] > 0ul){
                    const double numMols_accum = static_cast<double>(_numMolecules_accum[i]);
                    const double cellVolume = _binwidth*3.1415926536*std::pow(_radius,2);
                    numSamples = _countSamples[i];

                    numMolsPerStep = numMols_accum/numSamples;
                    rho         = numMolsPerStep  / cellVolume;
                }
                ofs << std::setw(24) << rho;   
            }
            ofs << std::endl;
            ofs.close();
        }

        // Reset vectors to zero
        resetVectors();
    }
}




// Resize vectors
void CoalescenceSampling::resizeVectors() {

    _numMolecules_accum.resize(_numBinsGlobal);
    _countSamples.resize(_numBinsGlobal);
}

// Fill vectors with zeros
void CoalescenceSampling::resetVectors() {
    std::fill(_numMolecules_accum.begin(), _numMolecules_accum.end(), 0ul);
    std::fill(_countSamples.begin(), _countSamples.end(), 0ul);
}
