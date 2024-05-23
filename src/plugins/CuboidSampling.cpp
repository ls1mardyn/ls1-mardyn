/*
 * CuboidSampling.cpp
 *
 *  Created on: Feb 2022
 *      Author: homes
 */

#include "CuboidSampling.h"

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


CuboidSampling::CuboidSampling() {}

void CuboidSampling::init(ParticleContainer* /* particleContainer */, DomainDecompBase* domainDecomp, Domain* domain) {

    _globalBoxLength[0] = domain->getGlobalLength(0);
    _globalBoxLength[1] = domain->getGlobalLength(1);
    _globalBoxLength[2] = domain->getGlobalLength(2);

    if (_numBinsX < 1 or _numBinsZ <1) {
        Log::global_log->error() << "["<< getPluginName()<<"] Please select valid vlaues for number of Bins in X and Z direction" << std::endl;
        Simulation::exit(-1);
    }
    _binwidthX = (_globalBoxLength[0]/_numBinsX);
    _binwidthZ = (_globalBoxLength[2]/_numBinsZ);
    _cellVolume = _binwidthX * _binwidthZ;
    Log::global_log->info() << "["<< getPluginName()<<"] _binwidthX:_binwidthZ : " <<  _binwidthX<<":"<<_binwidthZ << std::endl;

    // Entry per bin; all components sampled as one
    _lenVector = _numBinsX * _numBinsZ;

    resizeVectors();
    resetVectors();
}

void CuboidSampling::readXML(XMLfileUnits& xmlconfig) {

    xmlconfig.getNodeValue("numBinsX", _numBinsX);
    xmlconfig.getNodeValue("numBinsZ", _numBinsZ);
    xmlconfig.getNodeValue("start", _startSampling);
    xmlconfig.getNodeValue("writefrequency", _writeFrequency);
    xmlconfig.getNodeValue("stop", _stopSampling);
    xmlconfig.getNodeValue("yLower", _yLower);
    xmlconfig.getNodeValue("yUpper", _yUpper);


    Log::global_log->info() << "["<< getPluginName()<<"] Start:WriteFreq:Stop: " << _startSampling << " : " << _writeFrequency << " : " << _stopSampling << std::endl;
    Log::global_log->info() << "["<< getPluginName()<<"] All components treated as single one" << std::endl;
}

void CuboidSampling::afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) {

    // Sampling starts after _startSampling and is conducted up to _stopSampling
    if ((simstep < _startSampling) or (simstep > _stopSampling)) {
        return;
    }

    // Do not write or sample data directly after (re)start in the first step
    // if ( simstep == _simulation.getNumInitTimesteps() ) {
    //     return;
    // }

    // Variables per step
    // std::array<double, 3> regionLowCorner {};
    // std::array<double, 3> regionHighCorner {};
    // std::array<double, 3> regionSize {};
    // for (unsigned short d = 0; d < 3; d++) {
    //     regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
    //     regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
    //     regionSize[d] = regionHighCorner[d] - regionLowCorner[d];
    // }


    CommVar<std::vector<unsigned long>> numMolecules_step;
    CommVar<std::vector<double>> mass_step;
    CommVar<std::vector<double>> ekin_step;                        // Including drift energy
    std::array<CommVar<std::vector<double>>, 3> ekinVect_step;
    std::array<CommVar<std::vector<double>>, 3> velocityVect_step;
    std::array<CommVar<std::vector<double>>, 3> virialVect_step;

    numMolecules_step.local.resize(_lenVector);
    mass_step.local.resize(_lenVector);
    ekin_step.local.resize(_lenVector);

    numMolecules_step.global.resize(_lenVector);
    mass_step.global.resize(_lenVector);
    ekin_step.global.resize(_lenVector);

    for (unsigned short d = 0; d < 3; d++) {
        ekinVect_step[d].local.resize(_lenVector);
        velocityVect_step[d].local.resize(_lenVector);
        virialVect_step[d].local.resize(_lenVector);

        ekinVect_step[d].global.resize(_lenVector);
        velocityVect_step[d].global.resize(_lenVector);
        virialVect_step[d].global.resize(_lenVector);
    }

    std::fill(numMolecules_step.local.begin(), numMolecules_step.local.end(), 0ul);
    std::fill(mass_step.local.begin(),        mass_step.local.end(), 0.0f);
    std::fill(ekin_step.local.begin(),        ekin_step.local.end(), 0.0f);

    std::fill(numMolecules_step.global.begin(), numMolecules_step.global.end(), 0ul);
    std::fill(mass_step.global.begin(),        mass_step.global.end(), 0.0f);
    std::fill(ekin_step.global.begin(),        ekin_step.global.end(), 0.0f);

    for (unsigned short d = 0; d < 3; d++) {
        std::fill(ekinVect_step[d].local.begin(),      ekinVect_step[d].local.end(), 0.0f);
        std::fill(velocityVect_step[d].local.begin(),   velocityVect_step[d].local.end(), 0.0f);
        std::fill(virialVect_step[d].local.begin(),     virialVect_step[d].local.end(), 0.0f);

        std::fill(ekinVect_step[d].global.begin(),      ekinVect_step[d].global.end(), 0.0f);
        std::fill(velocityVect_step[d].global.begin(),   velocityVect_step[d].global.end(), 0.0f);
        std::fill(virialVect_step[d].global.begin(),     virialVect_step[d].global.end(), 0.0f);
    }

    for (auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit) {
        const double rx = pit->r(0);
        const double ry = pit->r(1);
        const double rz = pit->r(2);
        // Do not consider particles outside of most outer radius
        if ((ry < _yLower) or (ry > _yUpper)) { continue; }
        const unsigned int indexX = std::min(_numBinsX, static_cast<unsigned int>(rx/_binwidthX));  // Index of bin of height
        const unsigned int indexZ = std::min(_numBinsZ, static_cast<unsigned int>(rz/_binwidthZ));  // Index of bin of radius
        const unsigned int index = _numBinsX*indexX + indexZ;

        numMolecules_step.local[index] ++;

        const double u_x = pit->v(0);
        const double u_y = pit->v(1);
        const double u_z = pit->v(2);
        const double mass = pit->mass();
        const double vi_x = pit->Vi(0);
        const double vi_y = pit->Vi(1);
        const double vi_z = pit->Vi(2);

        mass_step.local[index] += mass;
        ekin_step.local[index] += pit->U_kin();

        velocityVect_step[0].local[index] += u_x;
        velocityVect_step[1].local[index] += u_y;
        velocityVect_step[2].local[index] += u_z;

        ekinVect_step[0].local[index] += 0.5*mass*u_x*u_x;
        ekinVect_step[1].local[index] += 0.5*mass*u_y*u_y;
        ekinVect_step[2].local[index] += 0.5*mass*u_z*u_z;

        virialVect_step[0].local[index] += vi_x;
        virialVect_step[1].local[index] += vi_y;
        virialVect_step[2].local[index] += vi_z;
    }

// Gather quantities. Note: MPI_Reduce instead of MPI_Allreduce! Therefore, only root has correct values
#ifdef ENABLE_MPI
    MPI_Reduce(numMolecules_step.local.data(), numMolecules_step.global.data(), _lenVector, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(mass_step.local.data(), mass_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekin_step.local.data(), ekin_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(velocityVect_step[0].local.data(), velocityVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(velocityVect_step[1].local.data(), velocityVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(velocityVect_step[2].local.data(), velocityVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[0].local.data(), virialVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[1].local.data(), virialVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[2].local.data(), virialVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekinVect_step[0].local.data(), ekinVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekinVect_step[1].local.data(), ekinVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekinVect_step[2].local.data(), ekinVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    for (unsigned long i = 0; i < _lenVector; i++) {
        numMolecules_step.global[i] = numMolecules_step.local[i];
        mass_step.global[i] = mass_step.local[i];
        ekin_step.global[i] = ekin_step.local[i];
        velocityVect_step[0].global[i] = velocityVect_step[0].local[i];
        velocityVect_step[1].global[i] = velocityVect_step[1].local[i];
        velocityVect_step[2].global[i] = velocityVect_step[2].local[i];
        virialVect_step[0].global[i] = virialVect_step[0].local[i];
        virialVect_step[1].global[i] = virialVect_step[1].local[i];
        virialVect_step[2].global[i] = virialVect_step[2].local[i];
        ekinVect_step[0].global[i] = ekinVect_step[0].local[i];
        ekinVect_step[1].global[i] = ekinVect_step[1].local[i];
        ekinVect_step[2].global[i] = ekinVect_step[2].local[i];
    }
#endif

    // Only root knows real quantities (MPI_Reduce instead of MPI_Allreduce)
    // Accumulate data
    if (domainDecomp->getRank() == 0) {
        for (unsigned long i = 0; i < _lenVector; i++) {
            const unsigned long numMols = numMolecules_step.global[i];
            unsigned int dof_rot {0};
            unsigned int dof_total {0};

            // For single component sampling, the rot. DOF of component 0 is taken
            dof_rot = _simulation.getEnsemble()->getComponent(0)->getRotationalDegreesOfFreedom();
            dof_total = (3 + dof_rot)*numMols;

            const double vi_x = virialVect_step[0].global[i];
            const double vi_y = virialVect_step[1].global[i];
            const double vi_z = virialVect_step[2].global[i];

            _doftotal_accum[i]                += dof_total;
            _numMolecules_accum[i]            += numMols;
            _mass_accum[i]                    += mass_step.global[i];
            _ekin_accum[i]                    += ekin_step.global[i];
            _virial_accum[i]                  += vi_x + vi_y + vi_z;

            _ekinVect_accum[0][i]             += ekinVect_step[0].global[i];
            _ekinVect_accum[1][i]             += ekinVect_step[1].global[i];
            _ekinVect_accum[2][i]             += ekinVect_step[2].global[i];
            _virialVect_accum[0][i]           += vi_x;
            _virialVect_accum[1][i]           += vi_y;
            _virialVect_accum[2][i]           += vi_z;

            if (numMols > 0ul) {
                _velocityVect_accum[0][i]         += velocityVect_step[0].global[i] / numMols;
                _velocityVect_accum[1][i]         += velocityVect_step[1].global[i] / numMols;
                _velocityVect_accum[2][i]         += velocityVect_step[2].global[i] / numMols;
            }

            _countSamples[i]++;
        }
    }

    // Write out data every _writeFrequency step
    if ( (simstep - _startSampling) % _writeFrequency == 0 ) {

        if (domainDecomp->getRank() == 0) {
            // Write output file
            std::stringstream ss;
            ss << std::setw(9) << std::setfill('0') << simstep;
            const std::string fname = "CuboidSampling_TS"+ss.str()+".dat";
            std::ofstream ofs;
            ofs.open(fname, std::ios::out);
            ofs << std::setw(24) << "rx"     // Bin position (x)
                << std::setw(24) << "rz"     // Bin position (z)
                << std::setw(24) << "numParts"   // Average number of molecules in bin per step
                << std::setw(24) << "rho"        // Density
                << std::setw(24) << "T"          // Temperature without drift (i.e. "real" temperature)
                << std::setw(24) << "ekin"       // Kinetic energy including drift
                << std::setw(24) << "p"          // Pressure
                << std::setw(24) << "T_r"        // Temperature in radial direction
                << std::setw(24) << "T_y"        // Temperature in y-direction
                << std::setw(24) << "T_t"        // Temperature in tangantial direction
                << std::setw(24) << "v_r"        // Drift velocity in radial direction
                << std::setw(24) << "v_y"        // Drift velocity in y-direction
                << std::setw(24) << "v_t"        // Drift velocity in tangantial direction
                << std::setw(24) << "p_r"        // Pressure in radial direction; the radial pressure is not easily accessible (see comment below)
                << std::setw(24) << "p_y"        // Pressure in y-direction
                << std::setw(24) << "p_t"        // Pressure in tangantial direction; the tangential pressure is not easily accessible (see comment below)
                << std::setw(24) << "numSamples";    // Number of samples (<= _writeFrequency)
            ofs << std::endl;

            for (unsigned long i = 0; i < _lenVector; i++) {
                const unsigned int idxX = i % _numBinsX;  // Remainder
                const unsigned int idxZ = i / _numBinsX; // Quotient
                ofs << FORMAT_SCI_MAX_DIGITS << (idxX+0.5)*_binwidthX;  // x bin
                ofs << FORMAT_SCI_MAX_DIGITS << (idxZ+0.5)*_binwidthZ;  // z bin
                unsigned long numSamples {0ul};
                double numMolsPerStep {std::nan("0")}; // Not an int as particles change bin during simulation
                double rho {0.0};
                double T {std::nan("0")};
                double ekin {std::nan("0")};
                double p {std::nan("0")};
                double T_r {std::nan("0")};
                double T_y {std::nan("0")};
                double T_t {std::nan("0")};
                double v_r {std::nan("0")};
                double v_y {std::nan("0")};
                double v_t {std::nan("0")};
                double p_r {std::nan("0")};
                double p_y {std::nan("0")};
                double p_t {std::nan("0")};
                if ((_countSamples[i] > 0ul) and (_doftotal_accum[i]) > 0ul) {
                    const double numMols_accum = static_cast<double>(_numMolecules_accum[i]);
                    numSamples = _countSamples[i];

                    numMolsPerStep = numMols_accum/numSamples;
                    rho         = numMolsPerStep           / _cellVolume;
                    v_r         = _velocityVect_accum[0][i]   / numSamples;
                    v_y         = _velocityVect_accum[1][i]   / numSamples;
                    v_t         = _velocityVect_accum[2][i]   / numSamples;

                    double v_drift_sqr = v_r*v_r + v_y*v_y + v_t*v_t;

                    T           = (2*_ekin_accum[i] - v_drift_sqr*_mass_accum[i]) / _doftotal_accum[i];
                    ekin        = _ekin_accum[i] / numMols_accum;
                    p           = rho * ( (_virial_accum[i])/(3.0*numMols_accum) + T);

                    T_r         = (2*_ekinVect_accum[0][i] - (v_r*v_r)*_mass_accum[i]) / numMols_accum;
                    T_y         = (2*_ekinVect_accum[1][i] - (v_y*v_y)*_mass_accum[i]) / numMols_accum;
                    T_t         = (2*_ekinVect_accum[2][i] - (v_t*v_t)*_mass_accum[i]) / numMols_accum;
                    // The radial and tangential virial/pressure are not easily accessible
                    // It must be computed/transformed during the force calculation for each particle pair individually
                    // see e.g. the spherical LRC
                    // p_r         = rho * ( _virialVect_accum[0][i]/numMols_accum + T);
                    p_y         = rho * ( _virialVect_accum[1][i]/numMols_accum + T);
                    // p_t         = rho * ( _virialVect_accum[2][i]/numMols_accum + T);
                }
                ofs << FORMAT_SCI_MAX_DIGITS << numMolsPerStep
                    << FORMAT_SCI_MAX_DIGITS << rho
                    << FORMAT_SCI_MAX_DIGITS << T
                    << FORMAT_SCI_MAX_DIGITS << ekin
                    << FORMAT_SCI_MAX_DIGITS << p
                    << FORMAT_SCI_MAX_DIGITS << T_r
                    << FORMAT_SCI_MAX_DIGITS << T_y
                    << FORMAT_SCI_MAX_DIGITS << T_t
                    << FORMAT_SCI_MAX_DIGITS << v_r
                    << FORMAT_SCI_MAX_DIGITS << v_y
                    << FORMAT_SCI_MAX_DIGITS << v_t
                    << FORMAT_SCI_MAX_DIGITS << p_r
                    << FORMAT_SCI_MAX_DIGITS << p_y
                    << FORMAT_SCI_MAX_DIGITS << p_t
                    << FORMAT_SCI_MAX_DIGITS << numSamples
                    << std::endl;
            }
            ofs.close();
        }

        // Reset vectors to zero
        resetVectors();
    }
}

// Resize vectors
void CuboidSampling::resizeVectors() {

    _numMolecules_accum.resize(_lenVector);
    _doftotal_accum.resize(_lenVector);
    _mass_accum.resize(_lenVector);
    _ekin_accum.resize(_lenVector);
    _virial_accum.resize(_lenVector);

    for (unsigned short d = 0; d < 3; d++) {
        _ekinVect_accum[d].resize(_lenVector);
        _velocityVect_accum[d].resize(_lenVector);
        _virialVect_accum[d].resize(_lenVector);
    }

    _countSamples.resize(_lenVector);
}

// Fill vectors with zeros
void CuboidSampling::resetVectors() {
    std::fill(_numMolecules_accum.begin(), _numMolecules_accum.end(), 0ul);
    std::fill(_doftotal_accum.begin(), _doftotal_accum.end(), 0ul);
    std::fill(_mass_accum.begin(), _mass_accum.end(), 0.0f);
    std::fill(_ekin_accum.begin(), _ekin_accum.end(), 0.0f);
    std::fill(_virial_accum.begin(), _virial_accum.end(), 0.0f);

    for (unsigned short d = 0; d < 3; d++) {
        std::fill(_ekinVect_accum[d].begin(), _ekinVect_accum[d].end(), 0.0f);
        std::fill(_velocityVect_accum[d].begin(), _velocityVect_accum[d].end(), 0.0f);
        std::fill(_virialVect_accum[d].begin(), _virialVect_accum[d].end(), 0.0f);
    }
    
    std::fill(_countSamples.begin(), _countSamples.end(), 0ul);
}
