/*
 * CylindricSampling.cpp
 *
 *  Created on: Feb 2022
 *      Author: homes
 */

#include "CylindricSampling.h"

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


CylindricSampling::CylindricSampling() {}

void CylindricSampling::init(ParticleContainer* /* particleContainer */, DomainDecompBase* domainDecomp, Domain* domain) {

    _globalBoxLength[0] = domain->getGlobalLength(0);
    _globalBoxLength[1] = domain->getGlobalLength(1);
    _globalBoxLength[2] = domain->getGlobalLength(2);

    _distMax = 0.5*min(_globalBoxLength[0],_globalBoxLength[2]);  // Only sample particles within maximum radius fitting into box

    _numBinsGlobalHeight = static_cast<unsigned int>(_globalBoxLength[1]/_binwidth);
    _numBinsGlobalRadius = static_cast<unsigned int>(_distMax/_binwidth);

    if (_globalBoxLength[1]/_binwidth != static_cast<float>(_numBinsGlobalHeight)) {
        global_log->error() << "[CylindricSampling] Can not divide domain without remainder in y-direction! Change binwidth" << std::endl;
        Simulation::exit(-1);
    }
    if (_distMax/_binwidth != static_cast<float>(_numBinsGlobalRadius)) {
        global_log->error() << "[CylindricSampling] Can not divide domain without remainder in x or z-direction! Change binwidth" << std::endl;
        Simulation::exit(-1);
    }

     // Entry per component and bin; 0 represents all components combined
    _lenVector = _numBinsGlobalRadius * _numBinsGlobalHeight;

    resizeVectors();
    resetVectors();
}

void CylindricSampling::readXML(XMLfileUnits& xmlconfig) {

    xmlconfig.getNodeValue("binwidth", _binwidth);
    xmlconfig.getNodeValue("start", _startSampling);
    xmlconfig.getNodeValue("writefrequency", _writeFrequency);
    xmlconfig.getNodeValue("stop", _stopSampling);


    global_log->info() << "[CylindricSampling] Start:WriteFreq:Stop: " << _startSampling << " : " << _writeFrequency << " : " << _stopSampling << std::endl;
    global_log->info() << "[CylindricSampling] Binwidth: " << _binwidth << std::endl;
    global_log->info() << "[CylindricSampling] All components treated as single one" << std::endl;
}

void CylindricSampling::afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) {

    // Sampling starts after _startSampling and is conducted up to _stopSampling
    if ((simstep <= _startSampling) or (simstep > _stopSampling)) {
        return;
    }

    // Do not write or sample data directly after (re)start in the first step
    if ( simstep == _simulation.getNumInitTimesteps() ) {
        return;
    }

    // Variables per step
    std::array<double, 3> regionLowCorner {};
    std::array<double, 3> regionHighCorner {};
    std::array<double, 3> regionSize {};
    for (unsigned short d = 0; d < 3; d++) {
        regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
        regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
        regionSize[d] = regionHighCorner[d] - regionLowCorner[d];
    }


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
        ekinVect_step.at(d).local.resize(_lenVector);
        velocityVect_step.at(d).local.resize(_lenVector);
        virialVect_step.at(d).local.resize(_lenVector);

        ekinVect_step.at(d).global.resize(_lenVector);
        velocityVect_step.at(d).global.resize(_lenVector);
        virialVect_step.at(d).global.resize(_lenVector);
    }

    std::fill(numMolecules_step.local.begin(), numMolecules_step.local.end(), 0ul);
    std::fill(mass_step.local.begin(),        mass_step.local.end(), 0.0f);
    std::fill(ekin_step.local.begin(),        ekin_step.local.end(), 0.0f);

    std::fill(numMolecules_step.global.begin(), numMolecules_step.global.end(), 0ul);
    std::fill(mass_step.global.begin(),        mass_step.global.end(), 0.0f);
    std::fill(ekin_step.global.begin(),        ekin_step.global.end(), 0.0f);

    for (unsigned short d = 0; d < 3; d++) {
        std::fill(ekinVect_step.at(d).local.begin(),      ekinVect_step.at(d).local.end(), 0.0f);
        std::fill(velocityVect_step.at(d).local.begin(),   velocityVect_step.at(d).local.end(), 0.0f);
        std::fill(virialVect_step.at(d).local.begin(),     virialVect_step.at(d).local.end(), 0.0f);

        std::fill(ekinVect_step.at(d).global.begin(),      ekinVect_step.at(d).global.end(), 0.0f);
        std::fill(velocityVect_step.at(d).global.begin(),   velocityVect_step.at(d).global.end(), 0.0f);
        std::fill(virialVect_step.at(d).global.begin(),     virialVect_step.at(d).global.end(), 0.0f);
    }

    // Calculate drift as it is needed first
    for (auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit) {
        const double ry = pit->r(1);
        const double distCenter_x = pit->r(0)-0.5*_globalBoxLength[0];
        const double distCenter_z = pit->r(2)-0.5*_globalBoxLength[2];
        const double distCenter = std::sqrt(std::pow(distCenter_x,2) + std::pow(distCenter_z,2));
        // Do not consider particles outside of most outer radius
        if (distCenter >= _distMax) { continue; }
        const unsigned int indexH = std::min(_numBinsGlobalHeight, static_cast<unsigned int>(ry/_binwidth));  // Index of bin of height
        const unsigned int indexR = std::min(_numBinsGlobalRadius, static_cast<unsigned int>(distCenter/_binwidth));  // Index of bin of radius
        const unsigned int index = _numBinsGlobalHeight*indexR + indexH;

        numMolecules_step.local.at(index) ++;

        const double u_x = pit->v(0);
        const double u_y = pit->v(1);
        const double u_z = pit->v(2);
        const double mass = pit->mass();
        const double vi_x = pit->Vi(0);
        const double vi_y = pit->Vi(1);
        const double vi_z = pit->Vi(2);

        // Transform to cylindric coordinates
        const double velo_r = (distCenter_x*u_x + distCenter_z*u_z)/distCenter;  // Radial
        const double velo_y = u_y;
        const double velo_t = (distCenter_x*u_z - distCenter_z*u_x)/distCenter;  // Tangential

        mass_step.local.at(index) += mass;
        ekin_step.local.at(index) += pit->U_kin();

        velocityVect_step.at(0).local.at(index) += velo_r;
        velocityVect_step.at(1).local.at(index) += velo_y;
        velocityVect_step.at(2).local.at(index) += velo_t;

        ekinVect_step[0].local.at(index) += 0.5*mass*velo_r*velo_r;
        ekinVect_step[1].local.at(index) += 0.5*mass*velo_y*velo_y;
        ekinVect_step[2].local.at(index) += 0.5*mass*velo_t*velo_t;

        // See comment below why the virial is not transformed to cylindric coordinates
        // short: must be done during force calculation for each particle pair individually
        virialVect_step[0].local.at(index) += vi_x;
        virialVect_step[1].local.at(index) += vi_y;
        virialVect_step[2].local.at(index) += vi_z;
    }

// Gather quantities needed by all processes
#ifdef ENABLE_MPI
    MPI_Allreduce(numMolecules_step.local.data(), numMolecules_step.global.data(), _lenVector, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(velocityVect_step[0].local.data(), velocityVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(velocityVect_step[1].local.data(), velocityVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(velocityVect_step[2].local.data(), velocityVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    for (unsigned long i = 0; i < _lenVector; i++) {
        numMolecules_step.global.at(i) = numMolecules_step.local.at(i);
        velocityVect_step[0].global.at(i) = velocityVect_step[0].local.at(i);
        velocityVect_step[1].global.at(i) = velocityVect_step[1].local.at(i);
        velocityVect_step[2].global.at(i) = velocityVect_step[2].local.at(i);
    }
#endif

// Gather other quantities. Note: MPI_Reduce instead of MPI_Allreduce!
#ifdef ENABLE_MPI
    MPI_Reduce(mass_step.local.data(), mass_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekin_step.local.data(), ekin_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[0].local.data(), virialVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[1].local.data(), virialVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[2].local.data(), virialVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekinVect_step[0].local.data(), ekinVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekinVect_step[1].local.data(), ekinVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekinVect_step[2].local.data(), ekinVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    for (unsigned long i = 0; i < _lenVector; i++) {
        mass_step.global.at(i) = mass_step.local.at(i);
        ekin_step.global.at(i) = ekin_step.local.at(i);
        virialVect_step[0].global.at(i) = virialVect_step[0].local.at(i);
        virialVect_step[1].global.at(i) = virialVect_step[1].local.at(i);
        virialVect_step[2].global.at(i) = virialVect_step[2].local.at(i);
        ekinVect_step[0].global.at(i) = ekinVect_step[0].local.at(i);
        ekinVect_step[1].global.at(i) = ekinVect_step[1].local.at(i);
        ekinVect_step[2].global.at(i) = ekinVect_step[2].local.at(i);
    }
#endif

    // Only root knows real quantities (MPI_Reduce instead of MPI_Allreduce)
    // Accumulate data
    if (domainDecomp->getRank() == 0) {
        for (unsigned long i = 0; i < _lenVector; i++) {
            const unsigned long numMols = numMolecules_step.global.at(i);
            unsigned int dof_rot {0};
            unsigned int dof_total {0};

            // For single component sampling, the rot. DOF of component 0 is taken
            dof_rot = _simulation.getEnsemble()->getComponent(0)->getRotationalDegreesOfFreedom();
            dof_total = (3 + dof_rot)*numMols;

            const double vi_x = virialVect_step[0].global.at(i);
            const double vi_y = virialVect_step[1].global.at(i);
            const double vi_z = virialVect_step[2].global.at(i);

            _doftotal_accum.at(i)                += dof_total;
            _numMolecules_accum.at(i)            += numMols;
            _mass_accum.at(i)                    += mass_step.global.at(i);
            _ekin_accum.at(i)                    += ekin_step.global.at(i);
            _virial_accum.at(i)                  += vi_x + vi_y + vi_z;

            _ekinVect_accum[0].at(i)             += ekinVect_step[0].global.at(i);
            _ekinVect_accum[1].at(i)             += ekinVect_step[1].global.at(i);
            _ekinVect_accum[2].at(i)             += ekinVect_step[2].global.at(i);
            _virialVect_accum[0].at(i)           += vi_x;
            _virialVect_accum[1].at(i)           += vi_y;
            _virialVect_accum[2].at(i)           += vi_z;

            if (numMols > 0ul) {
                _velocityVect_accum[0].at(i)         += velocityVect_step[0].global.at(i) / numMols;
                _velocityVect_accum[1].at(i)         += velocityVect_step[1].global.at(i) / numMols;
                _velocityVect_accum[2].at(i)         += velocityVect_step[2].global.at(i) / numMols;
            }

            _countSamples.at(i)++;
        }
    }

    // Write out data every _writeFrequency step
    if ( (simstep - _startSampling) % _writeFrequency == 0 ) {

        if (domainDecomp->getRank() == 0) {
            // Write output file
            std::stringstream ss;
            ss << std::setw(9) << std::setfill('0') << simstep;
            const std::string fname = "CylindricSampling_TS"+ss.str()+".dat";
            std::ofstream ofs;
            ofs.open(fname, std::ios::out);
            ofs << setw(24) << "height"     // Bin position (height)
                << setw(24) << "radius"     // Bin position (radius)
                << setw(24) << "numParts"   // Average number of molecules in bin per step
                << setw(24) << "rho"        // Density
                << setw(24) << "T"          // Temperature without drift (i.e. "real" temperature)
                << setw(24) << "ekin"       // Kinetic energy including drift
                << setw(24) << "p"          // Pressure
                << setw(24) << "T_r"        // Temperature in radial direction
                << setw(24) << "T_y"        // Temperature in y-direction
                << setw(24) << "T_t"        // Temperature in tangantial direction
                << setw(24) << "v_r"        // Drift velocity in radial direction
                << setw(24) << "v_y"        // Drift velocity in y-direction
                << setw(24) << "v_t"        // Drift velocity in tangantial direction
                << setw(24) << "p_r"        // Pressure in radial direction; the radial pressure is not easily accessible (see comment below)
                << setw(24) << "p_y"        // Pressure in y-direction
                << setw(24) << "p_t"        // Pressure in tangantial direction; the tangential pressure is not easily accessible (see comment below)
                << setw(24) << "numSamples";    // Number of samples (<= _writeFrequency)
            ofs << std::endl;

            for (unsigned long i = 0; i < _lenVector; i++) {
                const unsigned int idxH = i % _numBinsGlobalHeight;  // Remainder
                const unsigned int idxR = i / _numBinsGlobalHeight; // Quotient
                ofs << FORMAT_SCI_MAX_DIGITS << (idxH+0.5)*_binwidth;  // Height bin
                ofs << FORMAT_SCI_MAX_DIGITS << (idxR+0.5)*_binwidth;  // Radius bin
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
                if ((_countSamples.at(i) > 0ul) and (_doftotal_accum.at(i)) > 0ul) {
                    const double numMols_accum = static_cast<double>(_numMolecules_accum.at(i));
                    const double slabVolume = 3.1415926536*_binwidth*(std::pow((idxR+1)*_binwidth,2)-std::pow(idxR*_binwidth,2));
                    numSamples = _countSamples.at(i);

                    numMolsPerStep = numMols_accum/numSamples;
                    rho         = numMolsPerStep           / slabVolume;
                    v_r         = _velocityVect_accum[0].at(i)   / numSamples;
                    v_y         = _velocityVect_accum[1].at(i)   / numSamples;
                    v_t         = _velocityVect_accum[2].at(i)   / numSamples;

                    double v_drift_sqr = v_r*v_r + v_y*v_y + v_t*v_t;

                    T           = (2*_ekin_accum.at(i) - v_drift_sqr*_mass_accum.at(i)) / _doftotal_accum.at(i);
                    ekin        = _ekin_accum.at(i) / numMols_accum;
                    p           = rho * ( (_virial_accum.at(i))/(3.0*numMols_accum) + T);

                    T_r         = (2*_ekinVect_accum[0].at(i) - (v_r*v_r)*_mass_accum.at(i)) / numMols_accum;
                    T_y         = (2*_ekinVect_accum[1].at(i) - (v_y*v_y)*_mass_accum.at(i)) / numMols_accum;
                    T_t         = (2*_ekinVect_accum[2].at(i) - (v_t*v_t)*_mass_accum.at(i)) / numMols_accum;
                    // The radial and tangential virial/pressure are not easily accessible
                    // It must be computed/transformed during the force calculation for each particle pair individually
                    // see e.g. the spherical LRC
                    // p_r         = rho * ( _virialVect_accum[0].at(i)/numMols_accum + T);
                    p_y         = rho * ( _virialVect_accum[1].at(i)/numMols_accum + T);
                    // p_t         = rho * ( _virialVect_accum[2].at(i)/numMols_accum + T);
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
void CylindricSampling::resizeVectors() {

    _numMolecules_accum.resize(_lenVector);
    _doftotal_accum.resize(_lenVector);
    _mass_accum.resize(_lenVector);
    _ekin_accum.resize(_lenVector);
    _virial_accum.resize(_lenVector);

    for (unsigned short d = 0; d < 3; d++) {
        _ekinVect_accum.at(d).resize(_lenVector);
        _velocityVect_accum.at(d).resize(_lenVector);
        _virialVect_accum.at(d).resize(_lenVector);
    }

    _countSamples.resize(_lenVector);
}

// Fill vectors with zeros
void CylindricSampling::resetVectors() {
    std::fill(_numMolecules_accum.begin(), _numMolecules_accum.end(), 0ul);
    std::fill(_doftotal_accum.begin(), _doftotal_accum.end(), 0ul);
    std::fill(_mass_accum.begin(), _mass_accum.end(), 0.0f);
    std::fill(_ekin_accum.begin(), _ekin_accum.end(), 0.0f);
    std::fill(_virial_accum.begin(), _virial_accum.end(), 0.0f);

    for (unsigned short d = 0; d < 3; d++) {
        std::fill(_ekinVect_accum.at(d).begin(), _ekinVect_accum.at(d).end(), 0.0f);
        std::fill(_velocityVect_accum.at(d).begin(), _velocityVect_accum.at(d).end(), 0.0f);
        std::fill(_virialVect_accum.at(d).begin(), _virialVect_accum.at(d).end(), 0.0f);
    }

    std::fill(_countSamples.begin(), _countSamples.end(), 0ul);
}

