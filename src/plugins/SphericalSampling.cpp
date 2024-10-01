/*
 * SphericalSampling.cpp
 *
 *  Created on: July 2024
 *      Author: JakNiem
 */

#include "SphericalSampling.h"

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


SphericalSampling::SphericalSampling() {}

void SphericalSampling ::readXML(XMLfileUnits& xmlconfig) {

    xmlconfig.getNodeValue("nShells", _nShells);
    xmlconfig.getNodeValue("start", _startSampling);
    xmlconfig.getNodeValue("writefrequency", _writeFrequency);
    xmlconfig.getNodeValue("stop", _stopSampling);

    Log::global_log->info() << getPluginName() << "  Start:WriteFreq:Stop: " << _startSampling << " : " << _writeFrequency << " : " << _stopSampling << std::endl;
    Log::global_log->info() << getPluginName() << "  nShells: " << _nShells << std::endl;
    Log::global_log->info() << getPluginName() << "  all components treated as one" << std::endl;
}


void SphericalSampling::init(ParticleContainer* /* particleContainer */, DomainDecompBase* domainDecomp, Domain* domain) {

    // Entry per bin; all components sampled as one
    _lenVector = _nShells + 1;
    resizeVectors();
    resetAccumVectors();

    // calculate geometric properties of box and shells:
    _globalBoxLength[0] = domain->getGlobalLength(0);
    _globalBoxLength[1] = domain->getGlobalLength(1);
    _globalBoxLength[2] = domain->getGlobalLength(2);
    _distMax = *(std::min_element(_globalBoxLength.begin(), _globalBoxLength.end())) * 0.5;
    std::transform(_globalBoxLength.begin(), _globalBoxLength.end(), _globalCenter.begin(), [](auto& c){return c*.5;});

    _shellMappingDistMax = _distMax*_distMax; // for position-to-shell mapping by squared distance from center

    /* idea: shells are distributed not with constant thickness, but rather with the lower bound of shell i is proportional to sqrt(i). 
     * This makes shell volume less dependent on shell position (tiny innermost shells have bad statistics) 
     */ 
    //this implementation is problably bad. perhaps better using stl-algorithm? (#todo):
    _shell_lowerBound[0] = 0;
    const double nShellsInv = 1./_nShells;
    for(int i = 1; i < _lenVector; i++){
        _shell_lowerBound[i] = std::sqrt( static_cast<double>(i)*nShellsInv*_shellMappingDistMax );
        _shell_centralRadius[i-1] = (_shell_lowerBound[i-1] + _shell_lowerBound[i])/2.;
        _shell_volume[i-1] = 4./3. * M_PI * ( std::pow(_shell_lowerBound[i],3) - std::pow(_shell_lowerBound[i-1],3) );
    }
    _shell_centralRadius[_lenVector-1] = (_shell_lowerBound[_lenVector-1] + (_distMax/_nShells)); // not very precise
    _shell_volume[_lenVector-1] = _globalBoxLength[0]*_globalBoxLength[1]*_globalBoxLength[2] - 4./3.*M_PI*std::pow(_shell_lowerBound[_lenVector-1],3);
    // \end bad implementation
}


void SphericalSampling::afterForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) {

    // Sampling starts after _startSampling and is conducted up to _stopSampling
    if ((simstep <= _startSampling) or (simstep > _stopSampling)) {
        return;
    }

    // // Do not write or sample data directly after (re)start in the first step
    // if ( simstep == _simulation.getNumInitTimesteps() ) {
    //     return;
    // }

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
    CommVar<std::vector<double>> ekin_step;       // Including drift energy
    CommVar<std::vector<double>> virN_step;
    CommVar<std::vector<double>> virT_step;
    CommVar<std::vector<double>> velocityN_step;
    std::array<CommVar<std::vector<double>>, 3> ekinVect_step;
    std::array<CommVar<std::vector<double>>, 3> virialVect_step;

    numMolecules_step.local.resize(_lenVector); // shouldn't there be a cleaner way to do this?
    mass_step.local.resize(_lenVector);
    ekin_step.local.resize(_lenVector);
    virN_step.local.resize(_lenVector);
    virT_step.local.resize(_lenVector);
    velocityN_step.local.resize(_lenVector);

    numMolecules_step.global.resize(_lenVector);
    mass_step.global.resize(_lenVector);
    ekin_step.global.resize(_lenVector);
    virN_step.global.resize(_lenVector);
    virT_step.global.resize(_lenVector);
    velocityN_step.global.resize(_lenVector);

    for (unsigned short d = 0; d < 3; d++) {
        ekinVect_step[d].local.resize(_lenVector);
        virialVect_step[d].local.resize(_lenVector);
        ekinVect_step[d].global.resize(_lenVector);
        virialVect_step[d].global.resize(_lenVector);
    }

    std::fill(numMolecules_step.local.begin(), numMolecules_step.local.end(), 0ul);
    std::fill(mass_step.local.begin(),         mass_step.local.end(), 0.);
    std::fill(ekin_step.local.begin(),         ekin_step.local.end(), 0.);
    std::fill(virN_step.local.begin(),         virN_step.local.end(), 0.);
    std::fill(virT_step.local.begin(),         virT_step.local.end(), 0.);
    std::fill(velocityN_step.local.begin(),    velocityN_step.local.end(), 0.);

    std::fill(numMolecules_step.global.begin(), numMolecules_step.global.end(), 0ul);
    std::fill(mass_step.global.begin(),   mass_step.global.end(), 0.);
    std::fill(ekin_step.global.begin(),   ekin_step.global.end(), 0.);
    std::fill(virN_step.global.begin(),   virN_step.global.end(), 0.);
    std::fill(virT_step.global.begin(),   virT_step.global.end(), 0.);
    std::fill(velocityN_step.global.begin(),   velocityN_step.global.end(), 0.);

    for (unsigned short d = 0; d < 3; d++) {
        std::fill(ekinVect_step[d].local.begin(),      ekinVect_step[d].local.end(), 0.);
        std::fill(virialVect_step[d].local.begin(),   virialVect_step[d].local.end(), 0.);

        std::fill(ekinVect_step[d].global.begin(),      ekinVect_step[d].global.end(), 0.);
        std::fill(virialVect_step[d].global.begin(),   virialVect_step[d].global.end(), 0.);
    }

    for (auto pit = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); pit.isValid(); ++pit) {
        const double distCenter_x = pit->r(0) - (_globalCenter[0]);
        const double distCenter_y = pit->r(1) - (_globalCenter[1]);
        const double distCenter_z = pit->r(2) - (_globalCenter[2]);
        const double distCenterSquared = std::pow(distCenter_x,2) + std::pow(distCenter_y,2) + std::pow(distCenter_z,2);
        const double distCenter = std::sqrt(distCenterSquared);
        // // Do not consider particles outside of most outer radius
        // if (distCenterSquared >= _distMax*_distMax) { continue; }
        const double distMaxSquared = _distMax * _distMax;


        const unsigned int index = std::min(_nShells, static_cast<unsigned int>((distCenterSquared/_shellMappingDistMax)*_nShells));  // Index of bin of radius

        numMolecules_step.local[index]++;

        const double u_x = pit->v(0);
        const double u_y = pit->v(1);
        const double u_z = pit->v(2);
        const double mass = pit->mass();
        const double vi_x = pit->Vi(0);
        const double vi_y = pit->Vi(1);
        const double vi_z = pit->Vi(2);
        const double vi_n = pit->ViN();
        const double vi_t = pit->ViT();

        mass_step.local[index] += mass;
        ekin_step.local[index] += pit->U_kin();
        
        velocityN_step.local[index] = (distCenter_x*u_x + distCenter_y*u_y + distCenter_z*u_z)/distCenter;  // Radial 

        virialVect_step[0].local[index] += vi_x;
        virialVect_step[1].local[index] += vi_y;
        virialVect_step[2].local[index] += vi_z;
        virN_step.local[index] += vi_n;
        virT_step.local[index] += vi_t;
        
        ekinVect_step[0].local[index] += 0.5*mass*u_x*u_x;
        ekinVect_step[1].local[index] += 0.5*mass*u_y*u_y;
        ekinVect_step[2].local[index] += 0.5*mass*u_z*u_z;
    }

// Gather quantities. Note: MPI_Reduce instead of MPI_Allreduce! Therefore, only root has correct values
#ifdef ENABLE_MPI
    MPI_Reduce(numMolecules_step.local.data(), numMolecules_step.global.data(), _lenVector, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(mass_step.local.data(), mass_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekin_step.local.data(), ekin_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(velocityN_step.local.data(), velocityN_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekinVect_step[0].local.data(), ekinVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekinVect_step[1].local.data(), ekinVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(ekinVect_step[2].local.data(), ekinVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[0].local.data(), virialVect_step[0].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[1].local.data(), virialVect_step[1].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virialVect_step[2].local.data(), virialVect_step[2].global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virN_step.local.data(), virN_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(virT_step.local.data(), virT_step.global.data(), _lenVector, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    for (unsigned long i = 0; i < _lenVector; i++) {
        numMolecules_step.global[i] = numMolecules_step.local[i];
        mass_step.global[i] = mass_step.local[i];
        ekin_step.global[i] = ekin_step.local[i];
        velocityN_step.global[i] = velocityN_step.local[i];
        virialVect_step[0].global[i] = virialVect_step[0].local[i];
        virialVect_step[1].global[i] = virialVect_step[1].local[i];
        virialVect_step[2].global[i] = virialVect_step[2].local[i];
        virN_step.global[i] = virN_step.local[i];
        virT_step.global[i] = virT_step.local[i];
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
            _virN_accum[i]                    += .5 * virN_step.global[i];
            _virT_accum[i]                    += .5 * virT_step.global[i];

            _virialVect_accum[0][i]           += vi_x;
            _virialVect_accum[1][i]           += vi_y;
            _virialVect_accum[2][i]           += vi_z;
            
            _ekinVect_accum[0][i]             += ekinVect_step[0].global[i];
            _ekinVect_accum[1][i]             += ekinVect_step[1].global[i];
            _ekinVect_accum[2][i]             += ekinVect_step[2].global[i];

            if (numMols > 0ul) {
                _velocityN_accum[i]           += velocityN_step.global[i] / numMols;
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
            const std::string fname = "SphericalSampling_TS"+ss.str()+".dat";
            std::ofstream ofs;
            ofs.open(fname, std::ios::out);
            ofs << std::setw(24) << "radius_(mid)"     // Bin position (radius)
                << std::setw(24) << "rho"        // Density
                << std::setw(24) << "p_(xyzVir)"   // Pressure (using xyz vir)
                << std::setw(24) << "p_(sphVir)"  // Pressure (using spherical vir)
                << std::setw(24) << "p_n"        // Pressure in normal direction; 
                << std::setw(24) << "p_t"        // Pressure in tangantial direction; 
                << std::setw(24) << "VirN"          // virial in normal direction
                << std::setw(24) << "VirT"          // virial in tangential direction
                << std::setw(24) << "T"             // Temperature, assuming that there is no drift
                << std::setw(24) << "T_driftcorr"   // Temperature without radial drift (i.e. "real" temperature; given that there is no non-radial drift)
                << std::setw(24) << "shellVolume"   // Volume of Shell
                << std::setw(24) << "numParts"   // Average number of molecules in bin per step
                << std::setw(24) << "ekin"       // Kinetic energy including drift
                << std::setw(24) << "v_r"        // drift velocity in radial direction
                << std::setw(24) << "numSamples";    // Number of samples (<= _writeFrequency)
            ofs << std::endl;


            for (unsigned long i = 0; i < _lenVector; i++) {
                unsigned long numSamples = _countSamples[i];
                double numMolsPerStep {std::nan("0")}; // Not an int as particles change bin during simulation
                double rho {};
                double T {std::nan("0")};            //initialized to NaN, because should be NaN in case no particles are in the shell 
                double T_driftcorr {std::nan("0")};
                double ekin {std::nan("0")};
                double p_xyz {std::nan("0")};
                double p_sph {std::nan("0")};
                double v_r {std::nan("0")};
                double p_n {std::nan("0")};
                double p_t {std::nan("0")};
                double vir_n {std::nan("0")};
                double vir_t {std::nan("0")};
                 
                if ((numSamples > 0ul) and (_numMolecules_accum[i] > 0ul) and (_doftotal_accum[i]) > 0ul) {
                    const double numMols_accum = static_cast<double>(_numMolecules_accum[i]);

                    numMolsPerStep = numMols_accum    /numSamples;
                    rho         = numMolsPerStep      / _shell_volume[i];
                    v_r         = _velocityN_accum[i] / numSamples;

                    double v_drift_squared = v_r*v_r;  // considering only radial drift

                    vir_n = _virN_accum[i]/numMols_accum;
                    vir_t = _virT_accum[i]/numMols_accum;

                    T           = (2*_ekin_accum[i]) / _doftotal_accum[i];
                    T_driftcorr = (2*_ekin_accum[i] - v_drift_squared*_mass_accum[i]) / _doftotal_accum[i];
                    ekin        = _ekin_accum[i] / numMols_accum;
                    p_xyz       = rho * ( (_virialVect_accum[0][i]+_virialVect_accum[1][i]+_virialVect_accum[2][i])/(3.0*numMols_accum) + T);
                    p_sph       = rho * ( T + (vir_n + 2.*vir_t)/(3.) );

                    p_n         = rho * (T + vir_n);
                    p_t         = rho * (T + vir_t);
                }
                         
                ofs << FORMAT_SCI_MAX_DIGITS << _shell_centralRadius[i]  // Radius bin
                    << FORMAT_SCI_MAX_DIGITS << rho
                    << FORMAT_SCI_MAX_DIGITS << p_xyz
                    << FORMAT_SCI_MAX_DIGITS << p_sph
                    << FORMAT_SCI_MAX_DIGITS << p_n
                    << FORMAT_SCI_MAX_DIGITS << p_t
                    << FORMAT_SCI_MAX_DIGITS << vir_n
                    << FORMAT_SCI_MAX_DIGITS << vir_t
                    << FORMAT_SCI_MAX_DIGITS << T
                    << FORMAT_SCI_MAX_DIGITS << T_driftcorr
                    << FORMAT_SCI_MAX_DIGITS << _shell_volume[i]
                    << FORMAT_SCI_MAX_DIGITS << numMolsPerStep
                    << FORMAT_SCI_MAX_DIGITS << ekin
                    << FORMAT_SCI_MAX_DIGITS << v_r
                    << FORMAT_SCI_MAX_DIGITS << numSamples
                    << std::endl;
            }
            // ofs.close(); #automatically handled by std::ofstream's destructor.
        }
        // Reset accumulation-vectors to zero
        resetAccumVectors();
    }
}

// Resize vectors
void SphericalSampling ::resizeVectors() {

    _numMolecules_accum.resize(_lenVector);
    _doftotal_accum.resize(_lenVector);
    _mass_accum.resize(_lenVector);
    _ekin_accum.resize(_lenVector);
    _virN_accum.resize(_lenVector);              
    _virT_accum.resize(_lenVector);              
    _velocityN_accum.resize(_lenVector);

    for (unsigned short d = 0; d < 3; d++) {
        _ekinVect_accum[d].resize(_lenVector);
        _virialVect_accum[d].resize(_lenVector);
    }

    _shell_lowerBound.resize(_lenVector);
    _shell_centralRadius.resize(_lenVector);
    _shell_volume.resize(_lenVector);

    _countSamples.resize(_lenVector);
}

// Fill vectors with zeros. take care not to list vectors here, that are not supposed to be reset.
void SphericalSampling ::resetAccumVectors() {
    std::fill(_numMolecules_accum.begin(), _numMolecules_accum.end(), 0ul);
    std::fill(_doftotal_accum.begin(), _doftotal_accum.end(), 0ul);
    std::fill(_mass_accum.begin(), _mass_accum.end(), 0.0f);
    std::fill(_ekin_accum.begin(), _ekin_accum.end(), 0.0f);
    std::fill(_virN_accum.begin(), _virN_accum.end(), 0.0f);
    std::fill(_virT_accum.begin(), _virT_accum.end(), 0.0f);
    std::fill(_velocityN_accum.begin(), _velocityN_accum.end(), 0.0f);
    // std::fill(_virial_accum.begin(), _virial_accum.end(), 0.0f);

    for (unsigned short d = 0; d < 3; d++) {
        std::fill(_ekinVect_accum[d].begin(), _ekinVect_accum[d].end(), 0.0f);
        std::fill(_virialVect_accum[d].begin(), _virialVect_accum[d].end(), 0.0f);
    }

    std::fill(_countSamples.begin(), _countSamples.end(), 0ul);
}