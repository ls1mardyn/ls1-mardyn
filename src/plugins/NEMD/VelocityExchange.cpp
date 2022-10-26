
#include "plugins/NEMD/VelocityExchange.h"

#include <vector>
#include <cmath>
#include <algorithm>

#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"
#include "utils/CommVar.h"
#include "utils/xmlfileUnits.h"
#include "utils/FileUtils.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif


VelocityExchange::VelocityExchange() {}

void VelocityExchange::init(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* domain) {

    for (unsigned short d = 0; d < 3; d++) {
        _boxLength[d] = domain->getGlobalLength(d);
    }

    _numComp = global_simulation->getEnsemble()->getComponents()->size();
}

void VelocityExchange::readXML(XMLfileUnits& xmlconfig) {
    // Timestep control
    xmlconfig.getNodeValue("control/start", _control.start);
    xmlconfig.getNodeValue("control/frequency", _control.freq);
    xmlconfig.getNodeValue("control/stop", _control.stop);
    global_log->info() << "[VelocityExchange] takes place start:freq:stop = " << _control.start << ":" << _control.freq << ":" << _control.stop << std::endl;

    // regions
    _cold_region.xmax = _boxLength[0];
    _cold_region.ymax = _boxLength[1];
    _cold_region.zmax = _boxLength[2];
    _warm_region.xmax = _boxLength[0];
    _warm_region.ymax = _boxLength[1];
    _warm_region.zmax = _boxLength[2];
    _symmetry = false;

    xmlconfig.getNodeValue("coldregion/xmin", _cold_region.xmin);
    xmlconfig.getNodeValue("coldregion/ymin", _cold_region.ymin);
    xmlconfig.getNodeValue("coldregion/zmin", _cold_region.zmin);
    xmlconfig.getNodeValue("coldregion/xmax", _cold_region.xmax);
    xmlconfig.getNodeValue("coldregion/ymax", _cold_region.ymax);
    xmlconfig.getNodeValue("coldregion/zmax", _cold_region.zmax);

    xmlconfig.getNodeValue("warmregion/symmetric", _symmetry);
    xmlconfig.getNodeValue("warmregion/xmin", _warm_region.xmin);
    xmlconfig.getNodeValue("warmregion/zmin", _warm_region.zmin);
    xmlconfig.getNodeValue("warmregion/ymin", _warm_region.ymin);
    xmlconfig.getNodeValue("warmregion/xmax", _warm_region.xmax);
    xmlconfig.getNodeValue("warmregion/ymax", _warm_region.ymax);
    xmlconfig.getNodeValue("warmregion/zmax", _warm_region.zmax);

    global_log->info() << "[VelocityExchange] Cold region:"
                        << " x = " << _cold_region.xmin << " - " << _cold_region.xmax << " ;"
                        << " y = " << _cold_region.ymin << " - " << _cold_region.ymax << " ;"
                        << " z = " << _cold_region.zmin << " - " << _cold_region.zmax << endl;

    global_log->info() << "[VelocityExchange] Warm region" << ((_symmetry) ? " (left)" : "") << ":"
                        << " x = " << _cold_region.xmin << " - " << _cold_region.xmax << " ;"
                        << " y = " << _cold_region.ymin << " - " << _cold_region.ymax << " ;"
                        << " z = " << _cold_region.zmin << " - " << _cold_region.zmax << endl;

    if (_symmetry) {
        global_log->info() << "[VelocityExchange] Warm region (right):"
                        << " x = " << _cold_region.xmin << " - " << _cold_region.xmax << " ;"
                        << " y = " << _boxLength[1]-_cold_region.ymax << " - " << _boxLength[1]-_cold_region.ymin << " ;"
                        << " z = " << _cold_region.zmin << " - " << _cold_region.zmax << endl;
    }
}

void VelocityExchange::beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) {
    if (simstep < _control.start || simstep > _control.stop || simstep % _control.freq != 0)
        return;
    this->exchangeVelocities(particleContainer, domainDecomp);
}

void VelocityExchange::exchangeVelocities(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp) {

    // Absolute velocities of coldest (c) and warmest (w) particle per component
    CommVar<std::vector<double>> velocity_abs_coldP;
    CommVar<std::vector<double>> velocity_abs_warmP;

    // Molecule ID of coldest (c) and warmest (w) particle per component
    CommVar<std::vector<unsigned long>> molID_coldP;
    CommVar<std::vector<unsigned long>> molID_warmP;

    // Direction-wise velocities of coldest (c) and warmest (w) particle per component; vector contains [c1_x,c1_y,c1_z,c2_y,c2_z,c2_x,etc.]
    CommVar<std::array<std::vector<double>, 3>> velocity_coldP;
    CommVar<std::array<std::vector<double>, 3>> velocity_warmP;
    CommVar<std::array<std::vector<double>, 3>> rotVelo_coldP;
    CommVar<std::array<std::vector<double>, 3>> rotVelo_warmP;

    velocity_abs_coldP.local.resize(_numComp);
    velocity_abs_coldP.global.resize(_numComp);
    velocity_abs_warmP.local.resize(_numComp);
    velocity_abs_warmP.global.resize(_numComp);

    molID_coldP.local.resize(_numComp);
    molID_coldP.global.resize(_numComp);
    molID_warmP.local.resize(_numComp);
    molID_warmP.global.resize(_numComp);

    std::fill(velocity_abs_coldP.local.begin(),  velocity_abs_coldP.local.end(), 1000.0f);  // Set to high number for calculation reasons
    std::fill(velocity_abs_coldP.global.begin(), velocity_abs_coldP.global.end(),   0.0f);
    std::fill(velocity_abs_warmP.local.begin(),  velocity_abs_warmP.local.end(),    0.0f);
    std::fill(velocity_abs_warmP.global.begin(), velocity_abs_warmP.global.end(),   0.0f);

    std::fill(molID_coldP.local.begin(),  molID_coldP.local.end(),  0ul);
    std::fill(molID_coldP.global.begin(), molID_coldP.global.end(), 0ul);
    std::fill(molID_warmP.local.begin(),  molID_warmP.local.end(),  0ul);
    std::fill(molID_warmP.global.begin(), molID_warmP.global.end(), 0ul);

    for (unsigned short d = 0; d < 3; d++) {
        velocity_coldP.local[d].resize(_numComp);
        velocity_coldP.global[d].resize(_numComp);
        velocity_warmP.local[d].resize(_numComp);
        velocity_warmP.global[d].resize(_numComp);
        rotVelo_coldP.local[d].resize(_numComp);
        rotVelo_coldP.global[d].resize(_numComp);
        rotVelo_warmP.local[d].resize(_numComp);
        rotVelo_warmP.global[d].resize(_numComp);

        std::fill(velocity_coldP.local[d].begin(),  velocity_coldP.local[d].end(),  -100000.0f);
        std::fill(velocity_coldP.global[d].begin(), velocity_coldP.global[d].end(), -100000.0f);
        std::fill(velocity_warmP.local[d].begin(),  velocity_warmP.local[d].end(),  -100000.0f);
        std::fill(velocity_warmP.global[d].begin(), velocity_warmP.global[d].end(), -100000.0f);
        std::fill(rotVelo_coldP.local[d].begin(),  rotVelo_coldP.local[d].end(),  -100000.0f);
        std::fill(rotVelo_coldP.global[d].begin(), rotVelo_coldP.global[d].end(), -100000.0f);
        std::fill(rotVelo_warmP.local[d].begin(),  rotVelo_warmP.local[d].end(),  -100000.0f);
        std::fill(rotVelo_warmP.global[d].begin(), rotVelo_warmP.global[d].end(), -100000.0f);
    }


    // Corners of cold region
    double coldRegion_lowCorner[3], coldRegion_highCorner[3];
    // if linked cell in the region of interest
    for (unsigned short d = 0; d < 3; d++) {
        coldRegion_lowCorner[d] = particleContainer->getBoundingBoxMin(d);
        coldRegion_highCorner[d] = particleContainer->getBoundingBoxMax(d);
    }
    // ensure that we do not iterate over things outside of the container.
    coldRegion_lowCorner[0] = std::max(_cold_region.xmin, coldRegion_lowCorner[0]);
    coldRegion_lowCorner[1] = std::max(_cold_region.ymin, coldRegion_lowCorner[1]);
    coldRegion_lowCorner[2] = std::max(_cold_region.zmin, coldRegion_lowCorner[2]);
    coldRegion_highCorner[0] = std::min(_cold_region.xmax, coldRegion_highCorner[0]);
    coldRegion_highCorner[1] = std::min(_cold_region.ymax, coldRegion_highCorner[1]);
    coldRegion_highCorner[2] = std::min(_cold_region.zmax, coldRegion_highCorner[2]);
    auto begin_c = particleContainer->regionIterator(coldRegion_lowCorner, coldRegion_highCorner, ParticleIterator::ONLY_INNER_AND_BOUNDARY);


    // Corners of warm region
    double warmRegion_lowCorner[3], warmRegion_highCorner[3];
    // if linked cell in the region of interest
    for (unsigned short d = 0; d < 3; d++) {
        warmRegion_lowCorner[d] = particleContainer->getBoundingBoxMin(d);
        warmRegion_highCorner[d] = particleContainer->getBoundingBoxMax(d);
    }
    // ensure that we do not iterate over things outside of the container.
    warmRegion_lowCorner[0] = std::max(_warm_region.xmin, warmRegion_lowCorner[0]);
    warmRegion_lowCorner[1] = std::max(_warm_region.ymin, warmRegion_lowCorner[1]);
    warmRegion_lowCorner[2] = std::max(_warm_region.zmin, warmRegion_lowCorner[2]);
    warmRegion_highCorner[0] = std::min(_warm_region.xmax, warmRegion_highCorner[0]);
    warmRegion_highCorner[1] = std::min(_warm_region.ymax, warmRegion_highCorner[1]);
    warmRegion_highCorner[2] = std::min(_warm_region.zmax, warmRegion_highCorner[2]);
    auto begin_w = particleContainer->regionIterator(warmRegion_lowCorner, warmRegion_highCorner, ParticleIterator::ONLY_INNER_AND_BOUNDARY);


    // Corners of symmetric warm region
    double warmRegion_sym_lowCorner[3], warmRegion_sym_highCorner[3];
    double ymin_symm = _boxLength[1] - _warm_region.ymax;
    double ymax_symm = _boxLength[1] - _warm_region.ymin;
    // if linked cell in the region of interest
    for (unsigned short d = 0; d < 3; d++) {
        warmRegion_sym_lowCorner[d] = particleContainer->getBoundingBoxMin(d);
        warmRegion_sym_highCorner[d] = particleContainer->getBoundingBoxMax(d);
    }
    // ensure that we do not iterate over things outside of the container.
    warmRegion_sym_lowCorner[0] = std::max(_cold_region.xmin, warmRegion_sym_lowCorner[0]);
    warmRegion_sym_lowCorner[1] = std::max(ymin_symm, warmRegion_sym_lowCorner[1]);
    warmRegion_sym_lowCorner[2] = std::max(_cold_region.zmin, warmRegion_sym_lowCorner[2]);
    warmRegion_sym_highCorner[0] = std::min(_cold_region.xmax, warmRegion_sym_highCorner[0]);
    warmRegion_sym_highCorner[1] = std::min(ymax_symm, warmRegion_sym_highCorner[1]);
    warmRegion_sym_highCorner[2] = std::min(_cold_region.zmax, warmRegion_sym_highCorner[2]);
    auto begin_w_s = particleContainer->regionIterator(warmRegion_sym_lowCorner, warmRegion_sym_highCorner, ParticleIterator::ONLY_INNER_AND_BOUNDARY);


    double v_dummy[3] = {0.0f, 0.0f, 0.0f};

    std::vector<Molecule*> warm_mol_ptr(_numComp);
    std::vector<Molecule*> cold_mol_ptr(_numComp);

    // find warmest mol in cold region
    for (auto it = begin_c; it.isValid(); ++it) {
        uint32_t cid = it->componentid();
        v_dummy[0] = it->v(0);
        v_dummy[1] = it->v(1);
        v_dummy[2] = it->v(2);
        double v_dummy_abs = v_dummy[0]*v_dummy[0]+v_dummy[1]*v_dummy[1]+v_dummy[2]*v_dummy[2];
        if ( v_dummy_abs > velocity_abs_warmP.local[cid] ) {
            velocity_abs_warmP.local[cid] = v_dummy_abs;
            warm_mol_ptr[cid] = &(*it);
        }
    }

    // find coldest mol in warm region
    for (auto it = begin_w; it.isValid(); ++it) {
        uint32_t cid = it->componentid();
        v_dummy[0] = it->v(0);
        v_dummy[1] = it->v(1);
        v_dummy[2] = it->v(2);
        double v_dummy_abs = v_dummy[0]*v_dummy[0]+v_dummy[1]*v_dummy[1]+v_dummy[2]*v_dummy[2];
        if ( v_dummy_abs < velocity_abs_coldP.local[cid] ) {
            velocity_abs_coldP.local[cid] = v_dummy_abs;
            cold_mol_ptr[cid] = &(*it);
        }
    }

    if ( _symmetry ) {  // warm region is symmetric
        for (auto it = begin_w; it.isValid(); ++it) {
            uint32_t cid = it->componentid();
            v_dummy[0] = it->v(0);
            v_dummy[1] = it->v(1);
            v_dummy[2] = it->v(2);
            double v_dummy_abs = v_dummy[0]*v_dummy[0]+v_dummy[1]*v_dummy[1]+v_dummy[2]*v_dummy[2];
            if ( v_dummy_abs < velocity_abs_coldP.local[cid] ) {
                velocity_abs_coldP.local[cid] = v_dummy_abs;
                cold_mol_ptr[cid] = &(*it);
            }
        }
    }

#ifdef ENABLE_MPI

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allreduce(velocity_abs_coldP.local.data(), velocity_abs_coldP.global.data(), _numComp, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(velocity_abs_warmP.local.data(), velocity_abs_warmP.global.data(), _numComp, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

#else

    velocity_abs_coldP.global = velocity_abs_coldP.local;
    velocity_abs_warmP.global = velocity_abs_warmP.local;

#endif

#ifndef NDEBUG
    // find rank with coldest/warmest molecule and read velocities
    for (uint32_t cid = 0; cid < _numComp; cid++) {
        std::cout << "[VelocityExchange] component " << cid << " ; rank " << domainDecomp->getRank() << " ; velocity_abs_coldP.local = " << velocity_abs_coldP.local[cid] << std::endl;
        std::cout << "[VelocityExchange] component " << cid << " ; rank " << domainDecomp->getRank() << " ; velocity_abs_warmP.local = " << velocity_abs_warmP.local[cid] << std::endl;
        std::cout << "[VelocityExchange] component " << cid << " ; rank " << domainDecomp->getRank() << " ; velocity_abs_coldP.global = " << velocity_abs_coldP.global[cid] << std::endl;
        std::cout << "[VelocityExchange] component " << cid << " ; rank " << domainDecomp->getRank() << " ; velocity_abs_warmP.global = " << velocity_abs_warmP.global[cid] << std::endl;
    }
#endif

    // Seek the molecule with the highest abs. velocity and store its ID and direction-wise velos
    for (uint32_t cid = 0; cid < _numComp; cid++) {
        // Warm region
        if ( (abs(velocity_abs_coldP.local[cid]-velocity_abs_coldP.global[cid]) < 1e-9) && (velocity_abs_coldP.local[cid] < 1000.0) ) {
            molID_coldP.local[cid] = cold_mol_ptr[cid]->getID();
            for (unsigned short d = 0; d < 3; d++) {
                velocity_coldP.local[d].at(cid) = cold_mol_ptr[cid]->v(d);
                rotVelo_coldP.local[d].at(cid) = cold_mol_ptr[cid]->D(d);
            }
        }

        // Cold region
        if ( (abs(velocity_abs_warmP.local[cid]-velocity_abs_warmP.global[cid]) < 1e-9) && (velocity_abs_warmP.local[cid] > 0.0) ) {
            molID_warmP.local[cid] = warm_mol_ptr[cid]->getID();
            for (unsigned short d = 0; d < 3; d++) {
                velocity_warmP.local[d].at(cid) = warm_mol_ptr[cid]->v(d);
                rotVelo_warmP.local[d].at(cid) = warm_mol_ptr[cid]->D(d);
            }
        }
    }


// distribute IDs and velocities to all ranks using MPI_MAX as only one rank contains real information
#ifdef ENABLE_MPI
    MPI_Allreduce(molID_warmP.local.data(), molID_warmP.global.data(), _numComp, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(molID_coldP.local.data(), molID_coldP.global.data(), _numComp, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    for (unsigned short d = 0; d < 3; d++) {
        MPI_Allreduce(velocity_coldP.local[d].data(), velocity_coldP.global[d].data(), _numComp, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(velocity_warmP.local[d].data(), velocity_warmP.global[d].data(), _numComp, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(rotVelo_coldP.local[d].data(), rotVelo_coldP.global[d].data(), _numComp, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(rotVelo_warmP.local[d].data(), rotVelo_warmP.global[d].data(), _numComp, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
#else
    molID_warmP.global = molID_warmP.local;
    molID_coldP.global = molID_coldP.local;

    for (unsigned short d = 0; d < 3; d++) {
        velocity_coldP.local[d] = velocity_coldP.global[d];
        velocity_warmP.local[d] = velocity_warmP.global[d];
        rotVelo_coldP.local[d] = rotVelo_coldP.global[d];
        rotVelo_warmP.local[d] = rotVelo_warmP.global[d];
    }
#endif

    // Assign velocities of coldest particle in warm region to warmest particle in cold region
    for (auto it = begin_c; it.isValid(); ++it) {
        uint32_t cid = it->componentid();
        if ( it->getID() == molID_warmP.global[cid] ) {
            for (unsigned short d = 0; d < 3; d++) {
                it->setv(d, velocity_coldP.global[d].at(cid));
                it->setD(d, rotVelo_coldP.global[d].at(cid));
            }
        }
    }

    // Assign velocities of warmest particle in cold region to coldest particle in warm region
    for (auto it = begin_w; it.isValid(); ++it) {
        uint32_t cid = it->componentid();
        if ( it->getID() == molID_coldP.global[cid] ) {
            for (unsigned short d = 0; d < 3; d++) {
                it->setv(d, velocity_warmP.global[d].at(cid));
                it->setD(d, rotVelo_warmP.global[d].at(cid));
            }
        }
    }

    if ( _symmetry ) {
        for (auto it = begin_w_s; it.isValid(); ++it) {
            uint32_t cid = it->componentid();
            if ( it->getID() == molID_coldP.global[cid] ) {
                for (unsigned short d = 0; d < 3; d++) {
                    it->setv(d, velocity_warmP.global[d].at(cid));
                    it->setD(d, rotVelo_warmP.global[d].at(cid));
                }
            }
        }
    }

    for (uint32_t cid = 0; cid < _numComp; cid++) {
        global_log->info() << "[VelocityExchange] flipped velocities of molecules "
                            << molID_coldP.global[cid] << " (in warm region, v_new = " << velocity_abs_warmP.global[cid] << ") and "
                            << molID_warmP.global[cid] << " (in cold region v_new = " << velocity_abs_coldP.global[cid] << ") of component " << cid+1 << std::endl;
    }
}
