
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

    // range
    _cold_range.xmax = _boxLength[0];
    _cold_range.ymax = _boxLength[1];
    _cold_range.zmax = _boxLength[2];
    _warm_range.xmax = _boxLength[0];
    _warm_range.ymax = _boxLength[1];
    _warm_range.zmax = _boxLength[2];
    _symmetry = false;

    xmlconfig.getNodeValue("coldrange/xmin", _cold_range.xmin);
    xmlconfig.getNodeValue("coldrange/ymin", _cold_range.ymin);
    xmlconfig.getNodeValue("coldrange/zmin", _cold_range.zmin);
    xmlconfig.getNodeValue("coldrange/xmax", _cold_range.xmax);
    xmlconfig.getNodeValue("coldrange/ymax", _cold_range.ymax);
    xmlconfig.getNodeValue("coldrange/zmax", _cold_range.zmax);

    xmlconfig.getNodeValue("warmrange/symmetric", _symmetry);
    xmlconfig.getNodeValue("warmrange/xmin", _warm_range.xmin);
    xmlconfig.getNodeValue("warmrange/zmin", _warm_range.zmin);
    xmlconfig.getNodeValue("warmrange/ymin", _warm_range.ymin);
    xmlconfig.getNodeValue("warmrange/xmax", _warm_range.xmax);
    xmlconfig.getNodeValue("warmrange/ymax", _warm_range.ymax);
    xmlconfig.getNodeValue("warmrange/zmax", _warm_range.zmax);

    global_log->info() << "[VelocityExchange] Cold region:"
                        << " x = " << _cold_range.xmin << " - " << _cold_range.xmax << " ;"
                        << " y = " << _cold_range.ymin << " - " << _cold_range.ymax << " ;"
                        << " z = " << _cold_range.zmin << " - " << _cold_range.zmax << endl;

    global_log->info() << "[VelocityExchange] Warm region" << ((_symmetry) ? " (left)" : "") << ":"
                        << " x = " << _cold_range.xmin << " - " << _cold_range.xmax << " ;"
                        << " y = " << _cold_range.ymin << " - " << _cold_range.ymax << " ;"
                        << " z = " << _cold_range.zmin << " - " << _cold_range.zmax << endl;

    if (_symmetry) {
        global_log->info() << "[VelocityExchange] Warm region (right):"
                        << " x = " << _cold_range.xmin << " - " << _cold_range.xmax << " ;"
                        << " y = " << _boxLength[1]-_cold_range.ymax << " - " << _boxLength[1]-_cold_range.ymin << " ;"
                        << " z = " << _cold_range.zmin << " - " << _cold_range.zmax << endl;
    }
}

void VelocityExchange::beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep) {
    if (simstep < _control.start || simstep > _control.stop || simstep % _control.freq != 0)
        return;
    this->exchangeVelocities(particleContainer, domainDecomp);
}

void VelocityExchange::exchangeVelocities(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp) {

    // Absolute velocities of coldest (c) and warmest (w) particle per component
    CommVar<std::vector<double>> v_c_abs;
    CommVar<std::vector<double>> v_w_abs;

    // Molecule ID of coldest (c) and warmest (w) particle per component
    CommVar<std::vector<unsigned long>> mol_c;
    CommVar<std::vector<unsigned long>> mol_w;

    // Direction-wise velocities of coldest (c) and warmest (w) particle per component; vector contains [c1_x,c1_y,c1_z,c2_y,c2_z,c2_x,etc.]
    CommVar<std::array<std::vector<double>, 3>> v_c;
    CommVar<std::array<std::vector<double>, 3>> v_w;
    CommVar<std::array<std::vector<double>, 3>> D_c;
    CommVar<std::array<std::vector<double>, 3>> D_w;

    v_c_abs.local.resize(_numComp);
    v_c_abs.global.resize(_numComp);
    v_w_abs.local.resize(_numComp);
    v_w_abs.global.resize(_numComp);

    mol_c.local.resize(_numComp);
    mol_c.global.resize(_numComp);
    mol_w.local.resize(_numComp);
    mol_w.global.resize(_numComp);

    std::fill(v_c_abs.local.begin(),  v_c_abs.local.end(), 1000.0f);  // Set to high number for calculation reasons
    std::fill(v_c_abs.global.begin(), v_c_abs.global.end(),   0.0f);
    std::fill(v_w_abs.local.begin(),  v_w_abs.local.end(),    0.0f);
    std::fill(v_w_abs.global.begin(), v_w_abs.global.end(),   0.0f);

    std::fill(mol_c.local.begin(),  mol_c.local.end(),  0ul);
    std::fill(mol_c.global.begin(), mol_c.global.end(), 0ul);
    std::fill(mol_w.local.begin(),  mol_w.local.end(),  0ul);
    std::fill(mol_w.global.begin(), mol_w.global.end(), 0ul);

    for (unsigned short d = 0; d < 3; d++) {
        v_c.local[d].resize(_numComp);
        v_c.global[d].resize(_numComp);
        v_w.local[d].resize(_numComp);
        v_w.global[d].resize(_numComp);
        D_c.local[d].resize(_numComp);
        D_c.global[d].resize(_numComp);
        D_w.local[d].resize(_numComp);
        D_w.global[d].resize(_numComp);

        std::fill(v_c.local[d].begin(),  v_c.local[d].end(),  -100000.0f);
        std::fill(v_c.global[d].begin(), v_c.global[d].end(), -100000.0f);
        std::fill(v_w.local[d].begin(),  v_w.local[d].end(),  -100000.0f);
        std::fill(v_w.global[d].begin(), v_w.global[d].end(), -100000.0f);
        std::fill(D_c.local[d].begin(),  D_c.local[d].end(),  -100000.0f);
        std::fill(D_c.global[d].begin(), D_c.global[d].end(), -100000.0f);
        std::fill(D_w.local[d].begin(),  D_w.local[d].end(),  -100000.0f);
        std::fill(D_w.global[d].begin(), D_w.global[d].end(), -100000.0f);
    }


    // Corners of cold region
    double c_regionLowCorner[3], c_regionHighCorner[3];
    // if linked cell in the region of interest
    for (unsigned short d = 0; d < 3; d++) {
        c_regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
        c_regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
    }
    // ensure that we do not iterate over things outside of the container.
    c_regionLowCorner[0] = std::max(_cold_range.xmin, c_regionLowCorner[0]);
    c_regionLowCorner[1] = std::max(_cold_range.ymin, c_regionLowCorner[1]);
    c_regionLowCorner[2] = std::max(_cold_range.zmin, c_regionLowCorner[2]);
    c_regionHighCorner[0] = std::min(_cold_range.xmax, c_regionHighCorner[0]);
    c_regionHighCorner[1] = std::min(_cold_range.ymax, c_regionHighCorner[1]);
    c_regionHighCorner[2] = std::min(_cold_range.zmax, c_regionHighCorner[2]);
    auto begin_c = particleContainer->regionIterator(c_regionLowCorner, c_regionHighCorner, ParticleIterator::ONLY_INNER_AND_BOUNDARY);


    // Corners of warm region
    double w_regionLowCorner[3], w_regionHighCorner[3];
    // if linked cell in the region of interest
    for (unsigned short d = 0; d < 3; d++) {
        w_regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
        w_regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
    }
    // ensure that we do not iterate over things outside of the container.
    w_regionLowCorner[0] = std::max(_warm_range.xmin, w_regionLowCorner[0]);
    w_regionLowCorner[1] = std::max(_warm_range.ymin, w_regionLowCorner[1]);
    w_regionLowCorner[2] = std::max(_warm_range.zmin, w_regionLowCorner[2]);

    w_regionHighCorner[0] = std::min(_warm_range.xmax, w_regionHighCorner[0]);
    w_regionHighCorner[1] = std::min(_warm_range.ymax, w_regionHighCorner[1]);
    w_regionHighCorner[2] = std::min(_warm_range.zmax, w_regionHighCorner[2]);
    auto begin_w = particleContainer->regionIterator(w_regionLowCorner, w_regionHighCorner, ParticleIterator::ONLY_INNER_AND_BOUNDARY);


    // Corners of symmetric warm region
    double w_s_regionLowCorner[3], w_s_regionHighCorner[3];
    double ymin_symm = _boxLength[1] - _warm_range.ymax;
    double ymax_symm = _boxLength[1] - _warm_range.ymin;
    // if linked cell in the region of interest
    for (unsigned short d = 0; d < 3; d++) {
        w_s_regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
        w_s_regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
    }
    // ensure that we do not iterate over things outside of the container.
    w_s_regionLowCorner[0] = std::max(_cold_range.xmin, w_s_regionLowCorner[0]);
    w_s_regionLowCorner[1] = std::max(ymin_symm, w_s_regionLowCorner[1]);
    w_s_regionLowCorner[2] = std::max(_cold_range.zmin, w_s_regionLowCorner[2]);
    w_s_regionHighCorner[0] = std::min(_cold_range.xmax, w_s_regionHighCorner[0]);
    w_s_regionHighCorner[1] = std::min(ymax_symm, w_s_regionHighCorner[1]);
    w_s_regionHighCorner[2] = std::min(_cold_range.zmax, w_s_regionHighCorner[2]);
    auto begin_w_s = particleContainer->regionIterator(w_s_regionLowCorner, w_s_regionHighCorner, ParticleIterator::ONLY_INNER_AND_BOUNDARY);


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
        if ( v_dummy_abs > v_w_abs.local[cid] ) {
            v_w_abs.local[cid] = v_dummy_abs;
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
        if ( v_dummy_abs < v_c_abs.local[cid] ) {
            v_c_abs.local[cid] = v_dummy_abs;
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
            if ( v_dummy_abs < v_c_abs.local[cid] ) {
                v_c_abs.local[cid] = v_dummy_abs;
                cold_mol_ptr[cid] = &(*it);
            }
        }
    }

#ifdef ENABLE_MPI

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allreduce(v_c_abs.local.data(), v_c_abs.global.data(), _numComp, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(v_w_abs.local.data(), v_w_abs.global.data(), _numComp, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

#else

    v_c_abs.global = v_c_abs.local;
    v_w_abs.global = v_w_abs.local;

#endif

#ifndef NDEBUG
    // find rank with coldest/warmest molecule and read velocities
    for (uint32_t cid = 0; cid < _numComp; cid++) {
        std::cout << "[VelocityExchange] component " << cid << " ; rank " << domainDecomp->getRank() << " ; v_c_abs.local = " << v_c_abs.local[cid] << std::endl;
        std::cout << "[VelocityExchange] component " << cid << " ; rank " << domainDecomp->getRank() << " ; v_w_abs.local = " << v_w_abs.local[cid] << std::endl;
        std::cout << "[VelocityExchange] component " << cid << " ; rank " << domainDecomp->getRank() << " ; v_c_abs.global = " << v_c_abs.global[cid] << std::endl;
        std::cout << "[VelocityExchange] component " << cid << " ; rank " << domainDecomp->getRank() << " ; v_w_abs.global = " << v_w_abs.global[cid] << std::endl;
    }
#endif

    // Seek the molecule with the highest abs. velocity and store its ID and direction-wise velos
    for (uint32_t cid = 0; cid < _numComp; cid++) {
        // Warm region
        if ( (abs(v_c_abs.local[cid]-v_c_abs.global[cid]) < 1e-9) && (v_c_abs.local[cid] < 1000.0) ) {
            mol_c.local[cid] = cold_mol_ptr[cid]->getID();
            for (unsigned short d = 0; d < 3; d++) {
                v_c.local[d].at(cid) = cold_mol_ptr[cid]->v(d);
                D_c.local[d].at(cid) = cold_mol_ptr[cid]->D(d);
            }
        }

        // Cold region
        if ( (abs(v_w_abs.local[cid]-v_w_abs.global[cid]) < 1e-9) && (v_w_abs.local[cid] > 0.0) ) {
            mol_w.local[cid] = warm_mol_ptr[cid]->getID();
            for (unsigned short d = 0; d < 3; d++) {
                v_w.local[d].at(cid) = warm_mol_ptr[cid]->v(d);
                D_w.local[d].at(cid) = warm_mol_ptr[cid]->D(d);
            }
        }
    }


// distribute IDs and velocities to all ranks using MPI_MAX as only one rank contains real information
#ifdef ENABLE_MPI
    MPI_Allreduce(mol_w.local.data(), mol_w.global.data(), _numComp, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(mol_c.local.data(), mol_c.global.data(), _numComp, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    for (unsigned short d = 0; d < 3; d++) {
        MPI_Allreduce(v_c.local[d].data(), v_c.global[d].data(), _numComp, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(v_w.local[d].data(), v_w.global[d].data(), _numComp, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(D_c.local[d].data(), D_c.global[d].data(), _numComp, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(D_w.local[d].data(), D_w.global[d].data(), _numComp, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
#else
    mol_w.global = mol_w.local;
    mol_c.global = mol_c.local;

    for (unsigned short d = 0; d < 3; d++) {
        v_c.local[d] = v_c.global[d];
        v_w.local[d] = v_w.global[d];
        D_c.local[d] = D_c.global[d];
        D_w.local[d] = D_w.global[d];
    }
#endif

    // Assign velocities of coldest particle in warm region to warmest particle in cold region
    for (auto it = begin_c; it.isValid(); ++it) {
        uint32_t cid = it->componentid();
        if ( it->getID() == mol_w.global[cid] ) {
            for (unsigned short d = 0; d < 3; d++) {
                it->setv(d, v_c.global[d].at(cid));
                it->setD(d, D_c.global[d].at(cid));
            }
        }
    }

    // Assign velocities of warmest particle in cold region to coldest particle in warm region
    for (auto it = begin_w; it.isValid(); ++it) {
        uint32_t cid = it->componentid();
        if ( it->getID() == mol_c.global[cid] ) {
            for (unsigned short d = 0; d < 3; d++) {
                it->setv(d, v_w.global[d].at(cid));
                it->setD(d, D_w.global[d].at(cid));
            }
        }
    }

    if ( _symmetry ) {
        for (auto it = begin_w_s; it.isValid(); ++it) {
            uint32_t cid = it->componentid();
            if ( it->getID() == mol_c.global[cid] ) {
                for (unsigned short d = 0; d < 3; d++) {
                    it->setv(d, v_w.global[d].at(cid));
                    it->setD(d, D_w.global[d].at(cid));
                }
            }
        }
    }

    for (uint32_t cid = 0; cid < _numComp; cid++) {
        global_log->info() << "[VelocityExchange] flipped velocities of molecules "
                            << mol_c.global[cid] << " (in warm region, v_new = " << v_w_abs.global[cid] << ") and "
                            << mol_w.global[cid] << " (in cold region v_new = " << v_c_abs.global[cid] << ") of component " << cid+1 << std::endl;
    }
}
