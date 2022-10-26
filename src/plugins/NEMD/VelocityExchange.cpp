
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

void VelocityExchange::init(ParticleContainer* /* particleContainer */, DomainDecompBase* /* domainDecomp */, Domain* /* domain */) {
    _numComp = global_simulation->getEnsemble()->getComponents()->size();
}

void VelocityExchange::readXML(XMLfileUnits& xmlconfig) {
    // Timestep control
    xmlconfig.getNodeValue("control/start", _control.start);
    xmlconfig.getNodeValue("control/frequency", _control.freq);
    xmlconfig.getNodeValue("control/stop", _control.stop);
    global_log->info() << "[VelocityExchange] takes place start:freq:stop = " << _control.start << ":" << _control.freq << ":" << _control.stop << std::endl;

    // Needed here as readXML() is run before init()
    Domain* domain = global_simulation->getDomain();
    for (unsigned short d = 0; d < 3; d++) {
        _boxLength[d] = domain->getGlobalLength(d);
    }
    
    // set default in x and z direction to box size
    // default case is symmetric with region sizes arbitrarily set to 5.0
    // cold region in middle of simulation box
    _cold_region.max[0] = _boxLength[0];
    _cold_region.max[1] = _boxLength[1];
    _cold_region.max[2] = _boxLength[2];
    _warm_region.max[0] = _boxLength[0];
    _warm_region.max[1] = _boxLength[1];
    _warm_region.max[2] = _boxLength[2];
    _symmetry = true;

    std::string strVal[6];  // helper variable for accepting "box" as input

    xmlconfig.getNodeValue("coldregion/xmin", _cold_region.min[0]);
    xmlconfig.getNodeValue("coldregion/ymin", _cold_region.min[1]);
    xmlconfig.getNodeValue("coldregion/zmin", _cold_region.min[2]);
    xmlconfig.getNodeValue("coldregion/xmax", strVal[0]);
    xmlconfig.getNodeValue("coldregion/ymax", strVal[1]);
    xmlconfig.getNodeValue("coldregion/zmax", strVal[2]);

    xmlconfig.getNodeValue("warmregion/symmetric", _symmetry);
    xmlconfig.getNodeValue("warmregion/xmin", _warm_region.min[0]);
    xmlconfig.getNodeValue("warmregion/ymin", _warm_region.min[1]);
    xmlconfig.getNodeValue("warmregion/zmin", _warm_region.min[2]);
    xmlconfig.getNodeValue("warmregion/xmax", strVal[3]);
    xmlconfig.getNodeValue("warmregion/ymax", strVal[4]);
    xmlconfig.getNodeValue("warmregion/zmax", strVal[5]);

	// accept "box" as input
    _cold_region.max[0] = (strVal[0] == "box") ? _boxLength[0] : atof(strVal[0].c_str());
    _cold_region.max[1] = (strVal[1] == "box") ? _boxLength[1] : atof(strVal[1].c_str());
    _cold_region.max[2] = (strVal[2] == "box") ? _boxLength[2] : atof(strVal[2].c_str());
    _warm_region.max[0] = (strVal[3] == "box") ? _boxLength[0] : atof(strVal[3].c_str());
    _warm_region.max[1] = (strVal[4] == "box") ? _boxLength[1] : atof(strVal[4].c_str());
    _warm_region.max[2] = (strVal[5] == "box") ? _boxLength[2] : atof(strVal[5].c_str());

    global_log->info() << "[VelocityExchange] Cold region:"
                        << " x = " << _cold_region.min[0] << " - " << _cold_region.max[0] << " ;"
                        << " y = " << _cold_region.min[1] << " - " << _cold_region.max[1] << " ;"
                        << " z = " << _cold_region.min[2] << " - " << _cold_region.max[2] << endl;

    global_log->info() << "[VelocityExchange] Warm region" << ((_symmetry) ? " (left)" : "") << ":"
                        << " x = " << _warm_region.min[0] << " - " << _warm_region.max[0] << " ;"
                        << " y = " << _warm_region.min[1] << " - " << _warm_region.max[1] << " ;"
                        << " z = " << _warm_region.min[2] << " - " << _warm_region.max[2] << endl;

    if (_symmetry) {
        global_log->info() << "[VelocityExchange] Warm region (right):"
                        << " x = " << _warm_region.min[0] << " - " << _warm_region.max[0] << " ;"
                        << " y = " << _boxLength[1]-_warm_region.max[1] << " - " << _boxLength[1]-_warm_region.min[1] << " ;"
                        << " z = " << _warm_region.min[2] << " - " << _warm_region.max[2] << endl;
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
    coldRegion_lowCorner[0] = std::max(_cold_region.min[0], coldRegion_lowCorner[0]);
    coldRegion_lowCorner[1] = std::max(_cold_region.min[1], coldRegion_lowCorner[1]);
    coldRegion_lowCorner[2] = std::max(_cold_region.min[2], coldRegion_lowCorner[2]);
    coldRegion_highCorner[0] = std::min(_cold_region.max[0], coldRegion_highCorner[0]);
    coldRegion_highCorner[1] = std::min(_cold_region.max[1], coldRegion_highCorner[1]);
    coldRegion_highCorner[2] = std::min(_cold_region.max[2], coldRegion_highCorner[2]);
    const auto begin_c = particleContainer->regionIterator(coldRegion_lowCorner, coldRegion_highCorner, ParticleIterator::ONLY_INNER_AND_BOUNDARY);


    // Corners of warm region
    double warmRegion_lowCorner[3], warmRegion_highCorner[3];
    // if linked cell in the region of interest
    for (unsigned short d = 0; d < 3; d++) {
        warmRegion_lowCorner[d] = particleContainer->getBoundingBoxMin(d);
        warmRegion_highCorner[d] = particleContainer->getBoundingBoxMax(d);
    }
    // ensure that we do not iterate over things outside of the container.
    warmRegion_lowCorner[0] = std::max(_warm_region.min[0], warmRegion_lowCorner[0]);
    warmRegion_lowCorner[1] = std::max(_warm_region.min[1], warmRegion_lowCorner[1]);
    warmRegion_lowCorner[2] = std::max(_warm_region.min[2], warmRegion_lowCorner[2]);
    warmRegion_highCorner[0] = std::min(_warm_region.max[0], warmRegion_highCorner[0]);
    warmRegion_highCorner[1] = std::min(_warm_region.max[1], warmRegion_highCorner[1]);
    warmRegion_highCorner[2] = std::min(_warm_region.max[2], warmRegion_highCorner[2]);
    const auto begin_w = particleContainer->regionIterator(warmRegion_lowCorner, warmRegion_highCorner, ParticleIterator::ONLY_INNER_AND_BOUNDARY);


    // Corners of symmetric warm region
    double warmRegion_sym_lowCorner[3], warmRegion_sym_highCorner[3];
    double ymin_symm = _boxLength[1] - _warm_region.max[1];
    double ymax_symm = _boxLength[1] - _warm_region.min[1];
    // if linked cell in the region of interest
    for (unsigned short d = 0; d < 3; d++) {
        warmRegion_sym_lowCorner[d] = particleContainer->getBoundingBoxMin(d);
        warmRegion_sym_highCorner[d] = particleContainer->getBoundingBoxMax(d);
    }
    // ensure that we do not iterate over things outside of the container.
    warmRegion_sym_lowCorner[0] = std::max(_cold_region.min[0], warmRegion_sym_lowCorner[0]);
    warmRegion_sym_lowCorner[1] = std::max(ymin_symm, warmRegion_sym_lowCorner[1]);
    warmRegion_sym_lowCorner[2] = std::max(_cold_region.min[2], warmRegion_sym_lowCorner[2]);
    warmRegion_sym_highCorner[0] = std::min(_cold_region.max[0], warmRegion_sym_highCorner[0]);
    warmRegion_sym_highCorner[1] = std::min(ymax_symm, warmRegion_sym_highCorner[1]);
    warmRegion_sym_highCorner[2] = std::min(_cold_region.max[2], warmRegion_sym_highCorner[2]);
    const auto begin_w_s = particleContainer->regionIterator(warmRegion_sym_lowCorner, warmRegion_sym_highCorner, ParticleIterator::ONLY_INNER_AND_BOUNDARY);


    std::vector<Molecule*> warm_mol_ptr(_numComp);
    std::vector<Molecule*> cold_mol_ptr(_numComp);

    // find warmest mol in cold region
    findExtremeMols(begin_c, true, velocity_abs_warmP, warm_mol_ptr);

    // find coldest mol in warm region
    findExtremeMols(begin_w, false, velocity_abs_coldP, cold_mol_ptr);

    // if symmetric, look in both warm regions
    if ( _symmetry ) {
        findExtremeMols(begin_w_s, false, velocity_abs_coldP, cold_mol_ptr);
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
    // Needed to distribute the global values
    for (uint32_t cid = 0; cid < _numComp; cid++) {
        // Warm region
        if ( (abs(velocity_abs_coldP.local[cid]-velocity_abs_coldP.global[cid]) < 1e-9) && (velocity_abs_coldP.local[cid] < 1000.0) ) {
            // Only rank containing coldest particle
            molID_coldP.local[cid] = cold_mol_ptr[cid]->getID();
            for (unsigned short d = 0; d < 3; d++) {
                velocity_coldP.local[d].at(cid) = cold_mol_ptr[cid]->v(d);
                rotVelo_coldP.local[d].at(cid) = cold_mol_ptr[cid]->D(d);
            }
        }

        // Cold region
        if ( (abs(velocity_abs_warmP.local[cid]-velocity_abs_warmP.global[cid]) < 1e-9) && (velocity_abs_warmP.local[cid] > 0.0) ) {
            // Only rank containing warmest particle
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
    assignVelocities(begin_c, molID_warmP, velocity_coldP, rotVelo_coldP);

    // Assign velocities of warmest particle in cold region to coldest particle in warm region
    assignVelocities(begin_w, molID_coldP, velocity_warmP, rotVelo_warmP);

    if ( _symmetry ) {
        assignVelocities(begin_w_s, molID_coldP, velocity_warmP, rotVelo_warmP);
    }

    for (uint32_t cid = 0; cid < _numComp; cid++) {
        global_log->info() << domainDecomp->getRank() << " [VelocityExchange] flipped velocities of molecules "
                            << molID_coldP.global[cid] << " (in warm region, v_new = " << velocity_warmP.global[0].at(cid) << " " << velocity_warmP.global[1].at(cid) << " " << velocity_warmP.global[2].at(cid) << ") and "
                            << molID_warmP.global[cid] << " (in cold region v_new = " << velocity_coldP.global[0].at(cid) << " " << velocity_coldP.global[1].at(cid) << " " << velocity_coldP.global[2].at(cid) << ") of component " << cid+1 << std::endl;
    }
}

void VelocityExchange::findExtremeMols(const RegionParticleIterator& begin_iterator, const bool flgColdRegion,
                                       CommVar<std::vector<double>>& velocity_abs, std::vector<Molecule*>& mol_ptr) {
    for (auto it = begin_iterator; it.isValid(); ++it) {
        uint32_t cid = it->componentid();
        const double v_dummy[3] = {it->v(0), it->v(1), it->v(2)};
        const double v_dummy_abs = v_dummy[0]*v_dummy[0]+v_dummy[1]*v_dummy[1]+v_dummy[2]*v_dummy[2];
        double factor = 0.0f;
        if (flgColdRegion) {
            factor = v_dummy_abs - velocity_abs.local[cid];  // In cold region, warmest particle is searched for
        } else {
            factor = velocity_abs.local[cid] - v_dummy_abs;  // In warm region, coldest particle is searched for
        }
        if ( factor > 0.0f ) {
            velocity_abs.local[cid] = v_dummy_abs;
            mol_ptr[cid] = &(*it);
        }
    }
}

void VelocityExchange::assignVelocities(const RegionParticleIterator& begin_iterator, const  CommVar<std::vector<unsigned long>>& molID,
                                        const CommVar<std::array<std::vector<double>, 3>>& velocity, const CommVar<std::array<std::vector<double>, 3>>& rotVelo) {
    for (auto it = begin_iterator; it.isValid(); ++it) {
            uint32_t cid = it->componentid();
            if ( it->getID() == molID.global[cid] ) {
                for (unsigned short d = 0; d < 3; d++) {
                    it->setv(d, velocity.global[d].at(cid));
                    it->setD(d, rotVelo.global[d].at(cid));
                }
            }
    }
}
