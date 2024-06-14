//
// Created by alex on 09.05.23.
//

#include "AdResSForceAdapter.h"
#include "molecules/potforce.h"
#include "plugins/AdResS/AdResS.h"

AdResSForceAdapter::AdResSForceAdapter(Resolution::Handler& resolutionHandler) : _resolutionHandler(resolutionHandler),
																				 _mesoValues() {
    const int numThreads = mardyn_get_max_threads();
    Log::global_log->info() << "[AdResSForceAdapter]: allocate data for " << numThreads << " threads." << std::endl;

    _threadData.resize(numThreads);
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
    {
        auto* myown = new PP2PFAThreadData();
        const int myid = mardyn_get_thread_num();
        _threadData[myid] = myown;
    } // end pragma omp parallel
}

AdResSForceAdapter::~AdResSForceAdapter() noexcept {
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
    {
        const int myid = mardyn_get_thread_num();
        delete _threadData[myid];
    } // end pragma omp parallel
}

void AdResSForceAdapter::init() {
    Domain* domain = _simulation.getDomain();
	#if defined(_OPENMP)
	#pragma omp parallel
	#endif
	{
		const int myid = mardyn_get_thread_num();
		_threadData[myid]->initComp2Param(domain->getComp2Params());
		_threadData[myid]->clear();
	}
}

void AdResSForceAdapter::finish() {
    for(auto tLocal : _threadData) {
        _mesoValues._virial += tLocal->_virial;
		_mesoValues._upot6LJ += tLocal->_upot6LJ;
		_mesoValues._upotXpoles += tLocal->_upotXpoles;
		_mesoValues._myRF += tLocal->_myRF;
    }

	_mesoValues.setInDomain(_simulation.getDomain());
	_mesoValues.clear();
}

double
AdResSForceAdapter::processPair(Molecule& particle1, Molecule& particle2, double distanceVector[3], PairType pairType, double dd, bool calculateLJ) {
	auto& region = _resolutionHandler.getRegions()[0];
	auto& compMap = _resolutionHandler.getCompResMap();
    bool hasNoHybrid = compMap[particle1.componentid()] != Resolution::Hybrid && compMap[particle2.componentid()] != Resolution::Hybrid;
    return processPair(particle1, particle2, distanceVector, pairType, dd, calculateLJ, _resolutionHandler.getCompResMap(), hasNoHybrid, region);
}

double AdResSForceAdapter::processPair(Molecule &molecule1, Molecule &molecule2, double * distanceVector,
                                       PairType pairType,
                                       double dd, bool calculateLJ,
                                       const std::vector<Resolution::ResolutionType> &compResMap, bool noHybrid,
									   const Resolution::FPRegion &region) {
    const int tid = mardyn_get_thread_num();
    PP2PFAThreadData &my_threadData = *_threadData[tid];
    ParaStrm& params = (* my_threadData._comp2Param)(molecule1.componentid(), molecule2.componentid());
    ParaStrm& paramsInv = (* my_threadData._comp2Param)(molecule2.componentid(), molecule1.componentid());
    params.reset_read();
    paramsInv.reset_read();
    switch (pairType) {

        double dummy1, dummy2, dummy3, dummy4[3], Virial3[3];

        case MOLECULE_MOLECULE :
            potForce(molecule1, molecule2, params, paramsInv, distanceVector, my_threadData._upot6LJ, my_threadData._upotXpoles,
                     my_threadData._myRF, Virial3, calculateLJ, noHybrid, compResMap, region);
            my_threadData._virial += 2*(Virial3[0]+Virial3[1]+Virial3[2]);
            return my_threadData._upot6LJ + my_threadData._upotXpoles;
        case MOLECULE_HALOMOLECULE :
            potForce(molecule1, molecule2, params, paramsInv, distanceVector, dummy1, dummy2, dummy3, dummy4, calculateLJ, noHybrid,
                     compResMap, region);
            return 0.0;
        case MOLECULE_MOLECULE_FLUID :
            dummy1 = 0.0; // 6*U_LJ
            dummy2 = 0.0; // U_polarity
            dummy3 = 0.0; // U_dipole_reaction_field

            fluidPot(molecule1, molecule2, params, paramsInv, distanceVector, dummy1, dummy2, dummy3, calculateLJ, noHybrid, compResMap, region);
            return dummy1 / 6.0 + dummy2 + dummy3;
        default:
            Simulation::exit(670); // not implemented
    }
    return 0.0;
}

void
AdResSForceAdapter::potForce(Molecule &mi, Molecule &mj, ParaStrm &params, ParaStrm &paramInv, double * drm, double &Upot6LJ,
                             double &UpotXpoles,
                             double &MyRF, double Virial[3], bool calculateLJ, bool noHybrid,
                             const std::vector<Resolution::ResolutionType> &compResMap, const Resolution::FPRegion &region) {
    if(noHybrid) {
        PotForce(mi, mj, params, drm, Upot6LJ, UpotXpoles, MyRF, Virial, calculateLJ);
        return;
    }

    bool isHybridI, isHybridJ;
    isHybridI = compResMap[mi.componentid()] == Resolution::Hybrid;
    isHybridJ = compResMap[mj.componentid()] == Resolution::Hybrid;

    if(isHybridI && isHybridJ) {
        potForceFullHybrid(mi, mj, params, drm, Upot6LJ, UpotXpoles, MyRF, Virial, calculateLJ, region);
        return;
    }
    if(isHybridI) {
        potForceSingleHybrid(mi, mj, params, drm, Upot6LJ, UpotXpoles, MyRF, Virial, calculateLJ, region, compResMap[mj.componentid()]);
        return;
    }
    if(isHybridJ) {
        potForceSingleHybrid(mj, mi, paramInv, drm, Upot6LJ, UpotXpoles, MyRF, Virial, calculateLJ, region, compResMap[mi.componentid()]);
        return;
    }

    // we should never reach this point
    // only get to here if some particles component is not set correctly, can happen during initialization
    // fall back regular PotForce
    PotForce(mi, mj, params, drm, Upot6LJ, UpotXpoles, MyRF, Virial, calculateLJ);
    Log::global_log->debug() << "[AdResS] AdResSForceAdapter::potForce called with 2 non hybrid molecules" << std::endl;
}

void AdResSForceAdapter::fluidPot(Molecule &mi, Molecule &mj, ParaStrm &params, ParaStrm &paramInv, double *drm, double &Upot6LJ,
                                  double &UpotXpoles, double &MyRF, bool calculateLJ, bool noHybrid,
                                  const std::vector<Resolution::ResolutionType> &compResMap, const Resolution::FPRegion &region) {
    if(noHybrid) {
        FluidPot(mi, mj, params, drm, Upot6LJ, UpotXpoles, MyRF, calculateLJ);
        return;
    }

    bool isHybridI, isHybridJ;
    isHybridI = compResMap[mi.componentid()] == Resolution::Hybrid;
    isHybridJ = compResMap[mj.componentid()] == Resolution::Hybrid;

    if(isHybridI && isHybridJ) {
        fluidPotFullHybrid(mi, mj, params, drm, Upot6LJ, UpotXpoles, MyRF, calculateLJ, region);
        return;
    }
    if(isHybridI) {
        fluidPotSingleHybrid(mi, mj, params, drm, Upot6LJ, UpotXpoles, MyRF, calculateLJ, region, compResMap[mj.componentid()]);
        return;
    }
    if(isHybridJ) {
        fluidPotSingleHybrid(mj, mi, paramInv, drm, Upot6LJ, UpotXpoles, MyRF, calculateLJ, region, compResMap[mi.componentid()]);
        return;
    }

    // we should never reach this point
    //Simulation::exit(671);
    Log::global_log->warning() << "[AdResS] AdResSForceAdapter::fluidPot called with 2 non hybrid molecules" << std::endl;
}

void
AdResSForceAdapter::potForceFullHybrid(Molecule &mi, Molecule &mj, ParaStrm &params, double * drm,
                                       double &Upot6LJ,
                                       double &UpotXpoles, double &MyRF, double Virial[3], bool calculateLJ,
									   const Resolution::FPRegion &region) {
    auto& components = *_simulation.getEnsemble()->getComponents();
    // F_a,b = w(r_a)w(r_b)F_FP(a, b) + (1-w(r_a))(1-w(r_b))F_CG(a, b)
    // the first X sites of the component k with mass 0 are part of the CG model
    // only let CG sites interact with other CG sites etc...

    double wi, wj;
    double f[3];
    double u;
    double drs[3], dr2; // site distance vector & length^2
    Virial[0]=0.;
    Virial[1]=0.;
    Virial[2]=0.;
    // LJ centers
    // no LJ interaction between solid atoms of the same component
    wi = AdResS::weight(mi.r_arr(), region);
    wj = AdResS::weight(mj.r_arr(), region);

    //handle LJ sites
    {
        const unsigned int nCG_i = components[mi.componentid() + 1].numLJcenters();
        const unsigned int nCG_j = components[mj.componentid() + 1].numLJcenters();
        const unsigned int nc1 = mi.numLJcenters();
        const unsigned int nc2 = mj.numLJcenters();
        for (unsigned int si = 0; si < nc1; ++si) {
            const std::array<double,3> dii = mi.ljcenter_d_abs(si);
            bool isCGi = si < nCG_i;
            for (unsigned int sj = 0; sj < nc2; ++sj) {
                bool isCGj = sj < nCG_j;
                //both sites must be CG or FP but not mixed
                if(isCGi ^ isCGj) {
                    double tmp; params >> tmp; params >> tmp; params >> tmp;
                    continue;
                }

                const std::array<double,3> djj = mj.ljcenter_d_abs(sj);
                SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
                double eps24;
                params >> eps24;
                double sig2;
                params >> sig2;
                double shift6;
                params >> shift6; // must be 0.0 for full LJ
                if (calculateLJ) {
                    PotForceLJ(drs, dr2, eps24, sig2, f, u);
                    u += shift6;

                    //if mass 0 -> weight inv; if mass > 0 weight
                    if(isCGi) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj; }
                    else { for (double &d: f) d *= wi * wj; u *= wi * wj; }

                    mi.Fljcenteradd(si, f);
                    mj.Fljcentersub(sj, f);
                    Upot6LJ += u;
                    for (unsigned short d = 0; d < 3; ++d)
                        Virial[d] += 0.5*drm[d] * f[d];
                }
            }
        }
    }


    double m1[3], m2[3]; // angular momenta
    const unsigned nCG_Ci = components[mi.componentid() + 1].numCharges();
    const unsigned nCG_Cj = components[mj.componentid() + 1].numCharges();
    const unsigned nCG_Qi = components[mi.componentid() + 1].numQuadrupoles();
    const unsigned nCG_Qj = components[mj.componentid() + 1].numQuadrupoles();
    const unsigned nCG_Di = components[mi.componentid() + 1].numDipoles();
    const unsigned nCG_Dj = components[mj.componentid() + 1].numDipoles();
    const unsigned ne1 = mi.numCharges();
    const unsigned ne2 = mj.numCharges();
    const unsigned int nq1 = mi.numQuadrupoles();
    const unsigned int nq2 = mj.numQuadrupoles();
    const unsigned int nd1 = mi.numDipoles();
    const unsigned int nd2 = mj.numDipoles();
    for (unsigned si = 0; si < ne1; si++) {
        bool isCGi = si < nCG_Ci;
        const std::array<double,3> dii = mi.charge_d_abs(si);
        // Charge-Charge
        for (unsigned sj = 0; sj < ne2; sj++) {
            bool isCGj = sj < nCG_Cj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.charge_d_abs(sj);
            double q1q2per4pie0; // 4pie0 = 1 in reduced units
            params >> q1q2per4pie0;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            PotForce2Charge(drs, dr2, q1q2per4pie0, f, u);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCGi) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj; }
            else { for (double &d: f) d *= wi * wj; u *= wi * wj; }
            mi.Fchargeadd(si, f);
            mj.Fchargesub(sj, f);

            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; d++)
                Virial[d] += 0.5*drm[d] * f[d];
        }
        // Charge-Quadrupole
        for (unsigned sj = 0; sj < nq2; sj++) {
            bool isCGj = sj < nCG_Qj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
            double qQ05per4pie0; // 4pie0 = 1 in reduced units
            params >> qQ05per4pie0;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.quadrupole_e(sj);
            PotForceChargeQuadrupole(drs, dr2, ejj.data(), qQ05per4pie0, f, m2, u);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCGi) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
            else { for (double &d: f) d *= wi * wj; u *= wi * wj;}
            mi.Fchargeadd(si, f);
            mj.Fquadrupolesub(sj, f);
            mj.Madd(m2);

            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; d++)
                Virial[d] += 0.5*drm[d] * f[d];
        }
        // Charge-Dipole
        for (unsigned sj = 0; sj < nd2; sj++) {
            bool isCGj = sj < nCG_Dj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.dipole_d_abs(sj);
            double minusqmyper4pie0;
            params >> minusqmyper4pie0;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.dipole_e(sj);
            PotForceChargeDipole(drs, dr2, ejj.data(), minusqmyper4pie0, f, m2, u);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCGi) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
            else { for (double &d: f) d *= wi * wj; u *= wi * wj;}
            mi.Fchargeadd(si, f);
            mj.Fdipolesub(sj, f);
            mj.Madd(m2);

            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; d++)
                Virial[d] += 0.5*drm[d] * f[d];
        }
    }

    for (unsigned int si = 0; si < nq1; ++si) {
        const std::array<double,3> dii = mi.quadrupole_d_abs(si);
        const std::array<double,3> eii = mi.quadrupole_e(si);
        bool isCGi = si < nCG_Qi;
        // Quadrupole-Charge
        for (unsigned sj = 0; sj < ne2; sj++) {
            bool isCGj = sj < nCG_Cj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.charge_d_abs(sj);
            double qQ05per4pie0; // 4pie0 = 1 in reduced units
            params >> qQ05per4pie0;
            minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            PotForceChargeQuadrupole(drs, dr2, eii.data(), qQ05per4pie0, f, m1, u);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCGi) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
            else { for (double &d: f) d *= wi * wj; u *= wi * wj;}
            mi.Fquadrupolesub(si, f);
            mj.Fchargeadd(sj, f);
            mi.Madd(m1);

            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; d++)
                Virial[d] -= 0.5*drm[d] * f[d];
        }
        // Quadrupole-Quadrupole -------------------
        for (unsigned int sj = 0; sj < nq2; ++sj) {
            bool isCGj = sj < nCG_Qj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp;
                continue;
            }
            //double drs[3];
            const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
            double q2075;
            params >> q2075;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.quadrupole_e(sj);
            PotForce2Quadrupole(drs, dr2, eii.data(), ejj.data(), q2075, f, m1, m2, u);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCGi) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
            else { for (double &d: f) d *= wi * wj; u *= wi * wj;}
            mi.Fquadrupoleadd(si, f);
            mj.Fquadrupolesub(sj, f);
            mi.Madd(m1);
            mj.Madd(m2);

            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; ++d)
                Virial[d] += 0.5*drm[d] * f[d];
        }
        // Quadrupole-Dipole -----------------------
        for (unsigned int sj = 0; sj < nd2; ++sj) {
            bool isCGj = sj < nCG_Dj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp;
                continue;
            }
            //double drs[3];
            const std::array<double,3> djj = mj.dipole_d_abs(sj);
            double qmy15;
            params >> qmy15;
            minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.dipole_e(sj);
            PotForceDiQuadrupole(drs, dr2, ejj.data(), eii.data(), qmy15, f, m2, m1, u);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCGi) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
            else { for (double &d: f) d *= wi * wj; u *= wi * wj;}
            mi.Fquadrupolesub(si, f);
            mj.Fdipoleadd(sj, f);
            mi.Madd(m1);
            mj.Madd(m2);
            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; d++)
                Virial[d] -= 0.5*drm[d] * f[d];
        }
    }

    for (unsigned int si = 0; si < nd1; ++si) {
        bool isCGi = si < nCG_Di;
        const std::array<double,3> dii = mi.dipole_d_abs(si);
        const std::array<double,3> eii = mi.dipole_e(si);
        // Dipole-Charge
        for (unsigned sj = 0; sj < ne2; sj++) {
            bool isCGj = sj < nCG_Cj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.charge_d_abs(sj);
            double minusqmyper4pie0;
            params >> minusqmyper4pie0;
            minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            PotForceChargeDipole(drs, dr2, eii.data(), minusqmyper4pie0, f, m1, u);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCGi) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
            else { for (double &d: f) d *= wi * wj; u *= wi * wj;}
            mi.Fdipolesub(si, f);
            mj.Fchargeadd(sj, f);
            mi.Madd(m1);

            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; d++)
                Virial[d] -= 0.5*drm[d] * f[d];
        }
        // Dipole-Quadrupole -----------------------
        for (unsigned int sj = 0; sj < nq2; ++sj) {
            bool isCGj = sj < nCG_Qj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp;
                continue;
            }
            //double drs[3];
            const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
            double myq15;
            params >> myq15;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.quadrupole_e(sj);
            PotForceDiQuadrupole(drs, dr2, eii.data(), ejj.data(), myq15, f, m1, m2, u);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCGi) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
            else { for (double &d: f) d *= wi * wj; u *= wi * wj;}
            mi.Fdipoleadd(si, f);
            mj.Fquadrupolesub(sj, f);
            mi.Madd(m1);
            mj.Madd(m2);
            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; ++d)
                Virial[d] += 0.5*drm[d] * f[d];
        }
        // Dipole-Dipole ---------------------------
        for (unsigned int sj = 0; sj < nd2; ++sj) {
            bool isCGj = sj < nCG_Dj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.dipole_d_abs(sj);
            double my2;
            params >> my2;
            double rffac;
            params >> rffac;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.dipole_e(sj);
            PotForce2Dipole(drs, dr2, eii.data(), ejj.data(), my2, rffac, f, m1, m2, u, MyRF);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCGi) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
            else { for (double &d: f) d *= wi * wj; u *= wi * wj;}
            mi.Fdipoleadd(si, f);
            mj.Fdipolesub(sj, f);
            mi.Madd(m1);
            mj.Madd(m2);
            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; ++d)
                Virial[d] += 0.5*drm[d] * f[d];
        }
    }

    mi.Viadd(Virial);
    mj.Viadd(Virial);

    // check whether all parameters were used
    mardyn_assert(params.eos());
}

void
AdResSForceAdapter::potForceSingleHybrid(Molecule &mi, Molecule &mj, ParaStrm &params, double * drm,
                                         double &Upot6LJ,
                                         double &UpotXpoles, double &MyRF, double Virial[3], bool calculateLJ,
										 const Resolution::FPRegion &region,
										 Resolution::ResolutionType resolutionJ) {
    auto& components = *_simulation.getEnsemble()->getComponents();
    // F_a,b = w(r_a)w(r_b)F_FP(a, b) + (1-w(r_a))(1-w(r_b))F_CG(a, b)
    // the first X sites of the component k with mass 0 are part of the CG model
    // only let CG sites interact with other CG sites etc...

    double wi, wj;
    double f[3];
    double u;
    double drs[3], dr2; // site distance vector & length^2
    Virial[0]=0.;
    Virial[1]=0.;
    Virial[2]=0.;
    // LJ centers
    // no LJ interaction between solid atoms of the same component
    wi = AdResS::weight(mi.r_arr(), region);
    wj = AdResS::weight(mj.r_arr(), region);

    //handle LJ sites
    {
        const unsigned int nCG_LJ = components[mi.componentid() + 1].numLJcenters();
        const unsigned int nc1 = mi.numLJcenters();
        const unsigned int nc2 = mj.numLJcenters();
        for (unsigned int si = 0; si < nc1; ++si) {
            const std::array<double,3> dii = mi.ljcenter_d_abs(si);
            bool isCGSiteI = si < nCG_LJ;
            for (unsigned int sj = 0; sj < nc2; ++sj) {
                //both sites must be CG or FP but not mixed
                if((resolutionJ == Resolution::CoarseGrain && !isCGSiteI) ||
                   (resolutionJ == Resolution::FullParticle && isCGSiteI)) {
                    double tmp; params >> tmp; params >> tmp; params >> tmp;
                    continue;
                }

                const std::array<double,3> djj = mj.ljcenter_d_abs(sj);
                SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
                double eps24;
                params >> eps24;
                double sig2;
                params >> sig2;
                double shift6;
                params >> shift6; // must be 0.0 for full LJ
                if (calculateLJ) {
                    PotForceLJ(drs, dr2, eps24, sig2, f, u);
                    u += shift6;

                    //if mass 0 -> weight inv; if mass > 0 weight
                    if(isCGSiteI) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
                    else { for (double &d: f) d *= wi * wj; u *= wi * wj;}

                    mi.Fljcenteradd(si, f);
                    mj.Fljcentersub(sj, f);
                    Upot6LJ += u;
                    for (unsigned short d = 0; d < 3; ++d)
                        Virial[d] += 0.5*drm[d] * f[d];
                }
            }
        }
    }


    double m1[3], m2[3]; // angular momenta

    const unsigned nCG_C_i = components[mi.componentid() + 1].numCharges();
    const unsigned nCG_Q_i = components[mi.componentid() + 1].numQuadrupoles();
    const unsigned nCG_D_i = components[mi.componentid() + 1].numDipoles();
    const unsigned ne1 = mi.numCharges();
    const unsigned ne2 = mj.numCharges();
    const unsigned int nq1 = mi.numQuadrupoles();
    const unsigned int nq2 = mj.numQuadrupoles();
    const unsigned int nd1 = mi.numDipoles();
    const unsigned int nd2 = mj.numDipoles();
    for (unsigned si = 0; si < ne1; si++) {
        bool isCG_i = si < nCG_C_i;
        const std::array<double,3> dii = mi.charge_d_abs(si);
        // Charge-Charge
        for (unsigned sj = 0; sj < ne2; sj++) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.charge_d_abs(sj);
            double q1q2per4pie0; // 4pie0 = 1 in reduced units
            params >> q1q2per4pie0;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            PotForce2Charge(drs, dr2, q1q2per4pie0, f, u);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCG_i) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
            else { for (double &d: f) d *= wi * wj; u *= wi * wj;}
            mi.Fchargeadd(si, f);
            mj.Fchargesub(sj, f);

            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; d++)
                Virial[d] += 0.5*drm[d] * f[d];
        }
        // Charge-Quadrupole
        for (unsigned sj = 0; sj < nq2; sj++) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
            double qQ05per4pie0; // 4pie0 = 1 in reduced units
            params >> qQ05per4pie0;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.quadrupole_e(sj);
            PotForceChargeQuadrupole(drs, dr2, ejj.data(), qQ05per4pie0, f, m2, u);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCG_i) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
            else { for (double &d: f) d *= wi * wj; u *= wi * wj;}
            mi.Fchargeadd(si, f);
            mj.Fquadrupolesub(sj, f);
            mj.Madd(m2);

            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; d++)
                Virial[d] += 0.5*drm[d] * f[d];
        }
        // Charge-Dipole
        for (unsigned sj = 0; sj < nd2; sj++) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.dipole_d_abs(sj);
            double minusqmyper4pie0;
            params >> minusqmyper4pie0;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.dipole_e(sj);
            PotForceChargeDipole(drs, dr2, ejj.data(), minusqmyper4pie0, f, m2, u);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCG_i) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
            else { for (double &d: f) d *= wi * wj; u *= wi * wj;}
            mi.Fchargeadd(si, f);
            mj.Fdipolesub(sj, f);
            mj.Madd(m2);

            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; d++)
                Virial[d] += 0.5*drm[d] * f[d];
        }
    }

    for (unsigned int si = 0; si < nq1; ++si) {
        bool isCG_i = si < nCG_Q_i;
        const std::array<double,3> dii = mi.quadrupole_d_abs(si);
        const std::array<double,3> eii = mi.quadrupole_e(si);

        // Quadrupole-Charge
        for (unsigned sj = 0; sj < ne2; sj++) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.charge_d_abs(sj);
            double qQ05per4pie0; // 4pie0 = 1 in reduced units
            params >> qQ05per4pie0;
            minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            PotForceChargeQuadrupole(drs, dr2, eii.data(), qQ05per4pie0, f, m1, u);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCG_i) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
            else { for (double &d: f) d *= wi * wj; u *= wi * wj;}
            mi.Fquadrupolesub(si, f);
            mj.Fchargeadd(sj, f);
            mi.Madd(m1);

            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; d++)
                Virial[d] -= 0.5*drm[d] * f[d];
        }
        // Quadrupole-Quadrupole -------------------
        for (unsigned int sj = 0; sj < nq2; ++sj) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp;
                continue;
            }
            //double drs[3];
            const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
            double q2075;
            params >> q2075;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.quadrupole_e(sj);
            PotForce2Quadrupole(drs, dr2, eii.data(), ejj.data(), q2075, f, m1, m2, u);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCG_i) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
            else { for (double &d: f) d *= wi * wj; u *= wi * wj;}
            mi.Fquadrupoleadd(si, f);
            mj.Fquadrupolesub(sj, f);
            mi.Madd(m1);
            mj.Madd(m2);

            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; ++d)
                Virial[d] += 0.5*drm[d] * f[d];
        }
        // Quadrupole-Dipole -----------------------
        for (unsigned int sj = 0; sj < nd2; ++sj) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp;
                continue;
            }
            //double drs[3];
            const std::array<double,3> djj = mj.dipole_d_abs(sj);
            double qmy15;
            params >> qmy15;
            minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.dipole_e(sj);
            PotForceDiQuadrupole(drs, dr2, ejj.data(), eii.data(), qmy15, f, m2, m1, u);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCG_i) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
            else { for (double &d: f) d *= wi * wj; u *= wi * wj;}
            mi.Fquadrupolesub(si, f);
            mj.Fdipoleadd(sj, f);
            mi.Madd(m1);
            mj.Madd(m2);
            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; d++)
                Virial[d] -= 0.5*drm[d] * f[d];
        }
    }

    for (unsigned int si = 0; si < nd1; ++si) {
        bool isCG_i = si < nCG_D_i;
        const std::array<double,3> dii = mi.dipole_d_abs(si);
        const std::array<double,3> eii = mi.dipole_e(si);
        // Dipole-Charge
        for (unsigned sj = 0; sj < ne2; sj++) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.charge_d_abs(sj);
            double minusqmyper4pie0;
            params >> minusqmyper4pie0;
            minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            PotForceChargeDipole(drs, dr2, eii.data(), minusqmyper4pie0, f, m1, u);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCG_i) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
            else { for (double &d: f) d *= wi * wj; u *= wi * wj;}
            mi.Fdipolesub(si, f);
            mj.Fchargeadd(sj, f);
            mi.Madd(m1);

            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; d++)
                Virial[d] -= 0.5*drm[d] * f[d];
        }
        // Dipole-Quadrupole -----------------------
        for (unsigned int sj = 0; sj < nq2; ++sj) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp;
                continue;
            }
            //double drs[3];
            const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
            double myq15;
            params >> myq15;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.quadrupole_e(sj);
            PotForceDiQuadrupole(drs, dr2, eii.data(), ejj.data(), myq15, f, m1, m2, u);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCG_i) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
            else { for (double &d: f) d *= wi * wj; u *= wi * wj;}
            mi.Fdipoleadd(si, f);
            mj.Fquadrupolesub(sj, f);
            mi.Madd(m1);
            mj.Madd(m2);
            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; ++d)
                Virial[d] += 0.5*drm[d] * f[d];
        }
        // Dipole-Dipole ---------------------------
        for (unsigned int sj = 0; sj < nd2; ++sj) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.dipole_d_abs(sj);
            double my2;
            params >> my2;
            double rffac;
            params >> rffac;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.dipole_e(sj);
            PotForce2Dipole(drs, dr2, eii.data(), ejj.data(), my2, rffac, f, m1, m2, u, MyRF);

            //if mass 0 -> weight inv; if mass > 0 weight
            if(isCG_i) { for (double &d: f) d *= 1 - wi * wj; u *= 1 - wi * wj;}
            else { for (double &d: f) d *= wi * wj; u *= wi * wj;}
            mi.Fdipoleadd(si, f);
            mj.Fdipolesub(sj, f);
            mi.Madd(m1);
            mj.Madd(m2);
            UpotXpoles += u;
            for (unsigned short d = 0; d < 3; ++d)
                Virial[d] += 0.5*drm[d] * f[d];
        }
    }

    mi.Viadd(Virial);
    mj.Viadd(Virial);

    // check whether all parameters were used
    mardyn_assert(params.eos());
}

void AdResSForceAdapter::fluidPotFullHybrid(Molecule &mi, Molecule &mj, ParaStrm &params, double *drm, double &Upot6LJ,
                                            double &UpotXpoles, double &MyRF, bool calculateLJ, const Resolution::FPRegion &region) {
    auto& components = *_simulation.getEnsemble()->getComponents();
    // F_a,b = w(r_a)w(r_b)F_FP(a, b) + (1-w(r_a))(1-w(r_b))F_CG(a, b)
    // the first X sites of the component k with mass 0 are part of the CG model
    // only let CG sites interact with other CG sites etc...

    double wi, wj;
    double f[3];
    double u;
    double drs[3], dr2; // site distance vector & length^2
    // LJ centers
    // no LJ interaction between solid atoms of the same component
    wi = AdResS::weight(mi.r_arr(), region);
    wj = AdResS::weight(mj.r_arr(), region);

    //handle LJ sites
    {
        const unsigned int nCG_i = components[mi.componentid() + 1].numLJcenters();
        const unsigned int nCG_j = components[mj.componentid() + 1].numLJcenters();
        const unsigned int nc1 = mi.numLJcenters();
        const unsigned int nc2 = mj.numLJcenters();
        for (unsigned int si = 0; si < nc1; ++si) {
            const std::array<double,3> dii = mi.ljcenter_d_abs(si);
            bool isCGi = si < nCG_i;
            for (unsigned int sj = 0; sj < nc2; ++sj) {
                bool isCGj = sj < nCG_j;
                //both sites must be CG or FP but not mixed
                if(isCGi ^ isCGj) {
                    double tmp; params >> tmp; params >> tmp; params >> tmp;
                    continue;
                }

                const std::array<double,3> djj = mj.ljcenter_d_abs(sj);
                SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
                double eps24;
                params >> eps24;
                double sig2;
                params >> sig2;
                double shift6;
                params >> shift6; // must be 0.0 for full LJ
                if (calculateLJ) {
                    PotForceLJ(drs, dr2, eps24, sig2, f, u);
                    if(isCGi) { u *= 1 - wi * wj;}
                    else { u *= wi * wj;}
                    u += shift6;
                    Upot6LJ += u;
                }
            }
        }
    }


    double m1[3], m2[3]; // angular momenta
    const unsigned nCG_Ci = components[mi.componentid() + 1].numCharges();
    const unsigned nCG_Cj = components[mj.componentid() + 1].numCharges();
    const unsigned nCG_Qi = components[mi.componentid() + 1].numQuadrupoles();
    const unsigned nCG_Qj = components[mj.componentid() + 1].numQuadrupoles();
    const unsigned nCG_Di = components[mi.componentid() + 1].numDipoles();
    const unsigned nCG_Dj = components[mj.componentid() + 1].numDipoles();
    const unsigned ne1 = mi.numCharges();
    const unsigned ne2 = mj.numCharges();
    const unsigned int nq1 = mi.numQuadrupoles();
    const unsigned int nq2 = mj.numQuadrupoles();
    const unsigned int nd1 = mi.numDipoles();
    const unsigned int nd2 = mj.numDipoles();
    for (unsigned si = 0; si < ne1; si++) {
        bool isCGi = si < nCG_Ci;
        const std::array<double,3> dii = mi.charge_d_abs(si);
        // Charge-Charge
        for (unsigned sj = 0; sj < ne2; sj++) {
            bool isCGj = sj < nCG_Cj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.charge_d_abs(sj);
            double q1q2per4pie0; // 4pie0 = 1 in reduced units
            params >> q1q2per4pie0;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            PotForce2Charge(drs, dr2, q1q2per4pie0, f, u);
            if(isCGi) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
        // Charge-Quadrupole
        for (unsigned sj = 0; sj < nq2; sj++) {
            bool isCGj = sj < nCG_Qj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
            double qQ05per4pie0; // 4pie0 = 1 in reduced units
            params >> qQ05per4pie0;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.quadrupole_e(sj);
            PotForceChargeQuadrupole(drs, dr2, ejj.data(), qQ05per4pie0, f, m2, u);
            if(isCGi) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
        // Charge-Dipole
        for (unsigned sj = 0; sj < nd2; sj++) {
            bool isCGj = sj < nCG_Dj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.dipole_d_abs(sj);
            double minusqmyper4pie0;
            params >> minusqmyper4pie0;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.dipole_e(sj);
            PotForceChargeDipole(drs, dr2, ejj.data(), minusqmyper4pie0, f, m2, u);
            if(isCGi) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
    }

    for (unsigned int si = 0; si < nq1; ++si) {
        const std::array<double,3> dii = mi.quadrupole_d_abs(si);
        const std::array<double,3> eii = mi.quadrupole_e(si);
        bool isCGi = si < nCG_Qi;
        // Quadrupole-Charge
        for (unsigned sj = 0; sj < ne2; sj++) {
            bool isCGj = sj < nCG_Cj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.charge_d_abs(sj);
            double qQ05per4pie0; // 4pie0 = 1 in reduced units
            params >> qQ05per4pie0;
            minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            PotForceChargeQuadrupole(drs, dr2, eii.data(), qQ05per4pie0, f, m1, u);
            if(isCGi) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
        // Quadrupole-Quadrupole -------------------
        for (unsigned int sj = 0; sj < nq2; ++sj) {
            bool isCGj = sj < nCG_Qj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp;
                continue;
            }
            //double drs[3];
            const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
            double q2075;
            params >> q2075;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.quadrupole_e(sj);
            PotForce2Quadrupole(drs, dr2, eii.data(), ejj.data(), q2075, f, m1, m2, u);
            if(isCGi) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
        // Quadrupole-Dipole -----------------------
        for (unsigned int sj = 0; sj < nd2; ++sj) {
            bool isCGj = sj < nCG_Dj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp;
                continue;
            }
            //double drs[3];
            const std::array<double,3> djj = mj.dipole_d_abs(sj);
            double qmy15;
            params >> qmy15;
            minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.dipole_e(sj);
            PotForceDiQuadrupole(drs, dr2, ejj.data(), eii.data(), qmy15, f, m2, m1, u);
            if(isCGi) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
    }

    for (unsigned int si = 0; si < nd1; ++si) {
        bool isCGi = si < nCG_Di;
        const std::array<double,3> dii = mi.dipole_d_abs(si);
        const std::array<double,3> eii = mi.dipole_e(si);
        // Dipole-Charge
        for (unsigned sj = 0; sj < ne2; sj++) {
            bool isCGj = sj < nCG_Cj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.charge_d_abs(sj);
            double minusqmyper4pie0;
            params >> minusqmyper4pie0;
            minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            PotForceChargeDipole(drs, dr2, eii.data(), minusqmyper4pie0, f, m1, u);
            if(isCGi) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
        // Dipole-Quadrupole -----------------------
        for (unsigned int sj = 0; sj < nq2; ++sj) {
            bool isCGj = sj < nCG_Qj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp;
                continue;
            }
            //double drs[3];
            const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
            double myq15;
            params >> myq15;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.quadrupole_e(sj);
            PotForceDiQuadrupole(drs, dr2, eii.data(), ejj.data(), myq15, f, m1, m2, u);
            if(isCGi) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
        // Dipole-Dipole ---------------------------
        for (unsigned int sj = 0; sj < nd2; ++sj) {
            bool isCGj = sj < nCG_Dj;
            if(isCGi ^ isCGj) {
                double tmp; params >> tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.dipole_d_abs(sj);
            double my2;
            params >> my2;
            double rffac;
            params >> rffac;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.dipole_e(sj);
            PotForce2Dipole(drs, dr2, eii.data(), ejj.data(), my2, rffac, f, m1, m2, u, MyRF);
            if(isCGi) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
    }
    // check whether all parameters were used
    mardyn_assert(params.eos());
}

void
AdResSForceAdapter::fluidPotSingleHybrid(Molecule &mi, Molecule &mj, ParaStrm &params, double *drm, double &Upot6LJ,
                                         double &UpotXpoles, double &MyRF, bool calculateLJ, const Resolution::FPRegion &region,
										 Resolution::ResolutionType resolutionJ) {
    auto& components = *_simulation.getEnsemble()->getComponents();
    // F_a,b = w(r_a)w(r_b)F_FP(a, b) + (1-w(r_a))(1-w(r_b))F_CG(a, b)
    // the first X sites of the component k with mass 0 are part of the CG model
    // only let CG sites interact with other CG sites etc...

    double wi, wj;
    double f[3];
    double u;
    double drs[3], dr2; // site distance vector & length^2
    // LJ centers
    // no LJ interaction between solid atoms of the same component
    wi = AdResS::weight(mi.r_arr(), region);
    wj = AdResS::weight(mj.r_arr(), region);

    //handle LJ sites
    {
        const unsigned int nCG_LJ = components[mi.componentid() + 1].numLJcenters();
        const unsigned int nc1 = mi.numLJcenters();
        const unsigned int nc2 = mj.numLJcenters();
        for (unsigned int si = 0; si < nc1; ++si) {
            const std::array<double,3> dii = mi.ljcenter_d_abs(si);
            bool isCGSiteI = si < nCG_LJ;
            for (unsigned int sj = 0; sj < nc2; ++sj) {
                //both sites must be CG or FP but not mixed
                if((resolutionJ == Resolution::CoarseGrain && !isCGSiteI) ||
                   (resolutionJ == Resolution::FullParticle && isCGSiteI)) {
                    double tmp; params >> tmp; params >> tmp; params >> tmp;
                    continue;
                }

                const std::array<double,3> djj = mj.ljcenter_d_abs(sj);
                SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
                double eps24;
                params >> eps24;
                double sig2;
                params >> sig2;
                double shift6;
                params >> shift6; // must be 0.0 for full LJ
                if (calculateLJ) {
                    PotForceLJ(drs, dr2, eps24, sig2, f, u);
                    if(isCGSiteI) { u *= 1 - wi * wj;}
                    else { u *= wi * wj;}
                    u += shift6;
                    Upot6LJ += u;
                }
            }
        }
    }


    double m1[3], m2[3]; // angular momenta

    const unsigned nCG_C_i = components[mi.componentid() + 1].numCharges();
    const unsigned nCG_Q_i = components[mi.componentid() + 1].numQuadrupoles();
    const unsigned nCG_D_i = components[mi.componentid() + 1].numDipoles();
    const unsigned ne1 = mi.numCharges();
    const unsigned ne2 = mj.numCharges();
    const unsigned int nq1 = mi.numQuadrupoles();
    const unsigned int nq2 = mj.numQuadrupoles();
    const unsigned int nd1 = mi.numDipoles();
    const unsigned int nd2 = mj.numDipoles();
    for (unsigned si = 0; si < ne1; si++) {
        bool isCG_i = si < nCG_C_i;
        const std::array<double,3> dii = mi.charge_d_abs(si);
        // Charge-Charge
        for (unsigned sj = 0; sj < ne2; sj++) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.charge_d_abs(sj);
            double q1q2per4pie0; // 4pie0 = 1 in reduced units
            params >> q1q2per4pie0;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            PotForce2Charge(drs, dr2, q1q2per4pie0, f, u);
            if(isCG_i) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
        // Charge-Quadrupole
        for (unsigned sj = 0; sj < nq2; sj++) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
            double qQ05per4pie0; // 4pie0 = 1 in reduced units
            params >> qQ05per4pie0;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.quadrupole_e(sj);
            PotForceChargeQuadrupole(drs, dr2, ejj.data(), qQ05per4pie0, f, m2, u);
            if(isCG_i) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
        // Charge-Dipole
        for (unsigned sj = 0; sj < nd2; sj++) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.dipole_d_abs(sj);
            double minusqmyper4pie0;
            params >> minusqmyper4pie0;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.dipole_e(sj);
            PotForceChargeDipole(drs, dr2, ejj.data(), minusqmyper4pie0, f, m2, u);
            if(isCG_i) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
    }

    for (unsigned int si = 0; si < nq1; ++si) {
        bool isCG_i = si < nCG_Q_i;
        const std::array<double,3> dii = mi.quadrupole_d_abs(si);
        const std::array<double,3> eii = mi.quadrupole_e(si);

        // Quadrupole-Charge
        for (unsigned sj = 0; sj < ne2; sj++) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.charge_d_abs(sj);
            double qQ05per4pie0; // 4pie0 = 1 in reduced units
            params >> qQ05per4pie0;
            minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            PotForceChargeQuadrupole(drs, dr2, eii.data(), qQ05per4pie0, f, m1, u);
            if(isCG_i) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
        // Quadrupole-Quadrupole -------------------
        for (unsigned int sj = 0; sj < nq2; ++sj) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp;
                continue;
            }
            //double drs[3];
            const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
            double q2075;
            params >> q2075;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.quadrupole_e(sj);
            PotForce2Quadrupole(drs, dr2, eii.data(), ejj.data(), q2075, f, m1, m2, u);
            if(isCG_i) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
        // Quadrupole-Dipole -----------------------
        for (unsigned int sj = 0; sj < nd2; ++sj) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp;
                continue;
            }
            //double drs[3];
            const std::array<double,3> djj = mj.dipole_d_abs(sj);
            double qmy15;
            params >> qmy15;
            minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.dipole_e(sj);
            PotForceDiQuadrupole(drs, dr2, ejj.data(), eii.data(), qmy15, f, m2, m1, u);
            if(isCG_i) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
    }

    for (unsigned int si = 0; si < nd1; ++si) {
        bool isCG_i = si < nCG_D_i;
        const std::array<double,3> dii = mi.dipole_d_abs(si);
        const std::array<double,3> eii = mi.dipole_e(si);
        // Dipole-Charge
        for (unsigned sj = 0; sj < ne2; sj++) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.charge_d_abs(sj);
            double minusqmyper4pie0;
            params >> minusqmyper4pie0;
            minusSiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            PotForceChargeDipole(drs, dr2, eii.data(), minusqmyper4pie0, f, m1, u);
            if(isCG_i) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
        // Dipole-Quadrupole -----------------------
        for (unsigned int sj = 0; sj < nq2; ++sj) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp;
                continue;
            }
            //double drs[3];
            const std::array<double,3> djj = mj.quadrupole_d_abs(sj);
            double myq15;
            params >> myq15;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.quadrupole_e(sj);
            PotForceDiQuadrupole(drs, dr2, eii.data(), ejj.data(), myq15, f, m1, m2, u);
            if(isCG_i) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
        // Dipole-Dipole ---------------------------
        for (unsigned int sj = 0; sj < nd2; ++sj) {
            if((resolutionJ == Resolution::CoarseGrain && !isCG_i) ||
               (resolutionJ == Resolution::FullParticle && isCG_i)) {
                double tmp; params >> tmp; params >> tmp;
                continue;
            }
            const std::array<double,3> djj = mj.dipole_d_abs(sj);
            double my2;
            params >> my2;
            double rffac;
            params >> rffac;
            SiteSiteDistanceAbs(dii.data(), djj.data(), drs, dr2);
            const std::array<double,3> ejj = mj.dipole_e(sj);
            PotForce2Dipole(drs, dr2, eii.data(), ejj.data(), my2, rffac, f, m1, m2, u, MyRF);
            if(isCG_i) { u *= 1 - wi * wj;}
            else { u *= wi * wj;}
            UpotXpoles += u;
        }
    }
    // check whether all parameters were used
    mardyn_assert(params.eos());
}
