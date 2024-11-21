/*
 * Created on Thu Nov 21 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#include "Statistics.h"

void StatisticsAdResS::init(FPRegion& region){
    fp.low=region._low;
    fp.high = region._high;

    cg2.low= {region._highHybrid[0],0,0};
    cg2.high= {_simulation.getDomain()->getGlobalLength(0),_simulation.getDomain()->getGlobalLength(1),_simulation.getDomain()->getGlobalLength(2)};

    cg1.low={0,0,0};
    cg1.high={region._lowHybrid[0],region._lowHybrid[0],region._lowHybrid[0]};

    hy1.low=region._lowHybrid;
    hy1.high={region._low[0],region._low[0],region._low[0]};

    hy2.low={region._high[0],0,0};
    hy2.high=region._highHybrid;
}

void StatisticsAdResS::Output(){
    Log::global_log->info()<<"[AdResSStatistics] T_fp="<<temperatures[1]<<"\tT_cg="<<temperatures[0]<<"\tT_hy="<<temperatures[2]<<std::endl;
}

void StatisticsAdResS::MeasureCGTemperature(ParticleContainer* particleContainer){
    RegionParticleIterator it = particleContainer->regionIterator(cg1.low.data(),cg1.high.data(),ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    RegionParticleIterator it2 = particleContainer->regionIterator(cg2.low.data(),cg2.high.data(),ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    double E_kin=0;
    double T_cg1 =0;
    double mass =0;
    //Get total DOFs in CG1 region
    long long N = 0;
    long long N1 = 0;
    long long N2 = 0;

    mass = _simulation.getEnsemble()->getComponent(fp_component)->m();
    for(it;it.isValid();++it){
        E_kin += it->v2();
        ++N1;
    }

    for(it2;it2.isValid();++it2){
        E_kin += it2->v2();
        ++N2;
    }

    N = N1+N2;

    E_kin /= mass;
    T_cg1 = E_kin/(3.0*N);
    temperatures[0]=T_cg1;
}

void StatisticsAdResS::MeasureFPTemperature(ParticleContainer* particleContainer){

    RegionParticleIterator it = particleContainer->regionIterator(fp.low.data(),fp.high.data(),ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    double E_kin=0;
    double T_fp =0;
    double mass =0;
    //Get total DOFs in FP region
    long long N = _simulation.getEnsemble()->getComponent(fp_component)->getNumMolecules();
    N =0;
    mass = _simulation.getEnsemble()->getComponent(fp_component)->m();
    for(it;it.isValid();++it){
        E_kin += it->v2();
        ++N;
    }
    E_kin /= mass;
    T_fp = E_kin/(3.0*N);
    temperatures[1]=T_fp;

}

void StatisticsAdResS::MeasureHYTemperature(ParticleContainer* particleContainer){

    RegionParticleIterator it = particleContainer->regionIterator(hy1.low.data(),hy1.high.data(),ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    RegionParticleIterator it2 = particleContainer->regionIterator(hy2.low.data(),hy2.high.data(),ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    double E_kin=0;
    double T_hy =0;
    double mass =0;
    //Get total DOFs in HY region
    long long N = 0;
    long long N1 = 0;
    long long N2 = 0;

    mass = _simulation.getEnsemble()->getComponent(fp_component)->m();
    for(it;it.isValid();++it){
        E_kin += it->v2();
        ++N1;
    }

    for(it2;it2.isValid();++it2){
        E_kin += it2->v2();
        ++N2;
    }

    N = N1+N2;

    E_kin /= mass;
    T_hy = E_kin/(3.0*N);
    temperatures[2]=T_hy;

}