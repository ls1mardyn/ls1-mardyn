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
    fp.component_name="FP";

    cg2.low= {region._highHybrid[0],0,0};
    cg2.high= {_simulation.getDomain()->getGlobalLength(0),_simulation.getDomain()->getGlobalLength(1),_simulation.getDomain()->getGlobalLength(2)};
    cg2.component_name="CG";

    cg1.low={0,0,0};
    cg1.high={region._lowHybrid[0],_simulation.getDomain()->getGlobalLength(1),_simulation.getDomain()->getGlobalLength(2)};
    cg1.component_name="CG";

    hy1.low=region._lowHybrid;
    hy1.high={region._low[0],_simulation.getDomain()->getGlobalLength(1),_simulation.getDomain()->getGlobalLength(2)};
    hy1.component_name="HY";
    hy2.low={region._high[0],0,0};
    hy2.high=region._highHybrid;
    hy2.component_name="HY";
    //open for writing and clear
    statistics.open(file_name,std::ofstream::out|std::ofstream::trunc);
    statistics.close();
    statistics.open(file_name,std::ofstream::app);
    statistics<<"Step \t T_{fp}\t N_{fp} \t T_{hy1} \t N_{hy1} \t T_{hy2}\t N_{hy2}\t T_{hy}\t N_{hy}\t T_{cg1} \t N_{cg1} \t T_{cg2}\t N_{cg2}\t T_{cg}\t N_{cg}"<<std::endl;
    statistics.close();

    //setup profilers
    fp_profiler.init(_simulation.getcutoffRadius());
    fp_profiler.SetFilePrefix("fp_rdf_");
    cg1_profiler.init(_simulation.getcutoffRadius());
    cg1_profiler.SetFilePrefix("cg_rdf_");
    fp_theoreticalRdf = TheoreticalRDF(fp.low, fp.high, 100);
    cg_theoreticalRdf = TheoreticalRDF(cg1.low, cg1.high, 100);
}

void StatisticsAdResS::Output2File(long step){
    statistics.open(file_name,std::ofstream::app);
    if(step % output_stride ==0){
        size_t width = 4;
        statistics<<step<<"\t "
                  <<std::setw(width)<<std::left<<fp_temp<<"\t"
                  <<std::setw(width)<<std::left<<N_fp<<"\t"
                  <<std::setw(width)<<std::left<<hy1_temp<<"\t"
                  <<std::setw(width)<<std::left<<N_hy1<<"\t"
                  <<std::setw(width)<<std::left<<hy2_temp<<"\t"
                  <<std::setw(width)<<std::left<<N_hy2<<"\t"
                  <<std::setw(width)<<std::left<<hy_temp<<"\t"
                  <<std::setw(width)<<std::left<<(N_hy2+N_hy1)<<"\t"
                  <<std::setw(width)<<std::left<<cg1_temp<<"\t"
                  <<std::setw(width)<<std::left<<N_cg1<<"\t"
                  <<std::setw(width)<<std::left<<cg2_temp<<"\t"
                  <<std::setw(width)<<std::left<<N_cg2<<"\t"
                  <<std::setw(width)<<std::left<<cg_temp<<"\t"
                  <<std::setw(width)<<std::left<<(N_cg2+N_cg1)<<"\t"
                  <<std::endl;


    }
    statistics.close();
}

void StatisticsAdResS::MeasureCGTemperature(ParticleContainer* particleContainer){
    RegionParticleIterator it = particleContainer->regionIterator(cg1.low.data(),cg1.high.data(),ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    RegionParticleIterator it2 = particleContainer->regionIterator(cg2.low.data(),cg2.high.data(),ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    double E_kin1=0;
    double E_kin2=0;
    double mass =0;
    //Get total DOFs in CG1 region
    int N1 = 0;
    int N2 = 0;

    mass = _simulation.getEnsemble()->getComponents()->at(2).m();
    for(it;it.isValid();++it){
        E_kin1 += it->v2();
        ++N1;
    }
    
    E_kin1 *= mass;
    cg1_temp = E_kin1/(3.0*N1);
    N_cg1 = N1;

    for(it2;it2.isValid();++it2){
        E_kin2 += it2->v2();
        ++N2;
    }

    E_kin2 *= mass;
    cg2_temp = E_kin2/(3.0*N2);
    N_cg2 = N2;

    cg_temp = (E_kin1+E_kin2)/(3.0*(N1+N2));
}

void StatisticsAdResS::MeasureFPTemperature(ParticleContainer* particleContainer){

    RegionParticleIterator it = particleContainer->regionIterator(fp.low.data(),fp.high.data(),ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    double E_kin=0;
    double T_fp =0;
    double mass =0;
    //Get total DOFs in FP region
    int N =0;
    mass = _simulation.getEnsemble()->getComponents()->at(0).m();
    for(it;it.isValid();++it){
        E_kin += it->v2();
        ++N;
    }
    E_kin *= mass;
    T_fp = E_kin/(3.0*N);
    N_fp = N;
    fp_temp = T_fp;
    
}

void StatisticsAdResS::MeasureHYTemperature(ParticleContainer* particleContainer){

    RegionParticleIterator it = particleContainer->regionIterator(hy1.low.data(),hy1.high.data(),ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    RegionParticleIterator it2 = particleContainer->regionIterator(hy2.low.data(),hy2.high.data(),ParticleIterator::ONLY_INNER_AND_BOUNDARY);
    double E_kin1=0;
    double E_kin2=0;
    double T_hy =0;
    double mass =0;
    //Get total DOFs in HY region
    long long N1 = 0;
    long long N2 = 0;

    mass = _simulation.getEnsemble()->getComponents()->at(1).m();
    for(it;it.isValid();++it){
        E_kin1 += it->v2();
        ++N1;
    }
    E_kin1 *= mass;
    hy1_temp = E_kin1/(3.0*N1);
    N_hy1 = N1;

    for(it2;it2.isValid();++it2){
        E_kin2 += it2->v2();
        ++N2;
    }

    E_kin2 *= mass;
    hy2_temp = E_kin2/(3.0*N2);
    N_hy2 = N2;

    hy_temp = (E_kin1+E_kin2)/(3.0*(N1+N2));
}

void StatisticsAdResS::ClearAll(){
    fp_temp=0;
    cg1_temp=0;
    cg2_temp=0;
    hy1_temp=0;
    hy2_temp=0;
    N_cg1=0;
    N_cg2=0;
    N_hy1=0;
    N_hy2=0;
    N_fp=0;
}