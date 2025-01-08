#pragma once

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"

inline std::array<double, 3> ComputeCOM(Molecule& m1){
    std::array<double, 3> com{0,0,0};
    double total_mass = m1.component()->m();

    for(int lj=0;lj<m1.component()->numLJcenters();++lj){
        auto lj_site = m1.ljcenter_d_abs(lj);
        double site_mass = m1.component()->ljcenter(lj).m();
        for(int i=0;i<3;++i){
            com[i] += lj_site[i]*site_mass;
            // com[i] = com[i]/total_mass;
        }
    }


    for(int i=0;i<3;++i){
        com[i] = com[i]/total_mass;
    }

    return com;
}

inline void VectorAdd(std::vector<double>& v1, const std::vector<double>& v2){
    if(v1.size()!=v2.size()){
        Log::global_log->error()<<"[VectorAdd]Not same size add"<<std::endl;
    }
    for(int i=0;i<v1.size();++i){
        if(!std::isfinite(v2[i])) continue;
        v1[i] += v2[i];
    }
}

inline void NormalizeVector(std::array<double,3>& V){
    double norm = 0.0;
    for(int i=0;i<V.size();++i){
        norm += std::pow(V[i],2.0);
    }
    norm = std::sqrt(norm);
    //TODO:assert norm is not zero
    for(int i=0;i<V.size();++i){
        V[i] /= norm;
    }

}

inline void ExtrapolateVector(const std::vector<double>& Vx, std::vector<double>& Vy){


    for(int i= Vy.size()-1;i>=0;--i){

        if(std::isfinite(Vy[i])){
            continue;
        }

        double y_k = Vy[i+1];
        double y_k_1 = Vy[i+2];
        double x_k = Vx[i+1];
        double x_k_1 = Vx[i+2];
        double x_star = Vx[i];
        Vy[i] = y_k_1 + 1.0*(x_star-x_k_1)/(x_k-x_k_1)*(y_k-y_k_1);
        
    }








}