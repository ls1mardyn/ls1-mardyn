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

inline double Distance2BetweenCOMs(std::array<double,3>& c1, std::array<double,3>& c2){
    double r =0.0;
    std::array<double,3> diff={0.0,0.0,0.0};
    
    for(int i=0;i<diff.size();i++){
        diff[i]=c1[i]-c2[i];
    }

    r = diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2];
    return r;
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

inline void VectorSub(std::vector<double>& v1, const std::vector<double>& v2){
    if(v1.size()!=v2.size()){
        Log::global_log->error()<<"[VectorAdd]Not same size add"<<std::endl;
    }
    for(int i=0;i<v1.size();++i){
        if(!std::isfinite(v2[i])) continue;
        v1[i] -= v2[i];
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

inline void FivePointAverageExtrapolation(const std::vector<double>& Vx, std::vector<double>& Vy){
    for(int i=Vy.size()-1;i>=0;--i){
        if(std::isfinite(Vy[i])){
            continue;
        }

        double y_k = Vy[i+1];
        double y_k_1 = Vy[i+2];
        double y_k_2 = Vy[i+3];
        double y_k_3 = Vy[i+4];
        double y_k_4 = Vy[i+5];

        double x_k = Vx[i+1];
        double x_k_1 = Vx[i+2];
        double h = x_k-x_k_1;




        Vy[i] = (y_k+2.0*y_k_1+3.0*y_k_2+4.0*y_k_3+5.0*y_k_4)/5.0;


    }
}

inline void FivePointDifferenceExtrapolation(const std::vector<double>& Vx, std::vector<double>& Vy){
    for(int i=Vy.size()-1;i>=0;--i){
        if(std::isfinite(Vy[i])){
            continue;
        }

        double y_k = Vy[i+1];
        double y_k_1 = Vy[i+2];
        double y_k_2 = Vy[i+3];
        double y_k_3 = Vy[i+4];
        double y_k_4 = Vy[i+5];

        double x_k = Vx[i+1];
        double x_k_1 = Vx[i+2];
        double h = x_k-x_k_1;


        double dF1 = 25.0*y_k - 48.0*y_k_1+ 36.0*y_k_2- 16.0*y_k_3+ 3.0*y_k_4/(12.0*h);
        double dF2 = 35.0*y_k - 104.0*y_k_1+ 114.0*y_k_2- 56.0*y_k_3+ 11.0*y_k_4/(12.0*h*h);
        double dF3 = 15.0*y_k - 49.0*y_k_1+ 78.0*y_k_2- 52.0*y_k_3+ 10.0*y_k_4/(6.0*h*h*h);
        double dF4 = 5.0*y_k - 20.0*y_k_1+ 30.0*y_k_2- 20.0*y_k_3+ 5.0*y_k_4/(h*h*h*h);



        Vy[i] = y_k_1 + h*dF1 + 0.5*h*h*dF2 + h*h*h/6.0 *dF3+ h*h*h*h/24.0 *dF4;


    }
}

inline void LinearExtrapolation(const std::vector<double>& Vx, std::vector<double>& Vy){
    

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