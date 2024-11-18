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
            com[i] = com[i]/total_mass;
        }
    }

    return com;
}

