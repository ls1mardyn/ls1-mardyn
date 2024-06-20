#include "PMF.h"

void InteractionSite::AddForce(std::array<double,3> f){
    for(int i=0;i<f.size();i++){
        this->f_com[i] += f[i];
    }
}

void InteractionSite::AddPotential(double u){
    this->u_com += u;
}