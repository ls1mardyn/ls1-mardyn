#include "PMF.h"

void InteractionSite::SubForce(std::array<double,3> f){
    for(int i=0;i<f.size();i++){
        this->f_com[i] -= f[i];
    }
}

void InteractionSite::AddForce(std::array<double,3> f){
    for(int i=0;i<f.size();i++){
        this->f_com[i] += f[i];
    }
}

void InteractionSite::AddPotential(double u){
    this->u_com += u;
}

void InteractionSite::SetPosition(std::array<double, 3> pos){
    this->_r=pos;
}

void InteractionSite::SetVelocity(std::array<double, 3> vel){
    this->v_com = vel;
}

std::array<double,3>& InteractionSite::GetPosition(){
    return this->_r;
}

std::array<double,3>& InteractionSite::GetVelocity(){
    return this->v_com;
}
