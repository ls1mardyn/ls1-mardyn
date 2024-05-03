/*
 * Created on Wed Apr 24 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */


#pragma once
#include<array>
#include<cmath>

template<class T> 
std::array<int, 3> MapGlobalToLocal(int idx, T vec){
        std::array<int, 3> local_indeces;
        local_indeces[0]= std::fmod(idx, vec[0]);
        local_indeces[1]= std::fmod(idx/vec[0],vec[1]);
        local_indeces[2]= idx/(vec[0]*vec[1]);

        return local_indeces;
}

template<class T>
int MapLocalToGlobal(std::array<int, 3>& local_indeces, T vec){
    int global_index = local_indeces[0] + local_indeces[1]*vec[0] + local_indeces[2]*vec[0]*vec[1];

    return global_index;
} 

template<class T>
int MapLocalToGlobal(int i, int j, int k, T vec){
    T idcs = {i,j,k};
    return MapLocalToGlobal(idcs, vec);
}





