/*
 * Created on Thu Apr 25 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */


#pragma once

#include<string>
#include<numeric>
#include<fstream>
//#include<iostream>
template<class T>
class DataLog{

    public:

    void ColumnOutput();

    private:

    std::string file_name;
    std::ofstream file;

};