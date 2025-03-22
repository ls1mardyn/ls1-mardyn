/*
 * Created on Fri Oct 04 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once

#include<vector>
#include<math.h>
#include<iostream>
#include<fstream>
#include<iomanip>


class Convergence{

    private:
    int ibi_iteration=0;
    double tolerance;
    std::vector<double> ibi_convergence;
    std::vector<double> local_convergence;
    std::string name;


    public:

    Convergence(double tol, const std::string& name);
    bool CheckConvergence(std::vector<double>& rdf_ref, std::vector<double>& rdf_i);
    void PrintLocalConvergence2File();
    void PrintGlobalConvergence2File();
    void PrepareUpdate();
};