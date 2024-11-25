#pragma once

#include<vector>
#include<math.h>

struct ConvergenceTest{
        bool ConvergenceCheck(std::vector<double>& ref, std::vector<double>& crrnt){
            double conv=0.0;
            for(int i=0;i<ref.size();++i){
                conv += std::abs(ref[i]-crrnt[i]);
            }
            convergence_per_step.emplace_back(conv);
            ++ibi_iteration;
            if(conv<tolerance)
            return true;

            return false;
        }
        double tolerance=0.05;
        std::vector<double> convergence_per_step;
        int ibi_iteration=0;
};