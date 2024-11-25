#include "Convergence.h"

Convergence::Convergence(double t):tolerance{t}{

}

bool Convergence::CheckConvergence(std::vector<double>& rdf_ref, std::vector<double>& rdf_i){

    double conv=0.0;
    for(int i=0;i<rdf_ref.size();++i){
        conv += std::abs(rdf_ref[i]-rdf_i[i]);
    }
    local_convergence.emplace_back(conv);

    if(conv<tolerance){
        return true;
    }

    return false;

}

void Convergence::PrepareUpdate(){
    int local_steps = local_convergence.size();
    ibi_convergence.emplace_back(local_convergence[local_steps-1]);
    ++ibi_iteration;            
    local_convergence.clear();
}

bool Convergence::TriggerPotentialUpdate(){

    bool trigger = false;
    int local_steps = local_convergence.size();
    
    if(local_steps>100){
        if((local_convergence[local_steps-1]-local_convergence[local_steps-2])>0){
            trigger = true;
            ibi_convergence.emplace_back(local_convergence[local_steps-1]);
            ++ibi_iteration;            

            local_convergence.clear();
        }
    }

    

    return trigger;

}

void Convergence::PrintLocalConvergence2File(int step){
    std::string conv_name = "local_convergence_ibi_"+std::to_string(ibi_iteration)+".txt";
    std::ofstream conv{conv_name};
    for(int i=0;i<local_convergence.size();++i){
        conv<<std::setw(8)<<std::left
        <<i
        <<"\t"
        <<std::setw(8)<<std::left<<local_convergence[i]
        <<std::endl;
    }
    conv.close();
}

void Convergence::PrintGlobalConvergence2File(){
    std::string conv_name = "global_convergence_ibi.txt";
    std::ofstream conv{conv_name};
    for(int i=0;i<ibi_convergence.size();++i){
        conv<<std::setw(8)<<std::left
        <<i
        <<"\t"
        <<std::setw(8)<<std::left<<ibi_convergence[i]
        <<std::endl;
    }
    conv.close();
}