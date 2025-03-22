#include "Convergence.h"

Convergence::Convergence(double t, const std::string& n):tolerance{t}, name(n){

}

bool Convergence::CheckConvergence(std::vector<double>& rdf_ref, std::vector<double>& rdf_i){

    double conv=0.0;
    for(int i=0;i<rdf_ref.size();++i){
        conv += std::abs(rdf_ref[i]-rdf_i[i]);
    }
    /*double sum=0;
    double sub=0;
    for(int i=0;i<rdf_ref.size();++i){
        sum += (std::abs(rdf_i[i])+std::abs(rdf_ref[i]));
        sub += std::abs(rdf_i[i]-rdf_ref[i]);
    }

    conv = 1 - (sub/sum);*/
    
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


void Convergence::PrintLocalConvergence2File(){
    std::string conv_name = "local_convergence_"+name+"_"+std::to_string(ibi_iteration)+".txt";
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
    std::string conv_name = "global_convergence_"+name+".txt";
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