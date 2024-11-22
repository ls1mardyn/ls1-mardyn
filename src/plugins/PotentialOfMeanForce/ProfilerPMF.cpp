#include "ProfilerPMF.h"

InternalProfiler::InternalProfiler():cell_processor{nullptr},measured_steps{0}{


}

void InternalProfiler::init(ParticleContainer* pc, int bins, int freq){
    this->number_bins = bins;
    this->sample_frequency = freq;
    this->SetBinContainer(pc);
    this->InitRNodes();
    cell_processor=new InternalCellProcessor(global_simulation->getcutoffRadius(), bins,bin_width); 
    Log::global_log->info()<<"[PMF] Profiler enabled "<<std::endl;
}

void InternalProfiler::InitRNodes(){
    for(int i=0;i<number_bins;++i){
        double rmid;
        rmid = (i+0.5)*bin_width;
        r_nodes[i] = rmid;
    }  
}

void InternalProfiler::ResetBuffer(){

    std::fill(rdf_buffer.begin(),rdf_buffer.end(),0.0);
}


void InternalProfiler::SetBinContainer(ParticleContainer* pc){
    this->bin_width = pc->getCutoff()/number_bins;

    this->r_nodes.resize(number_bins);
    std::fill(r_nodes.begin(),r_nodes.end(),0.0);
    this->InitRNodes();

    this->rdf_buffer.resize(number_bins);
    std::fill(rdf_buffer.begin(),rdf_buffer.end(),0.0);;

    Log::global_log->info()<<"[PMF] Internal profiler bin width "<<bin_width<<"\n";
    measured_distance_squared = bin_width*bin_width*number_bins*number_bins;
    Log::global_log->info()<<"[PMF] Limit distance "<<measured_distance_squared<<"\n";

}

double InternalProfiler::GetMeasuredSteps(){
    return measured_steps;
}

std::vector<double>& InternalProfiler::GetBinCounts(){
    return this->rdf_buffer;
}

std::vector<double>& InternalProfiler::GetBinCenters(){
    return this->r_nodes;
}

void InternalProfiler::ProfileData(ParticleContainer* pc, unsigned long step){
    if(step%sample_frequency==0){
        measured_steps++;
        pc->traverseCells(*cell_processor);
    }
}

std::vector<double> InternalProfiler::GetInstantaneousData(Domain* domain){

    std::vector<double> instantaneous_rdf;
    instantaneous_rdf.resize(number_bins);
    
    for(int i=0;i<number_bins;++i){
        //Generate RDF g(r) data && U(r) instantaneous values
        double rmin, rmax, rmid, binvol, rmin3,rmax3, den;
        rmid = (i+0.5)*bin_width;
        rmin = i*bin_width;
        rmax =(i+1)*bin_width;
        rmin3 = rmin*rmin*rmin;
        rmax3 = rmax*rmax*rmax;
        binvol = (4.0/3.0)*M_PI*(rmax3-rmin3);
        den = 0.5*domain->getglobalNumMolecules()*(domain->getglobalNumMolecules()-1.0)*binvol/domain->getGlobalVolume();
        // rdf_buffer[i] /= den;
        instantaneous_rdf[i] = rdf_buffer[i]/den;

    }   
    ResetBuffer();
    return instantaneous_rdf;
}