#include "ProfilerPMF.h"

InternalProfiler::InternalProfiler():cell_processor{nullptr},measured_steps{0}{


}

void InternalProfiler::init(ParticleContainer* pc, int bins, int freq,int density_bins){
    this->number_bins = bins;
    this->sample_frequency = freq;
    this->SetBinContainer(pc);
    this->InitRNodes();
    density.InitBins(density_bins,_simulation.getDomain()->getGlobalLength(0));
    cell_processor=new InternalCellProcessor(global_simulation->getcutoffRadius(), bins,bin_width,density.bin_width); 


    Log::global_log->info()<<"[PMF] Profiler enabled "<<std::endl;
}

void InternalProfiler::InitRNodes(){
    for(int i=0;i<number_bins;++i){
        double rmid;
        rmid = (i+0.5)*bin_width;
        r_nodes[i] = rmid;
    }  
}



void InternalProfiler::SetBinContainer(ParticleContainer* pc){
    this->bin_width = pc->getCutoff()/number_bins;

    this->r_nodes.resize(number_bins);
    std::fill(r_nodes.begin(),r_nodes.end(),0.0);
    this->InitRNodes();

    Log::global_log->info()<<"[PMF] Internal profiler bin width "<<bin_width<<"\n";
    measured_distance_squared = bin_width*bin_width*number_bins*number_bins;
    Log::global_log->info()<<"[PMF] Limit distance "<<measured_distance_squared<<"\n";

}

double InternalProfiler::GetMeasuredSteps(){
    return measured_steps;
}

std::vector<double>& InternalProfiler::GetBinCounts(){
    return this->cell_processor->GetBuffer();
}

std::vector<double>& InternalProfiler::GetDensityCounts(){
    return this->cell_processor->GetDensityBuffer();
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
        instantaneous_rdf[i] = GetBinCounts()[i]/den;

    }   
    return instantaneous_rdf;
}

std::vector<double> InternalProfiler::GetDensity(Domain* domain){

    std::vector<double> buffer_density;
    buffer_density.resize(density.total_bins);
    double binvol = domain->getGlobalVolume()/density.total_bins;
    double den = den = binvol*GetMeasuredSteps();
    for(int i=0;i<density.total_bins;++i){

        // double binvol,den;
        // binvol = domain->getGlobalVolume()/number_bins;
        // den = binvol*GetMeasuredSteps();
        buffer_density[i] = GetDensityCounts()[i]/den;

    }   
    return buffer_density;
}

