#include "ProfilerPMF.h"

void InternalProfiler::init(ParticleContainer* pc, int bins, int freq){
    this->number_bins = bins;
    this->sample_frequency = freq;
    cell_processor=new InternalCellProcessor(global_simulation->getcutoffRadius(), this); 


}


void InternalProfiler::SetBinContainer(ParticleContainer* pc){
    this->bin_width = pc->getCutoff()/number_bins;

    this->bin_counts.resize(number_bins);
    std::fill(bin_counts.begin(),bin_counts.end(),0.0);

    this->pot_vals.resize(number_bins);
    std::fill(pot_vals.begin(),pot_vals.end(),0.0);

    this->pot_counts.resize(number_bins);
    std::fill(pot_counts.begin(),pot_counts.end(),0.0);

    Log::global_log->info()<<"[PMF] Internal profiler bin width "<<bin_width<<"\n";
    measured_distance_squared = bin_width*bin_width*number_bins*number_bins;
    Log::global_log->info()<<"[PMF] Limit distance "<<measured_distance_squared<<"\n";
}


std::array<double,3> InternalProfiler::GetCOM(Molecule* m){

    std::array<double,3> com{0.0,0.0,0.0};
    double total_mass;//TODO:Can be computed once for every component
    Component* comp = m->component();
    total_mass = comp->m();
    for(int lj=0;lj<comp->numLJcenters();lj++){
        LJcenter& lj_center = comp->ljcenter(lj);
        for(int idx=0;idx<lj_center.r().size();idx++){               
            com[idx] += lj_center.m()*m->ljcenter_d_abs(lj)[idx];
        }

    }
    for(int qs=0;qs<comp->numCharges();qs++){
        Charge& q_center = comp->charge(qs);
        for(int idx=0;idx<q_center.r().size();idx++){               
            com[idx] += q_center.m()*m->charge_d_abs(qs)[idx];
        }
    }

    for(int i=0;i<com.size();i++){
        com[i] /= total_mass;
    }

    return com;

}

void InternalProfiler::ProcessDistance(double r, double pot){  
    if(r > measured_distance_squared){ return;}

    int index = std::floor(std::sqrt(r)/bin_width);
    this->bin_counts[index]++;
    //here we compute the potential
    pot_counts[index] ++;
    pot_vals[index] += pot;

    
}  

void InternalProfiler::GenerateInstantaneousData(ParticleContainer* particleContainer, Domain* domain){
    //Generate RDF g(r) data
    for(int i=0;i<number_bins;++i){
        double rmin, rmax, rmid, binvol, rmin3,rmax3, den;
        rmid = (i+0.5)*bin_width;
        rmin = i*bin_width;
        rmax =(i+1)*bin_width;
        rmin3 = rmin*rmin*rmin;
        rmax3 = rmax*rmax*rmax;
        binvol = (4.0/3.0)*M_PI*(rmax3-rmin3);
        den = 0.5*domain->getglobalNumMolecules()*(domain->getglobalNumMolecules()-1.0)*binvol/domain->getGlobalVolume();
        bin_counts[i] /= measured_steps;
        bin_counts[i] /= den;
    }   
}