#include "ProfilerPMF.h"

InternalProfiler::InternalProfiler():cell_processor{nullptr},measured_steps{0}{

}

void InternalProfiler::init(ParticleContainer* pc, int bins, int freq){
    this->number_bins = bins;
    this->sample_frequency = freq;
    cell_processor=new InternalCellProcessor(global_simulation->getcutoffRadius(), this); 

    this->SetBinContainer(pc);
    this->InitRNodes();
    Log::global_log->info()<<"[PMF] Profiler enabled "<<std::endl;
}

void InternalProfiler::InitRNodes(){
    for(int i=0;i<number_bins;++i){
        double rmid;
        rmid = (i+0.5)*bin_width;
        r_nodes[i] = rmid;
    }  
}

void InternalProfiler::ResetBuffers(){

    std::fill(rdf_buffer.begin(),rdf_buffer.end(),0.0);

    std::fill(u_buffer.begin(),u_buffer.end(),0.0);

    std::fill(pairs_buffer.begin(),pairs_buffer.end(),0.0);
}


void InternalProfiler::SetBinContainer(ParticleContainer* pc){
    this->bin_width = pc->getCutoff()/number_bins;

    this->r_nodes.resize(number_bins);
    std::fill(r_nodes.begin(),r_nodes.end(),0.0);
    this->InitRNodes();

    this->rdf_buffer.resize(number_bins);
    std::fill(rdf_buffer.begin(),rdf_buffer.end(),0.0);

    this->u_buffer.resize(number_bins);
    std::fill(u_buffer.begin(),u_buffer.end(),0.0);

    this->pairs_buffer.resize(number_bins);
    std::fill(pairs_buffer.begin(),pairs_buffer.end(),0.0);

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
    //Add count for rdf
    this->rdf_buffer[index]++;
    //Add count for potential pairs and potential value itself
    this->pairs_buffer[index] ++;
    this->u_buffer[index] += pot;
    
}  

std::vector<double>& InternalProfiler::GetRDFValues(){
    return this->rdf_buffer;
}
std::vector<double>& InternalProfiler::GetPotentialValues(){
    return this->u_buffer;
}
std::vector<double>& InternalProfiler::GetRNodes(){
    return this->r_nodes;
}

void InternalProfiler::ProfileData(ParticleContainer* pc, unsigned long step){
    if(step%sample_frequency==0 && step > global_simulation->getInitStatistics()){
        measured_steps++;
        pc->traverseCells(*cell_processor);
    }
}

void InternalProfiler::GenerateInstantaneousData(ParticleContainer* particleContainer, Domain* domain){
    
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
        rdf_buffer[i] /= measured_steps;
        rdf_buffer[i] /= den;

        //Generate U(r)_{cg} values for interpolation
        double N = 0.5*pairs_buffer[i]*(pairs_buffer[i]-1.0);
        if(N>0){
            u_buffer[i]= u_buffer[i]/N;
        }else{
            u_buffer[i]=0;
        }
    }   
}