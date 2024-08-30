#include "PotentialProfiler.h"


PotentialProfiler::PotentialProfiler():cell_processor{nullptr},measured_steps{0}{

}

void PotentialProfiler::readXML(XMLfileUnits& file){
    file.getNodeValue("totalBins",number_bins);
    Log::global_log->info()<<"[Potential Profiler] Total bins "<<number_bins<<"\n";
    file.getNodeValue("sampleFrequency",sample_frequency);
    Log::global_log->info()<<"[Potential Profiler] Sample frequency "<<sample_frequency<<"\n";
    file.getNodeValue("averageFrequency", average_frequency);
    Log::global_log->info()<<"[Potential Profiler] Average frequency "<<average_frequency<<"\n";

}

void PotentialProfiler::init(ParticleContainer* pc, DomainDecompBase* dd, Domain* dom){
    SetBinContainer(pc);
    cell_processor=new CustomCellProcessor(global_simulation->getcutoffRadius(), this);
}

void PotentialProfiler::endStep(ParticleContainer* pc, DomainDecompBase* dd, Domain* dom, unsigned long simstep){

    //Generate && Print instantaneous data
    if(simstep%sample_frequency ==0 && 
    simstep > global_simulation->getInitStatistics()
    ){
        for(int i=0;i<pot_vals.size();++i){
            int N = 0.5*(double)pot_counts[i]*((double)pot_counts[i]-1.0);
            if(N>0 && !std::isinf(pot_vals[i])){
                // pot_vals[i] = pot_vals[i]/(0.5*(double)pot_counts[i]*((double)pot_counts[i]-1.0));
                pot_vals[i] = pot_vals[i]/pot_counts[i];
            }
            else{
                pot_vals[i]=0;
            }
        }
        pot_values_per_timestep.emplace_back(pot_vals);
        WriteDataToFile(pc, dom, simstep);
        
        std::fill(pot_vals.begin(),pot_vals.end(),0.0);
        std::fill(pot_counts.begin(),pot_counts.end(),0.0);

    }
        // for(int i=0;i<pot_vals.size();++i){
            // int N = 0.5*(double)pot_counts[i]*((double)pot_counts[i]-1.0);
            // if(N>0 && !std::isinf(pot_vals[i])){
                // pot_vals[i] = pot_vals[i]/(0.5*(double)pot_counts[i]*((double)pot_counts[i]-1.0));
            // }
            // else{
                // pot_vals[i]=0;
            // }
        // }
        // pot_values_per_timestep.emplace_back(pot_vals);
}

std::array<double,3> PotentialProfiler::GetCOM(Molecule* m){

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

void PotentialProfiler::SetBinContainer(ParticleContainer* pc){
    bin_width = pc->getCutoff()/(double)number_bins;
    
    this->bin_counts.resize(number_bins);
    std::fill(bin_counts.begin(),bin_counts.end(),0.0);

    this->pot_vals.resize(number_bins);
    std::fill(pot_vals.begin(),pot_vals.end(),0.0);

    this->pot_counts.resize(number_bins);
    std::fill(pot_counts.begin(),pot_counts.end(),0.0);

    this->current_potential_average.resize(number_bins);
    std::fill(current_potential_average.begin(),current_potential_average.end(),0.0);
    
    Log::global_log->info()<<"[PotentialProfiler] Bin  width "<<bin_width<<"\n";
    measured_distance_squared = bin_width*bin_width*number_bins*number_bins;
    Log::global_log->info()<<"[PotentialProfiler] Limit distance "<<measured_distance_squared<<"\n";
}

void PotentialProfiler::ProcessDistance(double r, double pot){  
    if(r > measured_distance_squared){ return;}

    int index = std::floor(std::sqrt(r)/bin_width);
    this->bin_counts[index]++;
    //here we compute the potential
    pot_counts[index] ++;
    pot_vals[index] += pot;
}  

void PotentialProfiler::WriteDataToFile(ParticleContainer* particleContainer, Domain* domain, unsigned long simstep){
    std::string filename="potential_data_"+std::to_string(simstep)+".txt";
    std::ofstream outfile(filename);

    //Averages over the stored data per timestep
    if(simstep%average_frequency==0){
        for(int d=0;d<pot_values_per_timestep.size();++d){
            for(int l=0;l<current_potential_average.size();++l){
                current_potential_average[l] += pot_values_per_timestep[d][l];
            }

        }

        pot_values_per_timestep.clear();
        
    }
    // std::vector<double> averaged_potential;
    // averaged_potential.resize(pot_vals.size());
    // for(int d=0;d<pot_values_per_timestep.size();++d){
        // for(int l=0;l<averaged_potential.size();++l){
            // averaged_potential[l] += pot_values_per_timestep[d][l];
        // }
    // }

    double data=0.0;
    for(int i=0;i<number_bins;i++){
        double rmid;
        rmid = (i+0.5)*bin_width;



        outfile<<std::setw(8)<<std::left<<rmid
        <<"\t"<<std::setw(8)<<std::left<<pot_vals[i];
        if(simstep%average_frequency==0){
            outfile<<"\t"<<std::setw(8)<<std::left<<current_potential_average[i]/measured_steps;
        }
        outfile<<std::endl;
    }

    outfile.close();
}



CustomCellProcessor::CustomCellProcessor(const double cr, PotentialProfiler* r):CellProcessor{cr,cr},my_profiler{r}{

}

double CustomCellProcessor::DistanceBetweenCOMs(std::array<double,3>& c1, std::array<double,3>& c2){
    double r =0.0;
    std::array<double,3> diff={0.0,0.0,0.0};
    
    for(int i=0;i<diff.size();i++){
        diff[i]=c1[i]-c2[i];
    }

    r = diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2];
    return r;

}

double CustomCellProcessor::PotentialCallBack(double e, double s, double r2){
    double pot=0;

    double inv_r2 = 1.0/r2;
    double lj6 = s*s*inv_r2;
    lj6 = lj6*lj6*lj6;
    double lj12 = lj6*lj6;
    double lj12m6 = lj12 - lj6;
    pot = lj12m6*e*4.0;
    if(std::isinf(pot)) return 0.0;
    return pot;
}

void CustomCellProcessor::processCell(ParticleCell& cell){
    if(cell.isInnerCell() || cell.isBoundaryCell()){
        auto begin = cell.iterator();
        double distance=0.0;
        for(auto it1 = begin;it1.isValid();++it1){
            std::array<double,3> com1={0.0,0.0,0.0};
            Molecule& m1 = *it1;
            double eps, sig, pot;
            sig = m1.component()->getSigma(0);
            eps = m1.component()->getEps(0);
            com1 = my_profiler->GetCOM(&m1);
            auto it2 = it1;
            ++it2;
            for(;it2.isValid();++it2){
                Molecule& m2 = *it2;
                std::array<double,3> com2={0.0,0.0,0.0};
                com2 = my_profiler->GetCOM(&m2);
                mardyn_assert(&m1 != &m2);
                //Now we compute the distance between the COMs
                distance = DistanceBetweenCOMs(com1,com2);
                if(distance < _cutoffRadiusSquare){

                    pot = PotentialCallBack(eps,sig,distance);
                    my_profiler->ProcessDistance(distance, pot);
                }
    
            }
    
        }
    }
}

void CustomCellProcessor::processCellPair(ParticleCell& c1, ParticleCell& c2, bool sumAll){

    auto begin1 = c1.iterator();
    auto begin2 = c2.iterator();
    double distance=0.0;
    std::array<double,3> com1={0.0,0.0,0.0};
    std::array<double,3> com2={0.0,0.0,0.0};
    if(sumAll){
        for(auto it1=begin1;it1.isValid();++it1){
            Molecule& m1 = *it1;
            double eps, sig, pot;
            sig = m1.component()->getSigma(0);
            eps = m1.component()->getEps(0);

            com1 = my_profiler->GetCOM(&m1);
            for(auto it2 =begin2;it2.isValid();++it2){
                Molecule& m2 = *it2;
                com2 = my_profiler->GetCOM(&m2);
                distance = DistanceBetweenCOMs(com1,com2);
                if(distance < _cutoffRadiusSquare){
                    pot = PotentialCallBack(eps,sig,distance);
                    my_profiler->ProcessDistance(distance,pot);
                }
            }
        }
    }
    else{
        if(c1.isInnerCell()){//no hallo cells at all

            for(auto it1=begin1;it1.isValid();++it1){
                Molecule& m1 = *it1;
                double eps, sig, pot;
                sig = m1.component()->getSigma(0);
                eps = m1.component()->getEps(0);
                com1 = my_profiler->GetCOM(&m1);
                for(auto it2=begin2;it2.isValid();++it2){
                    Molecule& m2 = *it2;
                    com2 = my_profiler->GetCOM(&m2);
                    distance = DistanceBetweenCOMs(com1,com2);
                    if(distance < _cutoffRadiusSquare){
                        pot = PotentialCallBack(eps,sig,distance);
                        my_profiler->ProcessDistance(distance,pot);
                    }

                }
            }

        }

        if(c1.isBoundaryCell()){//c1 is  boundary
            if(c2.isHaloCell() && !(c1.getCellIndex()<c2.getCellIndex())){
                return;
            }

            for(auto it1=begin1;it1.isValid();++it1){
                Molecule& m1 = *it1;
                double eps, sig, pot;
                sig = m1.component()->getSigma(0);
                eps = m1.component()->getEps(0);
                com1 = my_profiler->GetCOM(&m1);
                for(auto it2=begin2;it2.isValid();++it2){
                    Molecule& m2 = *it2;
                    com2 = my_profiler->GetCOM(&m2);
                    distance = DistanceBetweenCOMs(com1,com2);
                    
                    if(distance < _cutoffRadiusSquare){
                        pot = PotentialCallBack(eps,sig,distance);
                        my_profiler->ProcessDistance(distance,pot);
                    }
                }
            }
        }
    }

}