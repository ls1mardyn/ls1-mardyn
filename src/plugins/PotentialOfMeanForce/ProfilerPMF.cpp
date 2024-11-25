#include "ProfilerPMF.h"

InternalProfiler::InternalProfiler():cell_processor{nullptr},measured_steps{0}{


}

void InternalProfiler::ReadXML(XMLfileUnits& xmlconfig){
    //RDF config
    xmlconfig.getNodeValue("internalBins",number_bins);
    xmlconfig.getNodeValue("densityBins",density.total_bins);
    xmlconfig.getNodeValue("measureFreq",sample_frequency);
    xmlconfig.getNodeValue("outputFreq",output_frequency);

}

void InternalProfiler::init(double rc){

    this->SetBinContainer(rc);
    this->InitRNodes();

        density.InitBins(_simulation.getDomain()->getGlobalLength(0));
        cell_processor=new InternalCellProcessor(global_simulation->getcutoffRadius(), number_bins,bin_width,density.bin_width); 
    Log::global_log->info()<<"[PMF] Profiler enabled "<<std::endl;
}

void InternalProfiler::InitRNodes(){
    for(int i=0;i<number_bins;++i){
        double rmid;
        rmid = (i+0.5)*bin_width;
        r_nodes[i] = rmid;
    }  
}



void InternalProfiler::SetBinContainer(double rc){
    this->bin_width = rc/number_bins;

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
    return this->cell_processor->GetPairCountBuffer();
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

std::vector<double> InternalProfiler::GetRDF(){

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
        den = 0.5*_simulation.getDomain()->getglobalNumMolecules()*(_simulation.getDomain()->getglobalNumMolecules()-1.0)*binvol/_simulation.getDomain()->getGlobalVolume();
        instantaneous_rdf[i] = GetBinCounts()[i]/den;
        instantaneous_rdf[i] /= GetMeasuredSteps();
    }   
    return instantaneous_rdf;
}

std::vector<double> InternalProfiler::GetDensity(){

    std::vector<double> buffer_density;
    buffer_density.resize(density.total_bins);
    double binvol = _simulation.getDomain()->getGlobalVolume()/density.total_bins;
    double den = den = binvol*GetMeasuredSteps();
    for(int i=0;i<density.total_bins;++i){
        buffer_density[i] = GetDensityCounts()[i]/den;
    }   
    return buffer_density;
}

void InternalProfiler::PrintOutput2Files(unsigned long ss){
    if(ss>0 && ss%output_frequency == 0){
        std::string file_name=rdf_file_name_+std::to_string(ss)+".txt";
        std::ofstream output_file(file_name);

        for(int i=0;i<number_bins;++i){
            output_file<<std::setw(8)<<std::left<<r_nodes[i]<<"\t"<<std::setw(8)<<std::left<<GetRDF()[i]<<std::endl;
        }

        output_file.close();
        if(measure_density){
            file_name = "density_"+std::to_string(ss)+".txt";
            std::ofstream density_file(file_name);
            for(int i=0;i<density.bin_width;++i){
                density_file<<std::setw(8)<<std::left<<density.bin_centers[i]<<"\t"<<std::setw(8)<<std::left<<GetDensity()[i]<<std::endl;
            }

            density_file.close();

        }

    }
}