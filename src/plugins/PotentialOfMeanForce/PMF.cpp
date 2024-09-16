#include"PMF.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"


PMF::PMF():reference_rdf_interpolation{1.0},potential_interpolation{0.0},avg_rdf_interpolation{1.0}{
}


void PMF::init(ParticleContainer* pc, DomainDecompBase* domainDecomp, Domain* domain){
    
    this->ReadRDF();

    this->InitializePotentialValues();

    Log::global_log->info()<<"[PMF] RDF has been read successfully\n";
    pairs_handler = new InteractionForceAdapter(resolution_handler,this);
    _simulation.setParticlePairsHandler(pairs_handler);
    _simulation.setCellProcessor(new LegacyCellProcessor(_simulation.getcutoffRadius(), _simulation.getLJCutoff(), pairs_handler));

    Log::global_log->info()<<"[PMF]LegacyCellProcessor set\n";
    Log::global_log->info()<<"[PMF] ForcedAdapter Class being used\n";

    Log::global_log->info()<<"[PMF] Start the tracker sites\n";
    for(auto it= pc->iterator(ParticleIterator::ALL_CELLS);it.isValid();++it){
        unsigned long m_id = it->getID();
        std::array<double,3> com = profiler.GetCOM(&(*it));
        sites[m_id].first.SetPosition(com);
    }

    for(auto it= begin(regions);it!=end(regions);++it){
        Log::global_log->info()<<"[PMF] The hybrid region spans from: ("<<it->_lowHybrid[0]<<","<<it->_lowHybrid[1]<<","<<it->_lowHybrid[2]<<") to ("<<it->_highHybrid[0]<<","<<it->_highHybrid[1]<<","<<it->_highHybrid[2]<<")"<<std::endl;

        Log::global_log->info()<<"[PMF] The atomistic region spans from: ("<<it->_low[0]<<","<<it->_low[1]<<","<<it->_low[2]<<") to ("<<it->_high[0]<<","<<it->_high[1]<<","<<it->_high[2]<<")"<<std::endl;

        Log::global_log->info()<<"[PMF] The region  center is located at: ("<<it->_center[0]<<","<<it->_center[1]<<","<<it->_center[2]<<")"<<std::endl;
    }
    Log::global_log->info()<<"[PMF] Initializing the COM sites\n";
    resolution_handler.CheckResolution(pc,sites,regions);
    Log::global_log->info()<<"[PMF] Enabled "<<std::endl;

    Log::global_log->info()<<"[PMF] InternalProfiler uses "<<internal_bins
                           <<" bins, measures every "<<measure_frequency<<" steps\n";
    this->profiler.init(pc,internal_bins,measure_frequency);
    Log::global_log->info()<<"[PMF] InternalProfiler initialized\n";
    Log::global_log->info()<<"[PMF] Step size of "<<multiplier<<"\n";
    potential_interpolation.SetXValues(profiler.GetRNodes());
    avg_rdf_interpolation.SetXValues(profiler.GetRNodes());
    avg_rdf_interpolation.GetYValues().resize(internal_bins);
}

void PMF::readXML(XMLfileUnits& xmlfile){

    //create AT regions
    int num_regions =0;
    XMLfile::Query query = xmlfile.query("fpregions/region");
    num_regions = query.card();

    this->regions.resize(num_regions);

    XMLfile::Query::const_iterator region_iterator;
    std::string oldpath = xmlfile.getcurrentnodepath();
    for(region_iterator = query.begin();region_iterator;region_iterator++){
        xmlfile.changecurrentnode(region_iterator);
        unsigned int id=0;
        xmlfile.getNodeValue("@id",id);
        regions[id-1].readXML(xmlfile);
    }
    xmlfile.changecurrentnode(oldpath);
    xmlfile.getNodeValue("multiplier",multiplier);
    xmlfile.getNodeValue("internalBins",internal_bins);
    xmlfile.getNodeValue("measureFreq",measure_frequency);
    xmlfile.getNodeValue("output",output);
    Log::global_log->info()<<"[PMF] Target temperature = "<<_simulation.getEnsemble()->T()<<std::endl;
}

void PMF::beforeForces(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, unsigned long simstep){

}

void PMF::beforeEventNewTimestep(ParticleContainer* pc, DomainDecompBase* domainDecomp, unsigned long simstep){

    for(auto it= pc->iterator(ParticleIterator::ALL_CELLS);it.isValid();++it){
        unsigned long m_id = it->getID();
        std::array<double,3> com = profiler.GetCOM(&(*it));
        sites[m_id].first.SetPosition(com);
    }

    resolution_handler.CheckResolution(pc,sites,regions);
}

void PMF::afterForces(ParticleContainer* pc, DomainDecompBase* dd, unsigned long step){
    profiler.ProfileData(pc,step);
}

void PMF::endStep(ParticleContainer* pc, DomainDecompBase* dd, Domain* domain, unsigned long step){

    AccumulateRDF(pc,domain);
    
    if(output){
        std::string filename="avg_rdf_"+std::to_string(step)+".txt";
        std::ofstream rdf_file(filename);
        for(int i=0;i<avg_rdf_interpolation.GetYValues().size();++i){
            rdf_file<<std::setw(8)<<std::left<<avg_rdf_interpolation.GetXValues()[i]<<"\t"<<std::setw(8)<<std::left<<GetAverageRDF()[i]<<std::endl;
        }
        rdf_file.close();
        std::string potential_file = "potential_"+std::to_string(step)+".txt";
        std::ofstream potential(potential_file);
        for(int i=0;i<potential_interpolation.GetYValues().size();++i){
            potential<<std::setw(8)<<std::left
            <<potential_interpolation.GetXValues()[i]
            <<"\t"
            <<std::setw(8)<<std::left<<potential_interpolation.GetYValues()[i]<<std::endl;
        }
        potential.close();
    }

    AddPotentialCorrection(step);
    std::vector<double> current_rdf = GetAverageRDF();
    convergence_check.ConvergenceCheck(reference_rdf_interpolation.GetYValues(),current_rdf);
    std::string conv_name = "convergence.txt";
    std::ofstream conv{conv_name};
    for(int i=0;i<convergence_check.convergence_per_step.size();++i){
        conv<<std::setw(8)<<std::left
        <<i
        <<"\t"
        <<std::setw(8)<<std::left<<convergence_check.convergence_per_step[i]
        <<std::endl;
    }
    conv.close();

}

void PMF::siteWiseForces(ParticleContainer* pc, DomainDecompBase* dd, unsigned long step){

}
/********************
 * ****************** FUNCTIONS NOT FROM THE INTERFACE
 *******************/


void PMF::InitializePotentialValues(){
    std::vector<double> pot0 = reference_rdf_interpolation.GetYValues();
    for(int i=0;i<pot0.size();++i){
        pot0[i] = -1.0*_simulation.getEnsemble()->T()*std::log(pot0[i]);
    }

    potential_interpolation.SetYValues(pot0);
}

void PMF::AddPotentialCorrection(unsigned long step){
    std::vector<double> pot_i;
    std::vector<double> avg_rdf = GetAverageRDF();
    std::vector<double> current_correction;
    current_correction.resize(internal_bins);
    pot_i.resize(internal_bins);

    for(int i=0;i<pot_i.size();++i){
        double ratio = avg_rdf[i]/reference_rdf_interpolation.GetYValues()[i];
        current_correction[i] = std::log(ratio);
        pot_i[i] = potential_interpolation.GetYValues()[i] - multiplier* _simulation.getEnsemble()->T()*current_correction[i];
    }

    potential_interpolation.SetYValues(pot_i);

    if(output){
        std::string name="correction_step_"+std::to_string(step)+".txt";
        std::ofstream corr{name};
        for(int i=0;i<internal_bins;++i){
            corr<<reference_rdf_interpolation.GetXValues()[i]
                <<current_correction[i]
                <<std::endl;
        }

    }

}

void PMF::AccumulateRDF(ParticleContainer* pc, Domain* dom){
    this->profiler.GenerateInstantaneousData(pc,dom);

    std::vector<double>& current_rdf = profiler.GetRDFValues();
    std::vector<double>& accumulated_rdf = avg_rdf_interpolation.GetYValues();

    if(current_rdf.size()!=accumulated_rdf.size()){
        Log::global_log->error()<<"[PMF] Buffer sizes do not match"<<std::endl;
    }

    for(int i=0;i<current_rdf.size();++i){

        accumulated_rdf[i] += current_rdf[i];

    }
    avg_rdf_interpolation.SetYValues(accumulated_rdf);
    profiler.ResetBuffers();
}

double PMF::WeightValue(const std::array<double,3>& pos, FPRegion& region){
    return weight_function.WeightValue(pos,region);
}

std::vector<FPRegion>& PMF::GetRegions(){
    return this->regions;
}

ResolutionType PMF::GetMoleculeResolution(unsigned long idx){
    return sites[idx].second;
}

InteractionSite PMF::GetMoleculeCOMSite(unsigned long idx){
    return sites[idx].first;
}

Interpolate& PMF::GetRDFInterpolation(){
    return this->reference_rdf_interpolation;
}

Interpolate& PMF::GetAVGRDFInterpolation(){
    return this->avg_rdf_interpolation;
}

Interpolate& PMF::GetPotentialInterpolation(){
    return this->potential_interpolation;
}

void PMF::ReadRDF(){

    std::string filename;
    std::vector<double> x_values;
    std::vector<double> y_values;

    filename = "rdf.txt";
    std::ifstream file{filename};
    if(!file){
        Log::global_log->error()<<"[PMF] I could not read the rdf data file"<<std::endl;
    }
    double n1, n2;

    while(file >> n1 >> n2){
        x_values.push_back(n1);
        y_values.push_back(n2);
    }

    reference_rdf_interpolation.SetXValues(x_values);
    reference_rdf_interpolation.SetYValues(y_values);

}

void PMF::MapToAtomistic(std::array<double,3> f, Molecule& m1, Molecule& m2){
    //does something
    double mass = m1.mass();
    
    //a loop for every type of charge
    for(int i=0;i<m1.numLJcenters();i++){
        double site_ratio = m1.component()->ljcenter(i).m()/mass;
        for(int j=0;j<f.size();j++){
            f[j] *= site_ratio;
        }
        m1.Fljcenteradd(i,f.data());
        m2.Fljcentersub(i,f.data());
    }

    for(int i=0;i<m1.numCharges();i++){
        double site_ratio = m1.component()->charge(i).m()/mass;
        for(int j=0;j<f.size();j++){
            f[j] *= site_ratio;
        }
        m1.Fchargeadd(i,f.data());
        m2.Fchargesub(i,f.data());
    }

    for(int i=0;i<m1.numQuadrupoles();i++){
        double site_ratio = m1.component()->quadrupole(i).m()/mass;
        for(int j=0;j<f.size();j++){
            f[j] *= site_ratio;
        }
        m1.Fquadrupoleadd(i,f.data());
        m2.Fquadrupolesub(i,f.data());
    }

    for(int i=0;i<m1.numDipoles();i++){
        double site_ratio = m1.component()->dipole(i).m()/mass;
        for(int j=0;j<f.size();j++){
            f[j] *= site_ratio;
        }
        m1.Fdipoleadd(i,f.data());
        m2.Fdipolesub(i,f.data());
    }

}

std::vector<double> PMF::GetAverageRDF(){
    std::vector<double> average_rdf = avg_rdf_interpolation.GetYValues();
    
    for(int i=0;i<average_rdf.size();++i){
        average_rdf[i] /= (double)profiler.GetMeasuredSteps();
    }

    return average_rdf;
}

double PMF::ConvergenceCheck(){
    std::vector<double> difference;
    std::vector<double>& v_0 = reference_rdf_interpolation.GetYValues();
    std::vector<double> v_i = GetAverageRDF();

    difference.resize(v_0.size());
    double vec_norm=0.0;

    for(int i=0;i<v_0.size();++i){
        vec_norm += std::pow(v_0[i]-v_i[i],2.0);
    }

    return vec_norm;

}

