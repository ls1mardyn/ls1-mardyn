#include"PMF.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"


PMF::PMF():reference_rdf_interpolation{1.0,false},potential_interpolation{0.0,true},avg_rdf_interpolation{1.0,false},convergence{0.5}{


}


void PMF::init(ParticleContainer* pc, DomainDecompBase* domainDecomp, Domain* domain){

    for(auto it= begin(resolution_handler->GetRegions());it!=end(resolution_handler->GetRegions());++it){
        Log::global_log->info()<<"[PMF] The hybrid region spans from: ("<<it->_lowHybrid[0]<<","<<it->_lowHybrid[1]<<","<<it->_lowHybrid[2]<<") to ("<<it->_highHybrid[0]<<","<<it->_highHybrid[1]<<","<<it->_highHybrid[2]<<")"<<std::endl;

        Log::global_log->info()<<"[PMF] The atomistic region spans from: ("<<it->_low[0]<<","<<it->_low[1]<<","<<it->_low[2]<<") to ("<<it->_high[0]<<","<<it->_high[1]<<","<<it->_high[2]<<")"<<std::endl;

        Log::global_log->info()<<"[PMF] The region  center is located at: ("<<it->_center[0]<<","<<it->_center[1]<<","<<it->_center[2]<<")"<<std::endl;
    }



    adres_statistics.init(resolution_handler->GetRegions()[0]);
    Log::global_log->info()<<"[PMF] Statistical tool initialized"<<std::endl;

}

void PMF::readXML(XMLfileUnits& xmlfile){

    resolution_handler = new ResolutionComponentHandler();

    //create AT regions
    int num_regions =0;
    XMLfile::Query query = xmlfile.query("fpregions/region");
    num_regions = query.card();

    resolution_handler->GetRegions().resize(num_regions);

    XMLfile::Query::const_iterator region_iterator;
    std::string oldpath = xmlfile.getcurrentnodepath();
    for(region_iterator = query.begin();region_iterator;region_iterator++){
        xmlfile.changecurrentnode(region_iterator);
        unsigned int id=0;
        xmlfile.getNodeValue("@id",id);
        resolution_handler->GetRegions()[id-1].readXML(xmlfile);
    }
    xmlfile.changecurrentnode(oldpath);
    xmlfile.getNodeValue("multiplier",multiplier);
    xmlfile.getNodeValue("output",output);
    xmlfile.getNodeValue("outputFreq",output_freq);
    xmlfile.getNodeValue("updateStride",update_stride);
    xmlfile.getNodeValue("rdfPath",ref_rdf_path);
    xmlfile.getNodeValue("epPath",effective_potential_path);
    xmlfile.getNodeValue("doNothing",do_nothing);
    if(do_nothing){
        Log::global_log->info()<<"[PMF] Do Nothing set"<<std::endl;
    }

    std::string mode_name;
    xmlfile.getNodeValue("mode",mode_name);
    if(do_nothing && mode_name != "Production"){
        Log::global_log->info()<<"[PMF] Warning, do nothing but not in production mode, exiting"<<std::endl;
        Simulation::exit(1);
    }
    if(mode_name == "Equilibration"){
        mode = Mode::Equilibration;
        Log::global_log->info()<<"[PMF] Mode set to Equilibration"<<std::endl;
    }
    if(mode_name == "EffectivePotential"){
        mode = Mode::EffectivePotential;
        Log::global_log->info()<<"[PMF] Mode set to EffectivePotential"<<std::endl;
        this->ReadRDF();
        this->SetPotentialInitialGuess();
        Log::global_log->info()<<"[PMF] PMF set as initial guess\n";
        Log::global_log->info()<<"[PMF] RDF has been read successfully\n";
    }
    if(mode_name == "Production"){
        profiler.SetMeasureDensity(true);
        mode = Mode::Production;
        Log::global_log->info()<<"[PMF] Mode set to Production"<<std::endl;
        this->ReadEffectivePotential();
        Log::global_log->info()<<"[PMF] Effective potential has been read successfully\n";
    }


    profiler.ReadXML(xmlfile);
    Log::global_log->info()<<"[PMF] Target temperature = "<<_simulation.getEnsemble()->T()<<std::endl;    


    this->profiler.init(_simulation.getcutoffRadius());
    Log::global_log->info()<<"[PMF] InternalProfiler initialized\n";
    Log::global_log->info()<<"[PMF] Alpha value "<<multiplier<<"\n";

    if(mode == Mode::EffectivePotential){
        potential_interpolation.SetXValues(profiler.GetBinCenters());
        this->SetPotentialInitialGuess();
    }

    avg_rdf_interpolation.ResizeVectors(profiler.GetBinCenters().size());
    avg_rdf_interpolation.SetXValues(profiler.GetBinCenters());

    resolution_handler->init();
    Log::global_log->info()<<"[PMF]Initialized resolution manager\n";

    pairs_handler = new InteractionForceAdapter(resolution_handler,this);

    _simulation.setParticlePairsHandler(pairs_handler);
    _simulation.useLegacyCellProcessor();
    Log::global_log->info()<<"[PMF] InteractionForcedAdapter Class being used\n";
    _simulation.setCellProcessor(new LegacyCellProcessor(_simulation.getcutoffRadius(), _simulation.getLJCutoff(), pairs_handler));
    Log::global_log->info()<<"[PMF]LegacyCellProcessor set\n";



}

void PMF::beforeForces(ParticleContainer* pc, DomainDecompBase* domainDecomp, unsigned long simstep){

    resolution_handler->CheckContainerResolution(pc);



}

void PMF::beforeEventNewTimestep(ParticleContainer* pc, DomainDecompBase* domainDecomp, unsigned long simstep){

}

void PMF::afterForces(ParticleContainer* pc, DomainDecompBase* dd, unsigned long step){
    if(!do_nothing){
        profiler.ProfileData(pc,step);
        if(mode==Mode::Production){
            adres_statistics.MeasureLocalRDFs(pc,step);
        }  
    }
    
}

void PMF::endStep(ParticleContainer* pc, DomainDecompBase* dd, Domain* domain, unsigned long step){

    if(mode == Mode::EffectivePotential && !do_nothing){
        UpdateRDFInterpolation();
    }


    if(mode == Mode::Production && !do_nothing){
        adres_statistics.MeasureStatistics(pc);
        adres_statistics.Output2File(step); 
        adres_statistics.PrintRDFs2File(step, pc);
    }



    if(output && !do_nothing){

        profiler.PrintOutput2Files(step);

    }
    
    if(mode == Mode::EffectivePotential){

        std::vector<double> current_rdf = profiler.GetRDF();
        convergence.PrintGlobalConvergence2File();
        convergence.CheckConvergence(reference_rdf_interpolation.GetYValues(),current_rdf);
        convergence.PrintLocalConvergence2File();
    
        if(step%update_stride ==0 && step>0){
            Log::global_log->info()<<"[UpdatePotential]Update potential now"<<std::endl;
            convergence.PrepareUpdate();
            AddPotentialCorrection(step);
            
        }

        if(output && step%update_stride ==0){
            
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

    }
}

void PMF::siteWiseForces(ParticleContainer* pc, DomainDecompBase* dd, unsigned long step){

}
/********************
 * ****************** FUNCTIONS NOT FROM THE INTERFACE
 *******************/


void PMF::SetPotentialInitialGuess(){
    std::vector<double> pot0 = reference_rdf_interpolation.GetYValues();
    for(int i=0;i<pot0.size();++i){
        double val = pot0[i];
        if(val<0.01){
            val=0;
        }
        pot0[i] = -1.0*_simulation.getEnsemble()->T()*std::log(val);
    }

    potential_interpolation.SetYValues(pot0);
    FivePointAverageExtrapolation(potential_interpolation.GetXValues(),potential_interpolation.GetYValues());
}


void PMF::PrintPotentialCorrection(unsigned long step){
    std::vector<double> avg_rdf = profiler.GetRDF();
    int total_data_points = avg_rdf.size();
    std::vector<double> pot_i;
    std::vector<double> current_correction;
    current_correction.resize(total_data_points,0.0);
    pot_i.resize(total_data_points,0.0);

    for(int i=0;i<total_data_points;++i){
        double ratio;
        if(reference_rdf_interpolation.GetYValues()[i]==0){
            ratio=0;
        }else{
            ratio = avg_rdf[i]/reference_rdf_interpolation.GetYValues()[i];
        }
        
        current_correction[i] = 1.0*multiplier* _simulation.getEnsemble()->T()*std::log(ratio);
    }

    std::string name="correction_step_"+std::to_string(step)+".txt";
    std::ofstream corr{name};
    for(int i=0;i<total_data_points;++i){
        corr<<reference_rdf_interpolation.GetXValues()[i]
            <<"\t"
            <<current_correction[i]
            <<std::endl;
    }


}

void PMF::AddPotentialCorrection(unsigned long step){
    std::vector<double> avg_rdf = profiler.GetRDF();
    int total_data_points = avg_rdf.size();
    std::vector<double> pot_i;
    std::vector<double> current_correction;
    current_correction.resize(total_data_points,0.0);
    pot_i.resize(total_data_points,0.0);

    for(int i=0;i<total_data_points;++i){
        double ratio = avg_rdf[i]/reference_rdf_interpolation.GetYValues()[i];
        if(avg_rdf[i]<0.01){
            ratio=0;
        }
        current_correction[i] = -1.0*multiplier* _simulation.getEnsemble()->T()*std::log(ratio);
    }

    VectorAdd(potential_interpolation.GetYValues(),current_correction);
    // VectorSub(potential_interpolation.GetYValues(),current_correction);
    if(output){
        std::string name="correction_step_"+std::to_string(step)+".txt";
        std::ofstream corr{name};
        for(int i=0;i<total_data_points;++i){
            corr<<reference_rdf_interpolation.GetXValues()[i]
                <<"\t"
                <<current_correction[i]
                <<std::endl;
        }

    }

}

void PMF::UpdateRDFInterpolation(){
    std::vector<double> rdf_i = this->profiler.GetRDF();
    avg_rdf_interpolation.SetYValues(rdf_i);
}

double PMF::WeightValue(const std::array<double,3>& pos, FPRegion& region){
    return weight_function.WeightValue(pos,region);
}

std::vector<FPRegion>& PMF::GetRegions(){
    return resolution_handler->GetRegions();
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
    std::ifstream file{ref_rdf_path};
    if(!file){
        Log::global_log->error()<<"[PMF] I could not read the rdf data file"<<std::endl;
        Simulation::exit(1);
    }
    double n1, n2;

    while(file >> n1 >> n2){
        x_values.push_back(n1);
        y_values.push_back(n2);
    }

    reference_rdf_interpolation.SetXValues(x_values);
    reference_rdf_interpolation.SetYValues(y_values);

}

void PMF::ReadEffectivePotential(){

    std::vector<double> x_values;
    std::vector<double> y_values;

    std::ifstream file{effective_potential_path};
    if(!file){
        Log::global_log->error()<<"[PMF] I could not read the potential data file"<<std::endl;
        Simulation::exit(1);
    }
    double n1, n2;

    while(file >> n1 >> n2){
        x_values.push_back(n1);
        y_values.push_back(n2);
    }

    potential_interpolation.SetXValues(x_values);
    potential_interpolation.SetYValues(y_values);

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
