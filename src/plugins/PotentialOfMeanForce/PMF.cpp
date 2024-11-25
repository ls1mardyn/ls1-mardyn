#include"PMF.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"


PMF::PMF():reference_rdf_interpolation{1.0,false},potential_interpolation{0.0,true},acc_rdf_interpolation{1.0,false},convergence{0.5}{


}


void PMF::init(ParticleContainer* pc, DomainDecompBase* domainDecomp, Domain* domain){

    pairs_handler = new InteractionForceAdapter(resolution_handler,component_handler,this);
    _simulation.setParticlePairsHandler(pairs_handler);
    Log::global_log->info()<<"[PMF] InteractionForcedAdapter Class being used\n";
    _simulation.setCellProcessor(new LegacyCellProcessor(_simulation.getcutoffRadius(), _simulation.getLJCutoff(), pairs_handler));
    Log::global_log->info()<<"[PMF]LegacyCellProcessor set\n";

    for(auto it= pc->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);it.isValid();++it){
        unsigned long m_id = it->getID();
        std::array<double,3> com = ComputeCOM(*it);
        resolution_handler.SetMoleculeTrackerPosition(m_id,com);
    }


     Log::global_log->info()<<"[PMF] Started the tracker sites\n";

    for(auto it= begin(regions);it!=end(regions);++it){
        Log::global_log->info()<<"[PMF] The hybrid region spans from: ("<<it->_lowHybrid[0]<<","<<it->_lowHybrid[1]<<","<<it->_lowHybrid[2]<<") to ("<<it->_highHybrid[0]<<","<<it->_highHybrid[1]<<","<<it->_highHybrid[2]<<")"<<std::endl;

        Log::global_log->info()<<"[PMF] The atomistic region spans from: ("<<it->_low[0]<<","<<it->_low[1]<<","<<it->_low[2]<<") to ("<<it->_high[0]<<","<<it->_high[1]<<","<<it->_high[2]<<")"<<std::endl;

        Log::global_log->info()<<"[PMF] The region  center is located at: ("<<it->_center[0]<<","<<it->_center[1]<<","<<it->_center[2]<<")"<<std::endl;
    }
    resolution_handler.CheckResolution(pc,regions);
    Log::global_log->info()<<"[PMF] Resolutions started "<<std::endl;


    this->profiler.init(pc,internal_bins,measure_frequency,density_bins);
    Log::global_log->info()<<"[PMF] InternalProfiler uses "<<internal_bins
                           <<" bins, measures every "<<measure_frequency<<" steps\n";
    Log::global_log->info()<<"[PMF] InternalProfiler initialized\n";
    Log::global_log->info()<<"[PMF] Alpha value "<<multiplier<<"\n";



    potential_interpolation.SetXValues(profiler.GetBinCenters());
    acc_rdf_interpolation.SetXValues(profiler.GetBinCenters());
    acc_rdf_interpolation.GetYValues().resize(internal_bins);

    this->ReadRDF();
    Log::global_log->info()<<"[PMF] RDF has been read successfully\n";

    if(mode == Mode::Production){
        this->ReadEffectivePotential();
        Log::global_log->info()<<"[PMF] Effective potential has been read successfully\n";
    }

    if(mode ==Mode::OnlyCG){
        this->InitializePotentialValues();
        Log::global_log->info()<<"[PMF] PMF computed as initial guess\n";
    }

    component_handler.init();
    Log::global_log->info()<<"[PMF] Component handler initialized"<<std::endl;

    adres_statistics.init(regions[0]);
    Log::global_log->info()<<"[PMF] Statistical tool initialized"<<std::endl;

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
    xmlfile.getNodeValue("densityBins",density_bins);
    xmlfile.getNodeValue("measureFreq",measure_frequency);
    xmlfile.getNodeValue("output",output);
    xmlfile.getNodeValue("rdfOutputStride",output_frequency);
    xmlfile.getNodeValue("updateStride",update_stride);
    std::string mode_name;
    xmlfile.getNodeValue("mode",mode_name);
    if(mode_name == "Equilibration"){
        mode = Mode::Equilibration;
    }
    if(mode_name == "OnlyCG"){
        mode = Mode::OnlyCG;
    }
    if(mode_name == "Production"){
        mode = Mode::Production;
    }
    Log::global_log->info()<<"[PMF] Target temperature = "<<_simulation.getEnsemble()->T()<<std::endl;
    if(mode == Mode::Production){
        Log::global_log->info()<<"[PMF] Mode set to Production"<<std::endl;
    }
}

void PMF::beforeForces(ParticleContainer* pc, DomainDecompBase* domainDecomp, unsigned long simstep){
    // component_handler.CheckResolution(pc,regions);
}

void PMF::beforeEventNewTimestep(ParticleContainer* pc, DomainDecompBase* domainDecomp, unsigned long simstep){

    // component_handler.CheckResolution(pc,regions);
    
    for(auto it= pc->iterator(ParticleIterator::ALL_CELLS);it.isValid();++it){
        unsigned long m_id = it->getID();
        std::array<double,3> com = ComputeCOM(*it);
        resolution_handler.SetMoleculeTrackerPosition(m_id,com);
    }

    resolution_handler.CheckResolution(pc,regions);
}

void PMF::afterForces(ParticleContainer* pc, DomainDecompBase* dd, unsigned long step){
    profiler.ProfileData(pc,step);
}

void PMF::endStep(ParticleContainer* pc, DomainDecompBase* dd, Domain* domain, unsigned long step){
    adres_statistics.MeasureStatistics(pc);
    adres_statistics.Output();
    AccumulateRDF(pc,domain);

    if(step > 0 && step%output_frequency==0){
        std::string filename="avg_rdf_"+std::to_string(step)+".txt";
        std::ofstream rdf_file(filename);
        for(int i=0;i<acc_rdf_interpolation.GetYValues().size();++i){
            rdf_file<<std::setw(8)<<std::left<<acc_rdf_interpolation.GetXValues()[i]<<"\t"<<std::setw(8)<<std::left<<GetAverageRDF()[i]<<std::endl;
        }
        rdf_file.close();

        filename = "density_"+std::to_string(step)+".txt";
        std::ofstream density_file(filename);
        for(int i=0;i<profiler.GetDensity(domain).size();++i){
            density_file<<std::setw(8)<<std::left<<profiler.GetDensityBinCenters()[i]<<"\t"<<std::setw(8)<<std::left<<profiler.GetDensity(domain)[i]<<std::endl;
        }

    }
    
    if(mode == Mode::OnlyCG){
    Log::global_log->info()<<"[UpdatePotential]Computing effective potential"<<std::endl;
    std::vector<double> current_rdf = GetAverageRDF();

    convergence.PrintGlobalConvergence2File();
    convergence.CheckConvergence(reference_rdf_interpolation.GetYValues(),current_rdf);
    convergence.PrintLocalConvergence2File(step);

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


void PMF::InitializePotentialValues(){
    std::vector<double> pot0 = reference_rdf_interpolation.GetYValues();
    for(int i=0;i<pot0.size();++i){
        pot0[i] = -1.0*_simulation.getEnsemble()->T()*std::log(pot0[i]);
    }

    potential_interpolation.SetYValues(pot0);
    potential_interpolation.LinearExtrapolation();
}

void PMF::AddPotentialCorrection(unsigned long step){
    std::vector<double> pot_i;
    std::vector<double> avg_rdf = GetAverageRDF();
    std::vector<double> current_correction;
    current_correction.resize(internal_bins);
    pot_i.resize(internal_bins);

    for(int i=0;i<internal_bins;++i){
        double ratio = avg_rdf[i]/reference_rdf_interpolation.GetYValues()[i];
        current_correction[i] = 1.0*multiplier* _simulation.getEnsemble()->T()*std::log(ratio);
    }

    potential_interpolation.AddVector(current_correction);

    if(output && step%update_stride ==0 ){
        std::string name="correction_step_"+std::to_string(step)+".txt";
        std::ofstream corr{name};
        for(int i=0;i<internal_bins;++i){
            corr<<reference_rdf_interpolation.GetXValues()[i]
                <<"\t"
                <<current_correction[i]
                <<std::endl;
        }

    }

}

void PMF::AccumulateRDF(ParticleContainer* pc, Domain* dom){
    std::vector<double> rdf_i = this->profiler.GetInstantaneousData(dom);
    std::vector<double>& accumulated_rdf = acc_rdf_interpolation.GetYValues();

    for(int i=0;i<rdf_i.size();++i){

        accumulated_rdf[i] += rdf_i[i];

    }
    acc_rdf_interpolation.SetYValues(accumulated_rdf);
}

double PMF::WeightValue(const std::array<double,3>& pos, FPRegion& region){
    return weight_function.WeightValue(pos,region);
}

std::vector<FPRegion>& PMF::GetRegions(){
    return this->regions;
}

Interpolate& PMF::GetRDFInterpolation(){
    return this->reference_rdf_interpolation;
}

Interpolate& PMF::GetAVGRDFInterpolation(){
    return this->acc_rdf_interpolation;
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

void PMF::ReadEffectivePotential(){

    std::string filename;
    std::vector<double> x_values;
    std::vector<double> y_values;

    filename = "cg_potential.txt";
    std::ifstream file{filename};
    if(!file){
        Log::global_log->error()<<"[PMF] I could not read the potential data file"<<std::endl;
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

std::vector<double> PMF::GetAverageRDF(){
    std::vector<double> average_rdf = acc_rdf_interpolation.GetYValues();
    
    for(int i=0;i<average_rdf.size();++i){
        average_rdf[i] /= ((double)profiler.GetMeasuredSteps());
    }

    return average_rdf;
}
