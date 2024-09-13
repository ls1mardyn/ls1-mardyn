#include"PMF.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"


PMF::PMF():reference_rdf_interpolation{1.0},current_rdf_interpolation{1.0},potential_interpolation{0.0},avg_rdf_interpolation{1.0}{
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
    current_rdf_interpolation.SetXValues(profiler.GetRNodes());
    potential_interpolation.SetXValues(profiler.GetRNodes());
    avg_rdf_interpolation.SetXValues(profiler.GetRNodes());

    accumulate_rdf_buffer.resize(internal_bins);

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
    xmlfile.getNodeValue_double("multiplier",multiplier);
    xmlfile.getNodeValue("internalBins",internal_bins);
    xmlfile.getNodeValue("measureFreq",measure_frequency);
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

    //Generates instantaneous rdf values
    this->profiler.GenerateInstantaneousData(pc,domain);
    //Accumulates values
    AccumulateRDF(profiler.GetRDFValues());
    std::vector<double> avg_rdf = GetAverageRDF();
    //Pass to interpolation
    avg_rdf_interpolation.SetYValues(avg_rdf);

    if(global_simulation->getSimulationStep()>0){
        Log::global_log->info()<<"[PMF] Convergence check: "<<ConvergenceCheck()<<std::endl;
        std::string filename="avg_rdf_"+std::to_string(step)+".txt";
        std::ofstream rdf_file(filename);
        for(int i=0;i<avg_rdf_interpolation.GetYValues().size();++i){
            rdf_file<<std::setw(8)<<std::left<<avg_rdf_interpolation.GetXValues()[i]<<"\t"<<std::setw(8)<<std::left<<avg_rdf_interpolation.GetYValues()[i]<<std::endl;
        }
        rdf_file.close();
    }

    //need to check if correction required
    AddPotentialCorrection();
    profiler.ResetBuffers();

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

void PMF::AddPotentialCorrection(){
    std::vector<double> pot_i = potential_interpolation.GetYValues();
    for(int i=0;i<pot_i.size();++i){
        double correction = std::log(avg_rdf_interpolation.GetYValues()[i])/std::log(reference_rdf_interpolation.GetYValues()[i]);
        pot_i[i] +=  -1.0*_simulation.getEnsemble()->T()*correction;
    }
    potential_interpolation.SetYValues(pot_i);
}

void PMF::AccumulateRDF(std::vector<double>& current_rdf){

    for(int i=0;i<internal_bins;++i){
        accumulate_rdf_buffer[i] += current_rdf[i];
    }

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

Interpolate& PMF::GetCurrentRDFInterpolation(){
    return this->current_rdf_interpolation;
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
    std::vector<double> average_rdf = accumulate_rdf_buffer;
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

