#include"PMF.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"


PMF::PMF():reference_rdf_interpolation{1.0},current_rdf_interpolation{1.0},potential_interpolation{0.0}{
    adres_cell_processor = new InteractionCellProcessor(0,0);
}


void PMF::init(ParticleContainer* pc, DomainDecompBase* domainDecomp, Domain* domain){
    
    this->ReadRDF();
    Log::global_log->info()<<"[PMF] RDF has been read successfully\n.";
    pairs_handler = new InteractionForceAdapter(resolution_handler,this);
    _simulation.setParticlePairsHandler(pairs_handler);
    _simulation.setCellProcessor(new LegacyCellProcessor(_simulation.getcutoffRadius(), _simulation.getLJCutoff(), pairs_handler));

    Log::global_log->info()<<"[PMF]LegacyCellProcessor set\n";
    Log::global_log->info()<<"[PMF] ForcedAdapter Class being used\n";

    Log::global_log->info()<<"[PMF] Start the tracker sites\n";
    for(auto it= pc->iterator(ParticleIterator::ALL_CELLS);it.isValid();++it){
        unsigned long m_id = it->getID();
        std::array<double,3> com = rdf.GetCOM(&(*it));
        sites[m_id].first.SetPosition(com);
    }

    for(auto it= begin(regions);it!=end(regions);++it){
        Log::global_log->info()<<"[PMF] The hybrid region spans from: ("<<it->_lowHybrid[0]<<","<<it->_lowHybrid[1]<<","<<it->_lowHybrid[2]<<") to ("<<it->_highHybrid[0]<<","<<it->_highHybrid[1]<<","<<it->_highHybrid[2]<<")"<<std::endl;

        Log::global_log->info()<<"[PMF] The atomistic region spans from: ("<<it->_low[0]<<","<<it->_low[1]<<","<<it->_low[2]<<") to ("<<it->_high[0]<<","<<it->_high[1]<<","<<it->_high[2]<<")"<<std::endl;

        Log::global_log->info()<<"[PMF] The regio  center is located at: ("<<it->_center[0]<<","<<it->_center[1]<<","<<it->_center[2]<<")"<<std::endl;
    }
    Log::global_log->info()<<"[PMF] Initializing the COM sites\n";
    resolution_handler.CheckResolution(pc,sites,regions);
    Log::global_log->info()<<"[PMF] Enabled "<<std::endl;
    current_rdf_interpolation.SetXValues(profiler.GetRNodes());
    potential_interpolation.SetXValues(profiler.GetRNodes());
    this->profiler.init(pc,200,1);
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



}

void PMF::beforeEventNewTimestep(ParticleContainer* pc, DomainDecompBase* domainDecomp, unsigned long simstep){

    for(auto it= pc->iterator(ParticleIterator::ALL_CELLS);it.isValid();++it){
        unsigned long m_id = it->getID();
        std::array<double,3> com = rdf.GetCOM(&(*it));
        sites[m_id].first.SetPosition(com);
    }

    resolution_handler.CheckResolution(pc,sites,regions);
    //on every time step
    //transfer buffers to interpolation
    profiler.ResetBuffers();
}

void PMF::siteWiseForces(ParticleContainer* pc, DomainDecompBase* dd, unsigned long step){

}
/********************
 * ****************** FUNCTIONS NOT FROM THE INTERFACE
 *******************/

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

void PMF::ReadRDF(){

    reference_rdf_interpolation.ReadInRDF();

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

