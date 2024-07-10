#include"PMF.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"


PMF::PMF(){
    adres_cell_processor = new InteractionCellProcessor(0,0);
}


void PMF::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain){
    
    Log::global_log->info()<<"[PMF] Enabled "<<std::endl;
    this->ReadRDF();
    for(int i=0;i<rdf_interpolation.GetGValues().size();i++){
        std::cout<<rdf_interpolation.GetRValues()[i]<<"       "<<rdf_interpolation.GetGValues()[i]<<"\n";
    }

    pairs_handler = new InteractionForceAdapter(resolution_handler,this);
    _simulation.setParticlePairsHandler(pairs_handler);
    _simulation.setCellProcessor(new LegacyCellProcessor(_simulation.getcutoffRadius(), _simulation.getLJCutoff(), pairs_handler));

    Log::global_log->info()<<"[PMF]LegacyCellProcessor set\n";
    Log::global_log->info()<<"[PMF] AdResS ParticlePairsHandler being used\n";
    
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

void PMF::beforeEventNewTimestep(ParticleContainer* pc, DomainDecompBase* domainDecomp, unsigned long simstep)
{
    resolution_handler.CheckResolution(pc,sites,regions);

    for(auto it= pc->iterator(ParticleIterator::ALL_CELLS);it.isValid();++it){
        unsigned long m_id = it->getID();
        std::array<double,3> com = rdf.GetCOM(&(*it));
        sites[m_id].first.SetPosition(com);

        std::cout<<m_id<<"   "<<sites[m_id].second<<"\n";
    }
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
    return this->rdf_interpolation;
}

void PMF::ReadRDF(){

    rdf_interpolation.ReadInRDF();

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

