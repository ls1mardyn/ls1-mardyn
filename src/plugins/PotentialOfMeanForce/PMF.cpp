#include"PMF.h"
#include "particleContainer/adapter/VectorizedCellProcessor.h"

PMF::PMF(){
    adres_cell_processor = new InteractionCellProcessor(0,0);
}


void PMF::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain){
    
    Log::global_log->info()<<"[PMF] Enabled "<<std::endl;

    this->adres_cell_processor->SetCellProcessor(new VectorizedCellProcessor(*domain,_simulation.getcutoffRadius(),_simulation.getLJCutoff()));
    _simulation.setCellProcessor(adres_cell_processor);
    Log::global_log->info()<<"[PMF] Initialized AdResS CellProcessor"<<std::endl;

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
}
/********************
 * ****************** FUNCTIONS NOT FROM THE INTERFACE
 *******************/

double PMF::WeightValue(std::array<double,3>& pos, FPRegion& region){
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
