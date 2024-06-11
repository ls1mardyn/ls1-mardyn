#include"PMF.h"

PMF::PMF(){
    cell_processor = new InteractionCellProcessor(0,0);
}

void PMF::init(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain){
    
    Log::global_log->info()<<"[PMF] Enabled "<<std::endl;
    Log::global_log->info()<<"[PMF] Atomistic region: "<<std::endl;
    for(int i=0;i<regions.size();i++){
        Log::global_log->info()<<regions[i]._center[0]<<std::endl;
        Log::global_log->info()<<regions[i]._low[0]<<" "<<regions[i]._low[1]<<" "<<regions[i]._low[2]<<std::endl;
        Log::global_log->info()<<regions[i]._high[0]<<" "<<regions[i]._high[1]<<" "<<regions[i]._high[2]<<std::endl;
    }

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