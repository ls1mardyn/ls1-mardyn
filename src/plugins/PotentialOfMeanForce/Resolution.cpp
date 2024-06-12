#include "Resolution.h"
#include "PMF.h"

void ResolutionHandler::CheckResolution(ParticleContainer* pc, std::map<unsigned long, std::pair<InteractionSite,ResolutionType>>& resolution_map, std::vector<FPRegion>& regions){
    //check all molecules
    for(auto it = pc->iterator(ParticleIterator::ALL_CELLS);it.isValid();++it){
        unsigned long molecule_id = it->getID();
        bool stop=false;

        for(auto& reg:regions){
            if(reg.isInnerPointDomain(_simulation.getDomain(),FullParticle,it->r_arr())){
                //Modify if needed
                CheckAndModifyMoleculeResolution(resolution_map[molecule_id],FullParticle);
                stop=true;
                break;
            }
        }

        if(stop) continue;

        for(auto& reg:regions){
            if(reg.isInnerPointDomain(_simulation.getDomain(),Hybrid,it->r_arr())){
                //Modify if needed
                CheckAndModifyMoleculeResolution(resolution_map[molecule_id],Hybrid);
                stop=true;
                break;
            }
        }


    }
}

void ResolutionHandler::CheckAndModifyMoleculeResolution(std::pair<InteractionSite,ResolutionType>& val, ResolutionType res_type){

    if(val.second == res_type){
        return;
    }

    val.second=res_type;
}