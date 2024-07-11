#include "Resolution.h"
#include "PMF.h"

void ResolutionHandler::CheckResolution(ParticleContainer* pc, std::map<unsigned long, std::pair<InteractionSite,ResolutionType>>& resolution_map, std::vector<FPRegion>& regions){
    //check all molecules
    for(auto it = pc->iterator(ParticleIterator::ALL_CELLS);it.isValid();++it){
        unsigned long molecule_id = it->getID();
        //std::cout<<"Molecule "<<molecule_id<<" has resolution "<<resolution_map[molecule_id].second<<" \n";
        bool stop=false;

        for(auto& reg:regions){
            if(reg.IsInsideResolutionRegion(resolution_map[molecule_id].first.r(),FullParticle)){
                //Modify if needed
                CheckAndModifyMoleculeResolution(resolution_map[molecule_id],FullParticle);
                stop=true;
                break;
            }
        }

        if(stop) continue;

        for(auto& reg:regions){
            if(reg.IsInsideResolutionRegion(resolution_map[molecule_id].first.r(),Hybrid)){
                //Modify if needed
                CheckAndModifyMoleculeResolution(resolution_map[molecule_id],Hybrid);
                stop=true;
                break;
            }
        }

        if(stop) continue;

        CheckAndModifyMoleculeResolution(resolution_map[molecule_id], CoarseGrain);


    }
}

void ResolutionHandler::CheckAndModifyMoleculeResolution(std::pair<InteractionSite,ResolutionType>& tracker, ResolutionType res_type){

    if(tracker.second == res_type){
        return;
    }

    tracker.second=res_type;
}