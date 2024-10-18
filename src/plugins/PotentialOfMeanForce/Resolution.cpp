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








void ResolutionComponentHandler::init(){
    AddCGComponent();
    AddComponents2Ensemble();
}



void ResolutionComponentHandler::AddComponents2Ensemble(){
    _simulation.getEnsemble()->addComponent(cg);
}

void ResolutionComponentHandler::AddCGComponent(){
    double mass = _simulation.getEnsemble()->getComponent(1)->m();

    this->cg.addLJcenter(0,0,0,
                         mass,0,0,0);
    cg.setID(10);
    cg.setName("CG");

}

void ResolutionComponentHandler::CheckResolution(ParticleContainer* pc, std::vector<FPRegion>& regions){
    for(auto it = pc->iterator(ParticleIterator::ALL_CELLS);it.isValid();++it){
        
        bool stop=false;

        for(auto& reg:regions){
            if(reg.IsInsideResolutionRegion(it->r_arr(),FullParticle)){
                //Modify if needed
                CheckAndModifyMoleculeResolution(*it,FullParticle);
                stop=true;
                break;
            }
        }

        if(stop) continue;

        for(auto& reg:regions){
            if(reg.IsInsideResolutionRegion(it->r_arr(),Hybrid)){
                //Modify if needed
                CheckAndModifyMoleculeResolution(*it,Hybrid);
                stop=true;
                break;
            }
        }

        if(stop) continue;

        CheckAndModifyMoleculeResolution(*it, CoarseGrain);
    }
}

void ResolutionComponentHandler::CheckAndModifyMoleculeResolution(Molecule& mol,ResolutionType target){

    if(target == FullParticle){
        if(mol.componentid() != 1){
            mol.setComponent(_simulation.getEnsemble()->getComponent(1));
        }
    }

    if(target == CoarseGrain){
        if(mol.componentid() != 10){

            Coarsen(mol);
            
        }
    }

    if(target == Hybrid){
        Log::global_log->error()<<"[ResolutionHandler] Hybrid not implemented"<<std::endl;
    }

}

void ResolutionComponentHandler::Coarsen(Molecule& m){
    std::array<double, 3> com = ComputeCOM(m);
    m.setComponent(_simulation.getEnsemble()->getComponent(10));
    m.setr(0,com[0]);
    m.setr(1,com[1]);
    m.setr(2,com[2]);
}