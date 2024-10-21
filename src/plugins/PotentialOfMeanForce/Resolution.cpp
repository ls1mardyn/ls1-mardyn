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
    // AddComponents2Ensemble();
}


ResolutionType ResolutionComponentHandler::GetMoleculeResolution(const Molecule& m){
    if(m.component()->getName() == "CG")
        return ResolutionType::CoarseGrain;
    
    return ResolutionType::FullParticle;
}



void ResolutionComponentHandler::AddComponents2Ensemble(){
    _simulation.getEnsemble()->addComponent(cg);
}

void ResolutionComponentHandler::AddCGComponent(){
    double mass = _simulation.getEnsemble()->getComponent("LJ")->m();
    std::cout<<"Mass "<<mass<<std::endl;
    this->cg.addLJcenter(0,0,0,
                         mass,0,0,
                         0);
    cg.setID(5);
    cg.setName("CG");

    _simulation.getEnsemble()->addComponent(cg);
    std::ostream& out = std::cout;
    _simulation.getEnsemble()->getComponent("CG")->write(out);

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
            mol.setComponent(_simulation.getEnsemble()->getComponent("LJ"));
        }
    }

    if(target == CoarseGrain){
        if(mol.component()->getName() != "CG"){

            Coarsen(mol);
            
        }
    }

    if(target == Hybrid){
        Log::global_log->error()<<"[ResolutionHandler] Hybrid not implemented"<<std::endl;
    }

}

void ResolutionComponentHandler::Coarsen(Molecule& mol){
    std::array<double, 3> com = ComputeCOM(mol);
    mol.setComponent(_simulation.getEnsemble()->getComponent("CG"));
    mol.setr(0,com[0]);
    mol.setr(1,com[1]);
    mol.setr(2,com[2]);
}