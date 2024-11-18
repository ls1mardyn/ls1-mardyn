#include "Resolution.h"
#include "PMF.h"

void ResolutionHandler::SetMoleculeTrackerPosition(unsigned long idx, std::array<double, 3>& pos){
    sites[idx].first.SetPosition(pos);
}

std::array<double, 3>& ResolutionHandler::GetMoleculeTrackerPosition(unsigned long idx){
    return sites[idx].first.GetPosition();
}

ResolutionHandler::cg_map& ResolutionHandler::GetCGMap(){
    return sites;
}

ResolutionType ResolutionHandler::GetMoleculeResolution(unsigned long idx){
    return sites[idx].second;
}

InteractionSite ResolutionHandler::GetMoleculeCOMSite(unsigned long idx){
    return sites[idx].first;
}

ResolutionType ResolutionHandler::GetCOMResolution(std::array<double,3>& com, std::vector<FPRegion>& regions){

    for(auto& reg:regions){
        if(reg.IsInsideResolutionRegion(com,FullParticle)){
            return FullParticle;
        }

        if(reg.IsInsideResolutionRegion(com,Hybrid)){
            return Hybrid;
        }

        return CoarseGrain;

    }

}


void ResolutionHandler::CheckResolution(ParticleContainer* pc, std::vector<FPRegion>& regions){
    //check all molecules
    for(auto it = pc->iterator(ParticleIterator::ALL_CELLS);it.isValid();++it){
        unsigned long molecule_id = it->getID();
        //std::cout<<"Molecule "<<molecule_id<<" has resolution "<<resolution_map[molecule_id].second<<" \n";
        bool stop=false;

        for(auto& reg:regions){
            if(reg.IsInsideResolutionRegion(sites[molecule_id].first.GetPosition(),FullParticle)){
                //Modify if needed
                CheckAndModifyMoleculeResolution(sites[molecule_id],FullParticle);
                stop=true;
                break;
            }
        }

        if(stop) continue;

        for(auto& reg:regions){
            if(reg.IsInsideResolutionRegion(sites[molecule_id].first.GetPosition(),Hybrid)){
                //Modify if needed
                CheckAndModifyMoleculeResolution(sites[molecule_id],Hybrid);
                stop=true;
                break;
            }
        }

        if(stop) continue;

        CheckAndModifyMoleculeResolution(sites[molecule_id], CoarseGrain);


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
}


ResolutionType ResolutionComponentHandler::GetMoleculeResolution(const Molecule& m){
    if(m.component()->getName() == "CG")
        return ResolutionType::CoarseGrain;
    
    return ResolutionType::FullParticle;
}

void ResolutionComponentHandler::AddCGComponent(){
    cg.setID(2);
    cg.setName("CG");
    LJcenter single_site{};
    single_site.setEps(0);
    single_site.setSigma(1.0);
    std::string site_name{"IBI"};
    single_site.setName(site_name);
    single_site.setR(0,0.0);
    single_site.setR(1,0.0);
    single_site.setR(2,0.0);

    //sum all mass from the multi-site component
    Component* fp = _simulation.getEnsemble()->getComponent("LJ");
    double cg_mass = 0.0;
    for(int i=0;i<fp->numLJcenters();++i){
        cg_mass += fp->ljcenter(i).m();
    }
    single_site.setM(cg_mass);
    cg.addLJcenter(single_site);

    // _simulation.getEnsemble()->getComponents()->push_back(cg);
    // _simulation.getEnsemble()->setComponentLookUpIDs();//Why call this?

    // for(auto it=_simulation.getMoleculeContainer()->iterator(ParticleIterator::ALL_CELLS);it.isValid();++it){
    //     it->setComponent(fp);
    // }

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
    ResolutionType mol_resolution = GetMoleculeResolution(mol);
    Component* component;
    if(target != mol_resolution ){
        if(target == ResolutionType::FullParticle){
            component = _simulation.getEnsemble()->getComponent("LJ");
            mol.setComponent(component);
        }

        if(target == ResolutionType::CoarseGrain){
            
            std::array<double, 3> com{0,0,0};
            double total_mass = mol.component()->m();

            for(int lj=0;lj<mol.component()->numLJcenters();++lj){
                auto lj_site = mol.ljcenter_d_abs(lj);
                double site_mass = mol.component()->ljcenter(lj).m();
                for(int i=0;i<3;++i){
                    com[i] += lj_site[i]*site_mass;
                    com[i] = com[i]/total_mass;
                }
            }
            component = &cg;
            mol.setComponent(component);
        }
    }

}

void ResolutionComponentHandler::Coarsen(Molecule& mol){
    //TODO: implement
}