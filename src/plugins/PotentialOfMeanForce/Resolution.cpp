#include "Resolution.h"

void ResolutionHandler::init(){
    //left empty
}

ResolutionType ResolutionHandler::GetMoleculeResolution(Molecule& mol){
    std::array<double, 3> com = ComputeCOM(mol);

    return GetCOMResolution(com, regions);

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



void ResolutionComponentHandler::init(){
    AddHYComponent();
    AddCGComponent();
    _simulation.getEnsemble()->setComponentLookUpIDs();
    int n = _simulation.getEnsemble()->getComponents()->size();
    _simulation.getDomain()->getmixcoeff().resize(n*(n-1),1.0);
    for(int i=0;i<n;++i){
    Log::global_log->info()<<"NAME: "<<_simulation.getEnsemble()->getComponents()->at(i).getName()<<" ID :"<<_simulation.getEnsemble()->getComponents()->at(i).getLookUpId()<<std::endl;
    }
    Log::global_log->info()<<"[ResolutionHandler]Added "<<n<<" components"<<std::endl;
    Log::global_log->info()<<"[ResolutionHandler]Initialized"<<std::endl;
}


ResolutionType ResolutionComponentHandler::GetMoleculeResolution(Molecule& m){
    if(m.component()->getName() == "CG")
        return ResolutionType::CoarseGrain;

    if(m.component()->getName() == "HY")
        return ResolutionType::Hybrid;

    if(m.component()->getName()== "FP")
        return ResolutionType::FullParticle;

    //Ideally this is a mistake
    return ResolutionType::ResolutionCount;
}

void ResolutionComponentHandler::AddHYComponent(){

    Component hy{};
    hy.setID(1);
    hy.setName("HY");

    Component& fp = _simulation.getEnsemble()->getComponents()->at(0);
    double total_mass=0;
    total_mass = fp.m();
    std::array<double,3> com{0,0,0};
    for(int lj =0;lj<fp.numLJcenters();++lj){
        LJcenter& site = fp.ljcenter(lj);
        for(int i =0;i<com.size();++i){
            com[i] += site.m()*site.r()[i];
        }
    }

    for(int i=0;i<com.size();++i){
        com[i] /= total_mass;
    }
    LJcenter cg_site;
    cg_site.setR(0,com[0]);
    cg_site.setR(1,com[1]);
    cg_site.setR(2,com[2]);
    cg_site.setM(0);
    std::string cg_name="LJ126";cg_site.setName(cg_name);
    hy.addLJcenter(cg_site);

    //Copy FP LJ sites
    for(auto& lj: fp.ljcenters()){
        LJcenter site = lj;
        hy.addLJcenter(site);
    }

    hy.setI11(fp.I11());
    hy.setI22(fp.I22());
    hy.setI33(fp.I33());

    _simulation.getEnsemble()->getComponents()->push_back(hy);
}

void ResolutionComponentHandler::AddCGComponent(){
    Component cg;
    cg.setID(2);
    cg.setName("CG");
    Component& fp = _simulation.getEnsemble()->getComponents()->at(0);

    LJcenter single_site{};
    single_site.setEps(0);
    single_site.setSigma(1.0);
    std::string site_name{"LJ126"};
    single_site.setName(site_name);

    single_site.setR(0,0.0);
    single_site.setR(1,0.0);
    single_site.setR(2,0.0);

    cg.setI11(fp.I11());
    cg.setI22(fp.I22());
    cg.setI33(fp.I33());

    //sum all mass from the multi-site component
    double cg_mass = 0.0;
    for(int i=0;i<fp.numLJcenters();++i){
        cg_mass += fp.ljcenter(i).m();
    }
    single_site.setM(cg_mass);
    cg.addLJcenter(single_site);

    _simulation.getEnsemble()->getComponents()->push_back(cg);

}

void ResolutionComponentHandler::CheckContainerResolution(ParticleContainer* pc){
    #if defined _OPENMP
    #pragma omp parallel
    #endif
    for(auto it = pc->iterator(ParticleIterator::ALL_CELLS);it.isValid();++it){
        
        bool stop=false;

        for(auto& reg:regions){
            if(reg.IsInsideResolutionRegion(it->r_arr(),FullParticle)){
                //Modify if needed
                CheckMoleculeResolution(*it,FullParticle);
                stop=true;
                break;
            }
        }

        if(stop) continue;

        for(auto& reg:regions){
            if(reg.IsInsideResolutionRegion(it->r_arr(),Hybrid)){
                //Modify if needed
                CheckMoleculeResolution(*it,Hybrid);
                stop=true;
                break;
            }
        }

        if(stop) continue;

        CheckMoleculeResolution(*it, CoarseGrain);
    }
}




void ResolutionComponentHandler::CheckMoleculeResolution(Molecule& m, ResolutionType target){
    ResolutionType current = GetMoleculeResolution(m);
    if(current != target){
        ModifyMoleculeResolution(m,target);
    }
}

void ResolutionComponentHandler::ModifyMoleculeResolution(Molecule& m, ResolutionType target){
    if(target == FullParticle){
        Refine(m);
    }

    if(target == CoarseGrain){
        Coarsen(m);
    }

    if(target == Hybrid){
        Hybridize(m);
    }


}

void ResolutionComponentHandler::Coarsen(Molecule& mol){

    // std::array<double, 3> com{0,0,0};
    // double total_mass = mol.component()->m();
    // for(int lj=0;lj<mol.component()->numLJcenters();++lj){
    //     auto lj_site = mol.ljcenter_d_abs(lj);
    //     double site_mass = mol.component()->ljcenter(lj).m();
    //     for(int i=0;i<3;++i){
    //         com[i] += lj_site[i]*site_mass;
    //         com[i] = com[i]/total_mass;
    //     }
    // }
    Component* component = _simulation.getEnsemble()->getComponent(2);
    mol.setComponent(component);
    //TODO:Set position to COM....how?
}

void ResolutionComponentHandler::Refine(Molecule& mol){
    Component* component = _simulation.getEnsemble()->getComponent(0);
    mol.setComponent(component);
}

void ResolutionComponentHandler::Hybridize(Molecule& mol){
    Component* component = _simulation.getEnsemble()->getComponent(1);
    //TODO: also here, is this correct? Or a more elaborate thing is needed?
    mol.setComponent(component);
}