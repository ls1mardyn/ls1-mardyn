#include"RDFAtCOM.h"

void RadialDFCOM::readXML(XMLfileUnits& file){

}

void RadialDFCOM::init(ParticleContainer* pc, DomainDecompBase* dd, Domain* dom){

}

void RadialDFCOM::endStep(ParticleContainer* pc, DomainDecompBase* dd, Domain* dom, unsigned long simstep){
    //print all the molecule positions
    for(auto it = pc->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);it.isValid();++it){
        
        std::array<double,3> com{0.0,0.0,0.0};
        std::array<double,3> com2{0.0,0.0,0.0};
        com2 = GetCOM(&(*it));
        double total_mass;//TODO:Can be computed once for every component
        Component* comp = it->component();
        total_mass = comp->m();
        for(int lj=0;lj<comp->numLJcenters();lj++){
            LJcenter& lj_center = comp->ljcenter(lj);
            for(int idx=0;idx<lj_center.r().size();idx++){               
                com[idx] += lj_center.m()*it->ljcenter_d_abs(lj)[idx];
            }
  
        }

        for(int qs=0;qs<comp->numCharges();qs++){
            Charge& q_center = comp->charge(qs);
            for(int idx=0;idx<q_center.r().size();idx++){               
                com[idx] += q_center.m()*it->charge_d_abs(qs)[idx];
            }
        }
        std::cout<<"Molecule "<<it->getID()<<" r: (";
        for(int i=0;i<it->r_arr().size();i++){
            std::cout<<it->r(i)<<" ";
        }
        std::cout<<")\n";
        std::cout<<"Molecule "<<it->getID()<<" center of mass (";
        for(int i=0;i<com.size();i++){
            std::cout<<com[i]/total_mass<<",";
        }
        std::cout<<")\n(";
        for(int i=0;i<com2.size();i++){
            std::cout<<com2[i]<<" ";
        }
        std::cout<<")\n";
        for(int i=0;i<it->numSites();i++){
            std::cout<<"(";
            for(int j=0;j<3;j++){
                std::cout<<it->site_d_abs(i)[j]<<",";
            }
            std::cout<<")\n";
        }


    }
    
}

std::array<double,3> RadialDFCOM::GetCOM(Molecule* m){

    std::array<double,3> com{0.0,0.0,0.0};
    double total_mass;//TODO:Can be computed once for every component
    Component* comp = m->component();
    total_mass = comp->m();
    for(int lj=0;lj<comp->numLJcenters();lj++){
        LJcenter& lj_center = comp->ljcenter(lj);
        for(int idx=0;idx<lj_center.r().size();idx++){               
            com[idx] += lj_center.m()*m->ljcenter_d_abs(lj)[idx];
        }

    }
    for(int qs=0;qs<comp->numCharges();qs++){
        Charge& q_center = comp->charge(qs);
        for(int idx=0;idx<q_center.r().size();idx++){               
            com[idx] += q_center.m()*m->charge_d_abs(qs)[idx];
        }
    }

    for(int i=0;i<com.size();i++){
        com[i] /= total_mass;
    }

    return com;

}

/************************
 * **********************
 * Cell processor methods
 * **********************
 * **********************
 ***********************/

double COMDistanceCellProcessor::DistanceBetweenCOMs(std::array<double,3>& c1, std::array<double,3>& c2){
    double r =0.0;
    std::array<double,3> diff={0.0,0.0,0.0};
    
    for(int i=0;i<diff.size();i++){
        diff[i]=c1[i]-c2[i];
    }

    r = std::sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);

    return r;

}

void COMDistanceCellProcessor::processCell(ParticleCell& cell){
    auto begin = cell.iterator();
    
    for(auto it1 = begin;it1.isValid();++it1){
        std::array<double,3> com1={0.0,0.0,0.0};
        Molecule& m1 = *it1;
        com1 = rdf->GetCOM(&m1);
        auto it2 = it1;
        ++it2;
        for(;it2.isValid();++it2){
            Molecule& m2 = *it2;
            std::array<double,3> com2={0.0,0.0,0.0};
            com2 = rdf->GetCOM(&m2);
            mardyn_assert(&m1 != &m2);

            //Now we compute the distance between the COMs
            std::cout<<"COMS m1 (";
            for(int i=0;i<com1.size();i++){
                std::cout<<com1[i]<<" ";
            }
            std::cout<<")\n";
            std::cout<<"COMS m2 (";
            for(int i=0;i<com2.size();i++){
                std::cout<<com2[i]<<" ";
            }
            std::cout<<")\n";
            std::cout<<"Molecule m1 "<<m1.getID()<<" distance to m2 "<<m2.getID()<<" is "<<DistanceBetweenCOMs(com1,com2)<<"\n";
        }

    }
}

void COMDistanceCellProcessor::processCellPair(ParticleCell& c1, ParticleCell& c2, bool sumAll){

}