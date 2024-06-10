#include"RDFAtCOM.h"

RadialDFCOM::RadialDFCOM():cell_processor(new COMDistanceCellProcessor(this)),measured_steps{0}{

}

void RadialDFCOM::readXML(XMLfileUnits& file){
    file.getNodeValue("totalBins",number_bins);
    Log::global_log->info()<<"[RDF COM] Total bins "<<number_bins<<"\n";
    file.getNodeValue("sampleFrequency",sample_frequency);
    Log::global_log->info()<<"[RDF COM] Sample frequency "<<sample_frequency<<"\n";
}

void RadialDFCOM::init(ParticleContainer* pc, DomainDecompBase* dd, Domain* dom){
    SetBinContainer(pc);
}

void RadialDFCOM::endStep(ParticleContainer* pc, DomainDecompBase* dd, Domain* dom, unsigned long simstep){
    if(simstep%sample_frequency ==0){
        sample_frequency++;
        pc->traverseCells(*cell_processor);    
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

void RadialDFCOM::SetBinContainer(ParticleContainer* pc){
    bin_width = pc->getCutoff()/(double)number_bins;
    this->bin_counts.resize(number_bins);
    std::fill(bin_counts.begin(),bin_counts.end(),0.0);
    Log::global_log->info()<<"[RDF COM] Bin  width "<<bin_width<<"\n";

}

void RadialDFCOM::ProcessDistance(double r){  

    int index = std::floor(r/bin_width);
    this->bin_counts[index]++;

}   

void RadialDFCOM::WriteRDFToFile(){
    std::ofstream outfile("rdf.txt");
    outfile<<"#Total time steps: "<<_simulation.getNumTimesteps()<<"\n";
    double data=0.0;
    for(int i=0;i<number_bins;i++){
        double rmin, rmax, rmid, binvol, rmin3,rmax3;
        rmid = (i+0.5)*bin_width;
        rmin = i*bin_width;
        rmax =(i+1)*bin_width;
        rmin3 = rmin*rmin*rmin;
        rmax3 = rmax*rmax*rmax;
        binvol = (4.0/3.0)*M_PI*(rmax3-rmin3);
        data = (double)bin_counts[i]/(binvol*(double)_simulation.getNumTimesteps());
        outfile<<i<<"\t"<<bin_counts[i]<<"\t"<<rmid<<"\t"<<data<<"\t"<<binvol<<"\n";

    }

    outfile.close();
}

/************************
 * **********************
 * Cell processor methods
 * **********************
 * **********************
 ***********************/

COMDistanceCellProcessor::COMDistanceCellProcessor(RadialDFCOM* r):CellProcessor(0,0),rdf{r}{

}

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
    double distance=0.0;
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
            distance = DistanceBetweenCOMs(com1,com2);
            if(distance <= _simulation.getcutoffRadius()){
                rdf->ProcessDistance(distance);
            }

        }

    }
}

void COMDistanceCellProcessor::processCellPair(ParticleCell& c1, ParticleCell& c2, bool sumAll){

    auto begin1 = c1.iterator();
    auto begin2 = c2.iterator();
    double distance=0.0;
    std::array<double,3> com1={0.0,0.0,0.0};
    std::array<double,3> com2={0.0,0.0,0.0};
    if(sumAll){

    }
    else{
        if(c1.isInnerCell()){//no hallo cells at all

            for(auto it1=begin1;it1.isValid();++it1){
                Molecule& m1 = *it1;
                com1 = rdf->GetCOM(&m1);
                for(auto it2=begin2;it2.isValid();++it2){
                    Molecule& m2 = *it2;
                    com2 = rdf->GetCOM(&m2);

                    distance = DistanceBetweenCOMs(com1,com2);
                    if(distance <= _simulation.getcutoffRadius()){
                        rdf->ProcessDistance(distance);
                    }

                }
            }

        }

        if(c1.isBoundaryCell()){//c1 is  boundary
            if(c2.isHaloCell() && !(c1.getCellIndex()<c2.getCellIndex())){
                return;
            }

            for(auto it1=begin1;it1.isValid();++it1){
                Molecule& m1 = *it1;
                com1 = rdf->GetCOM(&m1);
                for(auto it2=begin2;it2.isValid();++it2){
                    Molecule& m2 = *it2;
                    com2 = rdf->GetCOM(&m2);
                    distance = DistanceBetweenCOMs(com1,com2);
                    if(distance <= _simulation.getcutoffRadius()){
                        rdf->ProcessDistance(distance);
                    }
                }
            }
        }
    }

}