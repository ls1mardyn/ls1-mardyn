#include "ProfilerPMF.h"

InternalCellProcessor::InternalCellProcessor(const double cr, InternalProfiler* r):CellProcessor{cr,cr},my_profiler{r}{

}

double InternalCellProcessor::DistanceBetweenCOMs(std::array<double,3>& c1, std::array<double,3>& c2){
    double r =0.0;
    std::array<double,3> diff={0.0,0.0,0.0};
    
    for(int i=0;i<diff.size();i++){
        diff[i]=c1[i]-c2[i];
    }

    r = diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2];
    return r;

}

double InternalCellProcessor::PotentialCallBack(double e, double s, double r2){
    double pot=0;

    double inv_r2 = 1.0/r2;
    double lj6 = s*s*inv_r2;
    lj6 = lj6*lj6*lj6;
    double lj12 = lj6*lj6;
    double lj12m6 = lj12 - lj6;
    pot = lj12m6*e*4.0;
    if(std::isinf(pot)) return 0.0;
    return pot;
}

void InternalCellProcessor::processCell(ParticleCell& cell){
    if(cell.isInnerCell() || cell.isBoundaryCell()){
        auto begin = cell.iterator();
        double distance=0.0;
        for(auto it1 = begin;it1.isValid();++it1){
            std::array<double,3> com1={0.0,0.0,0.0};
            Molecule& m1 = *it1;
            double eps, sig, pot;
            sig = m1.component()->getSigma(0);
            eps = m1.component()->getEps(0);
            com1 = my_profiler->GetCOM(&m1);
            auto it2 = it1;
            ++it2;
            for(;it2.isValid();++it2){
                Molecule& m2 = *it2;
                std::array<double,3> com2={0.0,0.0,0.0};
                com2 = my_profiler->GetCOM(&m2);
                mardyn_assert(&m1 != &m2);
                //Now we compute the distance between the COMs
                distance = DistanceBetweenCOMs(com1,com2);
                if(distance < _cutoffRadiusSquare){

                    pot = PotentialCallBack(eps,sig,distance);
                    my_profiler->ProcessDistance(distance, pot);
                }
    
            }
    
        }
    }
}

void InternalCellProcessor::processCellPair(ParticleCell& c1, ParticleCell& c2, bool sumAll){

    auto begin1 = c1.iterator();
    auto begin2 = c2.iterator();
    double distance=0.0;
    std::array<double,3> com1={0.0,0.0,0.0};
    std::array<double,3> com2={0.0,0.0,0.0};
    if(sumAll){
        for(auto it1=begin1;it1.isValid();++it1){
            Molecule& m1 = *it1;
            double eps, sig, pot;
            sig = m1.component()->getSigma(0);
            eps = m1.component()->getEps(0);

            com1 = my_profiler->GetCOM(&m1);
            for(auto it2 =begin2;it2.isValid();++it2){
                Molecule& m2 = *it2;
                com2 = my_profiler->GetCOM(&m2);
                distance = DistanceBetweenCOMs(com1,com2);
                if(distance < _cutoffRadiusSquare){
                    pot = PotentialCallBack(eps,sig,distance);
                    my_profiler->ProcessDistance(distance,pot);
                }
            }
        }
    }
    else{
        if(c1.isInnerCell()){//no hallo cells at all

            for(auto it1=begin1;it1.isValid();++it1){
                Molecule& m1 = *it1;
                double eps, sig, pot;
                sig = m1.component()->getSigma(0);
                eps = m1.component()->getEps(0);
                com1 = my_profiler->GetCOM(&m1);
                for(auto it2=begin2;it2.isValid();++it2){
                    Molecule& m2 = *it2;
                    com2 = my_profiler->GetCOM(&m2);
                    distance = DistanceBetweenCOMs(com1,com2);
                    if(distance < _cutoffRadiusSquare){
                        pot = PotentialCallBack(eps,sig,distance);
                        my_profiler->ProcessDistance(distance,pot);
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
                double eps, sig, pot;
                sig = m1.component()->getSigma(0);
                eps = m1.component()->getEps(0);
                com1 = my_profiler->GetCOM(&m1);
                for(auto it2=begin2;it2.isValid();++it2){
                    Molecule& m2 = *it2;
                    com2 = my_profiler->GetCOM(&m2);
                    distance = DistanceBetweenCOMs(com1,com2);
                    
                    if(distance < _cutoffRadiusSquare){
                        pot = PotentialCallBack(eps,sig,distance);
                        my_profiler->ProcessDistance(distance,pot);
                    }
                }
            }
        }
    }

}