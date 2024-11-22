#include "ProfilerPMF.h"

InternalCellProcessor::InternalCellProcessor(const double cr, int bins, double width):CellProcessor{cr,cr},bin_width{width}{

    thread_data.resize(mardyn_get_max_threads());
    global_buffer.resize(bins,0.0);
    #pragma omp parallel
    {
        thread_data[mardyn_get_thread_num()].resize(bins,0.0);
    }

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

void InternalCellProcessor::processCell(ParticleCell& cell){
    if(cell.isInnerCell() || cell.isBoundaryCell()){
        auto begin = cell.iterator();
        double distance=0.0;
        for(auto it1 = begin;it1.isValid();++it1){
            Molecule& m1 = *it1;
            auto it2 = it1;
            ++it2;
            for(;it2.isValid();++it2){
                Molecule& m2 = *it2;
                mardyn_assert(&m1 != &m2);
                this->ProcessPairData(m1,m2);
            }
    
        }
    }
}

void InternalCellProcessor::processCellPair(ParticleCell& c1, ParticleCell& c2, bool sumAll){

    auto begin1 = c1.iterator();
    auto begin2 = c2.iterator();
    if(sumAll){
        for(auto it1=begin1;it1.isValid();++it1){
            Molecule& m1 = *it1;
            for(auto it2 =begin2;it2.isValid();++it2){
                Molecule& m2 = *it2;
                this->ProcessPairData(m1,m2);
            }
        }
    }
    else{
        if(c1.isInnerCell()){//no hallo cells at all

            for(auto it1=begin1;it1.isValid();++it1){
                Molecule& m1 = *it1;
                for(auto it2=begin2;it2.isValid();++it2){
                    Molecule& m2 = *it2;
                    this->ProcessPairData(m1,m2);
                }
            }

        }

        if(c1.isBoundaryCell()){//c1 is  boundary
            if(c2.isHaloCell() && !(c1.getCellIndex()<c2.getCellIndex())){
                return;
            }

            for(auto it1=begin1;it1.isValid();++it1){
                Molecule& m1 = *it1;
                for(auto it2=begin2;it2.isValid();++it2){
                    Molecule& m2 = *it2;
                    this->ProcessPairData(m1,m2);
                }
            }
        }
    }

}

void InternalCellProcessor::endTraversal(){
    
    for(int b=0;b<global_buffer.size();++b){
        for(int t=0;t<thread_data.size();++t){
            global_buffer[b] += thread_data[t][b];
        }
    }
}

void InternalCellProcessor::ProcessPairData(Molecule& m1, Molecule& m2){
    double distance2=0.0;
    std::array<double,3> com1={0.0,0.0,0.0};
    std::array<double,3> com2={0.0,0.0,0.0};
    com1 = ComputeCOM(m1);
    com2 = ComputeCOM(m2);
    distance2 = DistanceBetweenCOMs(com1,com2);

    if(distance2 < _cutoffRadiusSquare){
        int index = std::floor(std::sqrt(distance2)/bin_width);
        thread_data[mardyn_get_thread_num()][index]++;
    }

}