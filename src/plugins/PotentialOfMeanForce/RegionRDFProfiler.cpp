#include "Statistics.h"

RegionRDFProfiler::RegionRDFProfiler(ResRegion& r):region{r}{

}


void RegionRDFProfiler::init(int rc){

    centers.resize(total_bins);
    bin_width = rc/(double)total_bins;
    InitCenters();

    Log::global_log->info()<<"[AdResSStatistics] Measuring local RDFs with "<<total_bins<<" bins"<<std::endl;
    Log::global_log->info()<<"[AdResSStatistics] Bin width"<<bin_width<<std::endl;
    processor = new RegionCellProcessor(rc,total_bins,bin_width,region);
    processor->SetComponentName(region.component_name.data());
}

void RegionRDFProfiler::MeasureRDF(ParticleContainer* particleContainer, unsigned long step){
    if(step%sample_frequency ==0){
        measured_steps++;
        particleContainer->traverseCells(*processor);
    }
}

int RegionRDFProfiler::CountMoleculesInRegion(ParticleContainer* pc){

    std::array<double,3> low = region.low;
    std::array<double,3> high = region.high;

    //correct for halo cells
    //low[1] = low[1]-_simulation.getcutoffRadius();
    //low[2] = low[2]-_simulation.getcutoffRadius();
    //high[1] = high[1]+_simulation.getcutoffRadius();
    //high[2] = high[2]+_simulation.getcutoffRadius();

    int N=0;
    #pragma omp parallel reduction(+:N)
    for(auto it = pc->regionIterator(low.data(),high.data(),ParticleIterator::ONLY_INNER_AND_BOUNDARY);it.isValid();++it) {
        N++;
    }   

    return N;
}

void RegionRDFProfiler::PrintOutput2Files(unsigned long simstep,ParticleContainer* pc){
    if(simstep>0 && simstep%output_frequency == 0){
        std::string file_name=file_prefix+std::to_string(simstep)+".txt";
        std::ofstream output_file(file_name);

        std::vector<double>& data = processor->GetPairCountBuffer();
        int N = CountMoleculesInRegion(pc);
        double region_volume = RegionVolume();
        for(int i=0;i<total_bins;++i){
            double rmin, rmax, binvol, rmin3,rmax3;
            rmin = i*bin_width;
            rmax =(i+1)*bin_width;
            rmin3 = rmin*rmin*rmin;
            rmax3 = rmax*rmax*rmax;
            binvol = (4.0/3.0)*M_PI*(rmax3-rmin3);

            double val = data[i]*region_volume/binvol;
            val /= (0.5*N*(N-1.0));
            val /= measured_steps;

            output_file<<std::setw(8)<<std::left<<centers[i]<<"\t"
            <<std::setw(8)<<std::left<<val<<"\t"
            <<std::setw(8)<<std::left<<data[i]/measured_steps<<"\t"
            <<std::setw(8)<<std::left<<binvol<<"\t"
            <<std::setw(8)<<std::left<<region_volume<<"\t"
            <<std::setw(8)<<std::left<<N<<"\t"
            <<std::endl;
        }

        output_file.close();
        
    }
}


/***************
 * 
 * 
 * CELL PROCESSOR FUNCTIONS
 * bunch of code repetition
 * 
 */

RegionCellProcessor::RegionCellProcessor(const double rc, int bins, double width,ResRegion& reg):CellProcessor{rc,rc},bin_width{width},region{reg}{

    thread_data.resize(mardyn_get_max_threads());
    global_buffer.resize(bins,0.0);
    #pragma omp parallel
    {
        thread_data[mardyn_get_thread_num()].resize(bins,0.0);
    }

}


void RegionCellProcessor::processCell(ParticleCell& cell){
    if(!CellInRegion(cell)){
        return;
    }
    if(cell.isInnerCell() || cell.isBoundaryCell()){
        auto begin = cell.iterator();
        double distance=0.0;
        for(auto it1 = begin;it1.isValid();++it1){
            Molecule& m1 = *it1;
            if(m1.component()->getName() != component_name)
            continue;
            auto it2 = it1;
            ++it2;
            for(;it2.isValid();++it2){
                Molecule& m2 = *it2;
                if(m2.component()->getName() != component_name)
                    continue;
                mardyn_assert(&m1 != &m2);
                this->ProcessPairData(m1,m2);
            }
    
        }
    }
}

void RegionCellProcessor::processCellPair(ParticleCell& c1, ParticleCell& c2, bool sumAll){
    if(!CellInRegion(c1) || !CellInRegion(c2)){
        return;
    }

    auto begin1 = c1.iterator();
    auto begin2 = c2.iterator();
    if(sumAll){
        for(auto it1=begin1;it1.isValid();++it1){
            Molecule& m1 = *it1;
            //if(m1.component()->getName() != component_name)
            //continue;
            for(auto it2 =begin2;it2.isValid();++it2){
                Molecule& m2 = *it2;
                //if(m2.component()->getName() != component_name)
                //continue;
                this->ProcessPairData(m1,m2);
            }
        }
    }
    else{
        if(c1.isInnerCell()){//no hallo cells at all

            for(auto it1=begin1;it1.isValid();++it1){
                Molecule& m1 = *it1;
                //if(m1.component()->getName() != component_name)
                //continue;
                for(auto it2=begin2;it2.isValid();++it2){
                    Molecule& m2 = *it2;
                    //if(m2.component()->getName() != component_name)
                    //continue;
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
                //if(m1.component()->getName() != component_name)
                //continue;
                for(auto it2=begin2;it2.isValid();++it2){
                    Molecule& m2 = *it2;
                    //if(m2.component()->getName() != component_name)
                    //continue;
                    this->ProcessPairData(m1,m2);
                }
            }
        }
    }

}

void RegionCellProcessor::ProcessPairData(Molecule& m1, Molecule& m2){
    double distance2=0.0;
    std::array<double,3> com1={0.0,0.0,0.0};
    std::array<double,3> com2={0.0,0.0,0.0};
    com1 = ComputeCOM(m1);
    com2 = ComputeCOM(m2);
    distance2 = Distance2BetweenCOMs(com1,com2);

    if(distance2 < _cutoffRadiusSquare){
        int index = std::floor(std::sqrt(distance2)/bin_width);
        thread_data[mardyn_get_thread_num()][index]++;
    }

}


void RegionCellProcessor::endTraversal(){
    
    for(int b=0;b<global_buffer.size();++b){
        for(int t=0;t<thread_data.size();++t){
            global_buffer[b] += thread_data[t][b];
        }
    }
}


bool RegionCellProcessor::CellInRegion(ParticleCell& cell){
    const double cell_min = cell.getBoxMin(0);
    const double cell_max = cell.getBoxMax(0);
    const double reg_min = region.low[0];
    const double reg_max = region.high[0];

    if(cell_max < reg_max && cell_min > reg_min) return true; //fully inside
    if(cell_max > reg_min && cell_max < reg_max) return true; //Cell upper is within region
    if(cell_min > reg_min && cell_min < reg_max) return true; //Cell lower is within region

    return false;
}