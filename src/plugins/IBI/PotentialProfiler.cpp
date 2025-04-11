#include "PotentialProfiler.h"
#include "IBI_Math.h"

void PotentialProfiler::readXML(XMLfileUnits& file) {
    file.getNodeValue("totalBins",number_bins);
    Log::global_log->info()<<"[Potential Profiler] Total bins "<<number_bins<<"\n";
    file.getNodeValue("sampleFrequency",sample_frequency);
    Log::global_log->info()<<"[Potential Profiler] Sample frequency "<<sample_frequency<<"\n";
    file.getNodeValue("avgFrequency", avg_frequency);
    if (avg_frequency % sample_frequency != 0) {
        Log::global_log->info() << "[Potential Profiler] Average frequency is not multiple of sample frequency."
                                << " Will use the next larger valid Average Frequency." << std::endl;
        avg_frequency = (static_cast<int>(avg_frequency / sample_frequency) + 1) * sample_frequency;
    }
    Log::global_log->info() << "[Potential Profiler] Average frequency " << avg_frequency << "\n";

}

void PotentialProfiler::init(ParticleContainer* pc, DomainDecompBase* dd, Domain* dom) {
    bin_width = pc->getCutoff() / number_bins;
    const auto max_r = bin_width * number_bins;
    cell_processor = new InternalCellProcessor(global_simulation->getcutoffRadius(), bin_width, number_bins, max_r);

    current_potential_average.resize(number_bins, 0.0);
    measured_steps = 0;
}

void PotentialProfiler::afterForces(ParticleContainer *pc, DomainDecompBase *domainDecomp, unsigned long simstep) {
    if (simstep % sample_frequency == 0 && simstep > global_simulation->getInitStatistics()) {
        pc->traverseCells(*cell_processor);

        auto& pots = cell_processor->get_sample_pots();
        auto& counts = cell_processor->get_sample_counts();
        for (int idx = 0; idx < pots.size(); idx++) current_potential_average[idx] += pots[idx] / counts[idx];

        measured_steps++;
    }
}

void PotentialProfiler::endStep(ParticleContainer* pc, DomainDecompBase* dd, Domain* dom, unsigned long simstep){
    if (simstep % sample_frequency == 0 && simstep > global_simulation->getInitStatistics()) {
        WriteDataToFile(simstep);
    }
}

void PotentialProfiler::WriteDataToFile(unsigned long simstep) {
    std::string filename = "potential_data_" + std::to_string(simstep) + ".txt";
    std::ofstream outfile(filename);

    auto& pots = cell_processor->get_sample_pots();
    auto& counts = cell_processor->get_sample_counts();
    for (int i = 0; i < number_bins; i++) {
        const double rmid = (i+0.5)*bin_width;

        outfile << std::setw(8) << std::left << rmid <<"\t"
                << std::setw(8) << std::left << pots[i] / counts[i];
        if (simstep % avg_frequency == 0) {
            outfile << "\t" << std::setw(8) << std::left << current_potential_average[i] / measured_steps;
        }
        outfile << std::endl;
    }

    if (simstep % avg_frequency == 0) {
        measured_steps = 0;
        std::fill(current_potential_average.begin(), current_potential_average.end(), 0.0);
    }

    outfile.close();
}

//=======================================================================
// InternalCellProcessor
//=======================================================================

PotentialProfiler::InternalCellProcessor::InternalCellProcessor(double cr, double bin_width, int num_bins, double max_r) :
        CellProcessor{cr, cr}, single_sample_counts(num_bins, 0.0), single_sample_pots(num_bins, 0.0), bin_width(bin_width), max_r(max_r) {
    const int numThreads = mardyn_get_max_threads();
    thread_single_sample_counts.resize(numThreads);
    thread_single_sample_pots.resize(numThreads);
    #if defined(_OPENMP)
    #pragma omp parallel
    #endif
    {
        const int myid = mardyn_get_thread_num();
        thread_single_sample_counts[myid].resize(num_bins, 0.0);
        thread_single_sample_pots[myid].resize(num_bins, 0.0);
    } // end pragma omp parallel
}

void PotentialProfiler::InternalCellProcessor::initTraversal() {
    #if defined(_OPENMP)
    #pragma omp parallel
    #endif
    {
        const int myid = mardyn_get_thread_num();
        std::fill(thread_single_sample_counts[myid].begin(), thread_single_sample_counts[myid].end(), 0.0);
        std::fill(thread_single_sample_pots[myid].begin(), thread_single_sample_pots[myid].end(), 0.0);
    } // end pragma omp parallel
}

void PotentialProfiler::InternalCellProcessor::endTraversal() {
    for (int bin = 0; bin < single_sample_counts.size(); bin++) {
        single_sample_counts[bin] = 0;
        single_sample_pots[bin] = 0;
        for (int thread = 0; thread < thread_single_sample_counts.size(); thread++) {
            single_sample_counts[bin] += thread_single_sample_counts[thread][bin];
            single_sample_pots[bin] += thread_single_sample_pots[thread][bin];
        }
    }
}

void PotentialProfiler::InternalCellProcessor::processCell(ParticleCell& cell){
    if(cell.isInnerCell() || cell.isBoundaryCell()) {
        const int thread = mardyn_get_thread_num();
        for(auto it1 = cell.iterator();it1.isValid();++it1){
            Molecule& m1 = *it1;

            auto it2 = it1; ++it2;
            for (; it2.isValid(); ++it2) {
                Molecule& m2 = *it2;
                mardyn_assert(&m1 != &m2);
                handleMoleculePair(m1, m2, thread);
            }
        }
    }
}

void PotentialProfiler::InternalCellProcessor::processCellPair(ParticleCell& c1, ParticleCell& c2, bool sumAll){
    const int thread = mardyn_get_thread_num();

    if(sumAll){
        for (auto it1 = c1.iterator(); it1.isValid(); ++it1) {
            Molecule& m1 = *it1;
            for (auto it2 = c2.iterator(); it2.isValid(); ++it2) {
                Molecule& m2 = *it2;
                handleMoleculePair(m1, m2, thread);
            }
        }
    }

    else {
        if (c1.isInnerCell()) { //no hallo cells at all
            for(auto it1 = c1.iterator(); it1.isValid(); ++it1) {
                Molecule& m1 = *it1;
                for(auto it2 = c2.iterator(); it2.isValid(); ++it2) {
                    Molecule& m2 = *it2;
                    handleMoleculePair(m1, m2, thread);
                }
            }

        }

        if(c1.isBoundaryCell()){//c1 is  boundary
            if(c2.isHaloCell() && c1.getCellIndex() >= c2.getCellIndex()){
                return;
            }

            for (auto it1 = c1.iterator(); it1.isValid(); ++it1) {
                Molecule& m1 = *it1;
                for (auto it2 = c2.iterator(); it2.isValid(); ++it2) {
                    Molecule& m2 = *it2;
                    handleMoleculePair(m1, m2, thread);
                }
            }
        }
    }

}

void PotentialProfiler::InternalCellProcessor::handleMoleculePair(Molecule &m1, Molecule &m2, int thread) {
    const std::array<double,3> com1 = GetCOM(m1);
    const std::array<double,3> com2 = GetCOM(m2);

    //Now we compute the distance between the COMs
    const double r2 = Distance2(com1,com2);
    if(r2 < _cutoffRadiusSquare) {
        const double r = std::sqrt(r2);
        if(r > max_r) return;

        int index = std::floor(r/bin_width);
        const double sigma = m1.component()->getSigma(0); // TODO: this might need some rethinking
        const double epsilon = m1.component()->getEps(0); // TODO: this might need some rethinking
        thread_single_sample_counts[thread][index] += 1.0;
        thread_single_sample_pots[thread][index] += calcPot(sigma, epsilon, r2);
    }
}

double PotentialProfiler::InternalCellProcessor::calcPot(double sigma, double epsilon, double r2) {
    const double inv_r2 = 1.0 / r2;
    const double lj6 = std::pow(sigma * sigma * inv_r2, 3.0);
    const double lj12 = lj6 * lj6;
    const double lj12m6 = lj12 - lj6;
    const double pot = lj12m6 * epsilon * 4.0;

    if (std::isinf(pot)) return 0.0;
    return pot;
}
