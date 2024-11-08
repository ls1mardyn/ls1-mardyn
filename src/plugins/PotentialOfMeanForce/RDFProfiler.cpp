#include "RDFProfiler.h"
#include "IBI_Math.h"

RDFProfiler::RDFProfiler(): cell_processor{nullptr}, measured_steps{0} { }

void RDFProfiler::init(int bins) {
    number_bins = bins;
    bin_width = _simulation.getcutoffRadius() / number_bins;
    max_r = bin_width * number_bins;
    cell_processor = new InternalCellProcessor(global_simulation->getcutoffRadius(), bin_width, number_bins, max_r);

    r_nodes.resize(number_bins);
    for (int i = 0; i < number_bins; ++i) r_nodes[i] = (i + 0.5) * bin_width;

    rdf_buffer.resize(number_bins);
    std::fill(rdf_buffer.begin(),rdf_buffer.end(),0.0);;

    Log::global_log->info() << "[PMF] Internal profiler bin width " << bin_width << "\n";
    Log::global_log->info() << "[PMF] Limit distance " << max_r << "\n";
    Log::global_log->info() << "[PMF] Profiler enabled " << std::endl;
}

void RDFProfiler::ProfileData(ParticleContainer* pc) {
    pc->traverseCells(*cell_processor);

    auto& rdf_sample = GetRDFPrevStepCounts();
    for (int bin = 0; bin < number_bins; bin++) {
        rdf_buffer[bin] += rdf_sample[bin];
    }

    measured_steps++;
}

void RDFProfiler::GetRDFTotal(std::vector<double> &buffer) {
    buffer.resize(number_bins);
    std::copy(rdf_buffer.begin(), rdf_buffer.end(), buffer.begin());
    normalize(buffer, measured_steps);
}

void RDFProfiler::GetRDFPrevStep(std::vector<double> &buffer) {
    buffer.resize(number_bins);
    std::copy(GetRDFPrevStepCounts().begin(), GetRDFPrevStepCounts().end(), buffer.begin());
    normalize(buffer, 1);
}

void RDFProfiler::normalize(std::vector<double> &buffer, int steps) {
    auto* domain = _simulation.getDomain();
    const auto N = domain->getglobalNumMolecules();
    const auto V = domain->getGlobalVolume();

    for (int bin = 0; bin < number_bins; ++bin) {
        //Generate RDF g(r) data
        const auto  rmin = bin * bin_width;
        const auto  rmax =(bin + 1) * bin_width;
        const auto  rmin3 = rmin * rmin * rmin;
        const auto  rmax3 = rmax * rmax * rmax;
        const auto  binvol = (4.0 / 3.0) * M_PI * (rmax3 - rmin3);
        const auto  den = 0.5 * N * (N - 1.0) * binvol / V;
        buffer[bin] /= steps;
        buffer[bin] /= den;
    }
}

//=======================================================================
// InternalCellProcessor
//=======================================================================

RDFProfiler::InternalCellProcessor::InternalCellProcessor(double cr, double bin_width, int num_bins, double max_r) :
CellProcessor{cr, cr}, single_sample_counts(num_bins, 0.0), bin_width(bin_width), max_r(max_r) {
    const int numThreads = mardyn_get_max_threads();
    thread_single_sample_counts.resize(numThreads);
    #if defined(_OPENMP)
    #pragma omp parallel
    #endif
    {
        const int myid = mardyn_get_thread_num();
        thread_single_sample_counts[myid].resize(num_bins, 0.0);
    } // end pragma omp parallel
}

void RDFProfiler::InternalCellProcessor::initTraversal() {
    #if defined(_OPENMP)
    #pragma omp parallel
    #endif
    {
        const int myid = mardyn_get_thread_num();
        std::fill(thread_single_sample_counts[myid].begin(), thread_single_sample_counts[myid].end(), 0.0);
    } // end pragma omp parallel
}

void RDFProfiler::InternalCellProcessor::endTraversal() {
    for (int bin = 0; bin < single_sample_counts.size(); bin++) {
        single_sample_counts[bin] = 0;
        for (int thread = 0; thread < thread_single_sample_counts.size(); thread++) {
            single_sample_counts[bin] += thread_single_sample_counts[thread][bin];
        }
    }
}

void RDFProfiler::InternalCellProcessor::processCell(ParticleCell& cell){
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

void RDFProfiler::InternalCellProcessor::processCellPair(ParticleCell& c1, ParticleCell& c2, bool sumAll){
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

void RDFProfiler::InternalCellProcessor::handleMoleculePair(Molecule &m1, Molecule &m2, int thread) {
    const std::array<double,3> com1 = GetCOM(m1);
    const std::array<double,3> com2 = GetCOM(m2);

    //Now we compute the distance between the COMs
    const double r2 = Distance2(com1,com2);
    if(r2 < _cutoffRadiusSquare) {
        const double r = std::sqrt(r2);
        if(r > max_r) return;
        int index = std::floor(r/bin_width);
        thread_single_sample_counts[thread][index]++;
    }
}
