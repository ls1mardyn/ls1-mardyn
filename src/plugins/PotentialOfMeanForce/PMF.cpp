#include"PMF.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"


void PMF::readXML(XMLfileUnits& xmlfile) {
    xmlfile.getNodeValue("singleRun", mode_single_run);
    // using this double tmp buffer here to allow for easier configuration in xml file
    double tmp = 1e+6;
    xmlfile.getNodeValue("stepsEquil", tmp);
    steps_equilibration = static_cast<int>(tmp);
    tmp = 1e+5;
    xmlfile.getNodeValue("stepsMeasure", tmp);
    steps_measurement = static_cast<int>(tmp);

    xmlfile.getNodeValue("alpha", alpha);

    int bins;
    xmlfile.getNodeValue("bins", bins);
    profiler.init(bins);

    xmlfile.getNodeValue("equilibrate", mode_equilibrate);

    xmlfile.getNodeValue("rdfPath", rdf_path);
    xmlfile.getNodeValue("potPath", pot_path);

    if (rdf_path.empty()) mode_initial_rdf = true;
    if (pot_path.empty()) ibi_iteration = 0;
    else {
        auto idx = pot_path.find_last_of('_');
        if (idx == std::string::npos) throw std::runtime_error("pot file has no iteration suffix.");
        const auto begin = idx + 1;

        idx = pot_path.find_last_of('.');
        if (idx == std::string::npos) throw std::runtime_error("pot file has no file ending.");
        const auto end = idx;
        const auto it_str = pot_path.substr(begin, end - begin);
        ibi_iteration = std::stoi(it_str);
    }

    if (mode_initial_rdf) mode_equilibrate = false;
}

void PMF::init(ParticleContainer* pc, DomainDecompBase* domainDecomp, Domain* domain) {
    if (mode_initial_rdf || mode_single_run) return;

    T = _simulation.getEnsemble()->T();
    Log::global_log->info() << "[PMF] Target temperature = " << T << std::endl;

    reference_rdf.read(rdf_path);
    Log::global_log->info() << "[PMF] RDF has been read successfully" << std::endl;
    InitializePotentialValues();

    pairs_handler = new IBIPairsHandler();
    _simulation.setParticlePairsHandler(pairs_handler);
    _simulation.setCellProcessor(new LegacyCellProcessor(_simulation.getcutoffRadius(), _simulation.getLJCutoff(), pairs_handler));

    if (pot_path.empty()) pairs_handler->getPotentialFunction() = reference_potential;
    else pairs_handler->getPotentialFunction().read(pot_path);
    DerivativeOfPotential();

    Log::global_log->info() << "[PMF] LegacyCellProcessor set" << std::endl;
    Log::global_log->info() << "[PMF] ForcedAdapter Class being used" << std::endl;
    Log::global_log->info() << "[PMF] Enabled " << std::endl;
    Log::global_log->info() << "[PMF] Damping factor of "<< alpha << std::endl;
}

void PMF::afterForces(ParticleContainer* pc, DomainDecompBase* dd, unsigned long step) {
    if (mode_single_run) {
        current_steps++;

        switch (ibi_phase) {
            case FIRST_INIT: {
                if (current_steps < steps_measurement) { profiler.ProfileData(pc); break; }

                ibi_phase = EQUILIBRATE;
                current_steps = 0;
                ibi_iteration = 0;

                // create rdf_0
                std::vector<double> avg_rdf;
                profiler.GetRDFTotal(avg_rdf);

                reference_rdf.SetXValues(profiler.GetRNodes());
                reference_rdf.SetYValues(avg_rdf);
                reference_rdf.write(createFilepath("rdf"));

                profiler.ResetBuffers();

                // set up remaining fields
                T = _simulation.getEnsemble()->T();
                InitializePotentialValues();

                pairs_handler = new IBIPairsHandler();
                _simulation.setParticlePairsHandler(pairs_handler);
                _simulation.setCellProcessor(new LegacyCellProcessor(_simulation.getcutoffRadius(), _simulation.getLJCutoff(), pairs_handler));
                pairs_handler->getPotentialFunction() = reference_potential;
                DerivativeOfPotential();

                Log::global_log->info() << "[PMF] Transitioning from first initialization to equilibration with PMF" << std::endl;
                break;
            }
            case EQUILIBRATE: {
                if (current_steps < steps_equilibration) break;
                ibi_phase = MEASURE;
                current_steps = 0;
                // nothing else to do in this phase
                Log::global_log->info() << "[PMF] Transitioning from equilibration to measurement" << std::endl;
                break;
            }
            case MEASURE: {
                if (current_steps < steps_measurement) { profiler.ProfileData(pc); break; }

                ibi_phase = EQUILIBRATE;
                current_steps = 0;

                // update pot -> derivative
                AddPotentialCorrection();
                ibi_iteration++;
                pairs_handler->getPotentialFunction().write(createFilepath("pot"));
                DerivativeOfPotential();
                WriteRDF();
                Log::global_log->info() << "[PMF] Convergence: " << ConvergenceCheck() << std::endl;
                Log::global_log->info() << "[PMF] Transitioning from measurement to equilibration" << std::endl;
                profiler.ResetBuffers();
                break;
            }
        }

        return;
    }


    if (!mode_equilibrate) profiler.ProfileData(pc);
}

void PMF::finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {
    // everything else was already done during afterForces
    if (mode_single_run) {
        pairs_handler->getPotentialFunction().write("final_pot.txt");

        std::stringstream ss;
        for (double conv : convergence_steps) {
            ss << conv << " ";
        }
        Log::global_log->info() << "[PMF] Convergence steps: " << ss.str() << std::endl;
        return;
    }

    // only want to create RDF after this run
    if (mode_initial_rdf) {
        WriteRDF(); return;
    }

    // only want to update potential during this run
    if (!mode_equilibrate) {
        AddPotentialCorrection();

        ibi_iteration++;
        WriteRDF();
        pairs_handler->getPotentialFunction().write(createFilepath("pot"));
        Log::global_log->info() << "[PMF] Convergence: " << ConvergenceCheck() << std::endl;
    }
}

/********************
 * ****************** FUNCTIONS NOT FROM THE INTERFACE
 *******************/

void PMF::InitializePotentialValues() {
    reference_potential.SetXValues(reference_rdf.GetXValues());
    reference_potential.SetYValues(reference_rdf.GetYValues());

    std::vector<double>& u0 = reference_potential.GetYValues();
    for (double& u : u0) u = -T * std::log(u);

    // extrapolate function values (linearly) for which log was undefined in range of [0, r_min)
    // first find r_min
    int r_min_idx = -1;
    for (int idx = 0; idx < reference_potential.GetXValues().size(); idx++) {
        if (std::isfinite(reference_potential.GetYValues()[idx])) { r_min_idx = idx; break; }
    }
    // there are inf values -> need to do something
    // r_min_idx has the first position with a valid y value
    if (r_min_idx != -1) {
        const double dy = reference_potential.GetYValues()[r_min_idx + 1] - reference_potential.GetYValues()[r_min_idx];
        double y_val = reference_potential.GetYValues()[r_min_idx];
        for (int idx = r_min_idx; idx >= 0; idx--) {
            reference_potential.GetYValues()[idx] = y_val;
            y_val -= dy;
        }
    }

    reference_potential.write("pot_0.txt");
}

void PMF::AddPotentialCorrection() {
    std::vector<double> avg_rdf;
    profiler.GetRDFTotal(avg_rdf);

    std::vector<double>& pot = pairs_handler->getPotentialFunction().GetYValues();
    FunctionPL updateFunction {};
    updateFunction.SetXValues(pairs_handler->getPotentialFunction().GetXValues());
    updateFunction.GetYValues().resize(pot.size(), 0.0);
    for (int idx = 0; idx < pot.size(); idx++) {
        const double update = alpha * T * std::log(avg_rdf[idx] / reference_rdf.GetYValues()[idx]);
        if (std::isfinite(update)) {
            pot[idx] -= update;
            updateFunction.GetYValues()[idx] = update;
        }
    }
    updateFunction.write(createFilepath("update"));

    // extrapolate function values (linearly) for which log was undefined in range of [0, r_min)
    // first find r_min
    int r_min_idx = -1;
    for (int idx = 0; idx < reference_rdf.GetXValues().size(); idx++) {
        if (reference_rdf.GetYValues()[idx] != 0.0) { r_min_idx = idx; break; }
    }
    // r_min_idx has the first position with a valid y value
    if (r_min_idx != -1) {
        const double dy = pot[r_min_idx + 1] - pot[r_min_idx];
        double y_val = pot[r_min_idx];
        for (int idx = r_min_idx; idx >= 0; idx--) {
            pot[idx] = y_val;
            y_val -= dy;
        }
    }
}

void PMF::DerivativeOfPotential() {
    pairs_handler->getForceFunction() = pairs_handler->getPotentialFunction().Derivative(1e+6, 0.0);
    pairs_handler->getForceFunction() *= -1;
    pairs_handler->getForceFunction().write(createFilepath("force"));
}

double PMF::ConvergenceCheck() {
    std::vector<double>& g_0 = reference_rdf.GetYValues();
    std::vector<double> g_i;
    profiler.GetRDFTotal(g_i);

    double vec_norm = 0.0;
    for (int i = 0; i < g_0.size(); ++i) {
        vec_norm += std::pow(g_0[i] - g_i[i], 2.0);
    }

    convergence_steps.push_back(vec_norm);
    return vec_norm;
}

void PMF::WriteRDF() {
    std::vector<double> avg_rdf;
    profiler.GetRDFTotal(avg_rdf);

    FunctionPL rdfFunction {};
    rdfFunction.SetXValues(profiler.GetRNodes());
    rdfFunction.SetYValues(avg_rdf);
    rdfFunction.write(createFilepath("rdf"));
}

std::string PMF::createFilepath(const std::string &prefix) const {
    std::stringstream path;
    path << prefix << "_" << ibi_iteration << ".txt";
    return path.str();
}

