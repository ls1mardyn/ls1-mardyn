#include"IBI.h"
#include "particleContainer/adapter/LegacyCellProcessor.h"


void IBI::readXML(XMLfileUnits& xmlfile) {
    xmlfile.getNodeValue("singleRun", mode_single_run);
    // using this double tmp buffer here to allow for easier configuration in xml file
    double tmp = 1e+6;
    xmlfile.getNodeValue("stepsEquil", tmp);
    steps_equilibration = static_cast<int>(tmp);
    tmp = 1e+5;
    xmlfile.getNodeValue("stepsMeasure", tmp);
    steps_measurement = static_cast<int>(tmp);

    std::array<const std::string, 2> opt_types {"default", "adam"};
    std::string opt_type = opt_types[0];
    double opt_alpha = 0.001, opt_beta1 = 0.9, opt_beta2 = 0.999, opt_eps = 1e-8;
    xmlfile.getNodeValue("optimizer/type", opt_type);
    xmlfile.getNodeValue("optimizer/alpha", opt_alpha);
    xmlfile.getNodeValue("optimizer/beta1", opt_beta1);
    xmlfile.getNodeValue("optimizer/beta2", opt_beta2);
    xmlfile.getNodeValue("optimizer/eps", opt_eps);
    if (opt_type == opt_types[1]) optConfig.type = OptConfig::ADAM;
    else optConfig.type = OptConfig::DEFAULT;
    optConfig.alpha = opt_alpha;
    optConfig.beta1 = opt_beta1;
    optConfig.beta2 = opt_beta2;
    optConfig.eps = opt_eps;

    int bins;
    xmlfile.getNodeValue("bins", bins);
    profiler.init(bins);

    xmlfile.getNodeValue("equilibrate", mode_equilibrate);

    xmlfile.getNodeValue("rdfPath", rdf_path);
    xmlfile.getNodeValue("potPath", pot_path);

    double conv_threshold = 0.01;
    std::string conv_method = "l2";
    std::string stop_method = "worse";
    int window_size = 10;
    int ignore = 0;
    xmlfile.getNodeValue("convergence/threshold", conv_threshold);
    xmlfile.getNodeValue("convergence/method", conv_method);
    xmlfile.getNodeValue("convergence/stop-method", stop_method);
    xmlfile.getNodeValue("convergence/window-size", window_size);
    xmlfile.getNodeValue("convergence/offset", ignore);
    ConvergenceCheck.init(conv_threshold, conv_method, stop_method, window_size, ignore);

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

void IBI::init(ParticleContainer* pc, DomainDecompBase* domainDecomp, Domain* domain) {
    if (mode_initial_rdf || mode_single_run) return;

    // check if we have more than one component
    if (_simulation.getEnsemble()->getComponents()->size() != 1) throw std::runtime_error("IBI does not support multi-component systems.");

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
    //Log::global_log->info() << "[PMF] Damping factor of "<< alpha << std::endl;
}

void IBI::afterForces(ParticleContainer* pc, DomainDecompBase* dd, unsigned long step) {
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
                CreateSwapComponent();
                CreateOptimizer();
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

                auto conv = ConvergenceCheck(reference_rdf, profiler);
                Log::global_log->info() << "[PMF] Convergence: target_reached=" << conv.first << " value=" << conv.second << std::endl;
                const bool should_stop = ConvergenceCheck.ShouldStop();
                if (should_stop) {
                    Log::global_log->info() << "[PMF] Convergence: stopping criterion reached. Stopping now." << std::endl;
                    _simulation.setNumTimesteps(_simulation.getSimulationStep());
                    break;
                }

                Log::global_log->info() << "[PMF] Transitioning from measurement to equilibration" << std::endl;
                profiler.ResetBuffers();
                break;
            }
        }

        return;
    }


    if (!mode_equilibrate) profiler.ProfileData(pc);
}

void IBI::finish(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain) {
    // everything else was already done during afterForces
    if (mode_single_run) {
        pairs_handler->getPotentialFunction().write("final_pot.txt");

        std::stringstream ss;
        ConvergenceCheck.LogValues(ss);
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
        auto conv = ConvergenceCheck(reference_rdf, profiler);
        Log::global_log->info() << "[PMF] Convergence: target_reached=" << conv.first << "value=" << conv.second << std::endl;
    }
}

/********************
 * ****************** FUNCTIONS NOT FROM THE INTERFACE
 *******************/

void IBI::InitializePotentialValues() {
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

void IBI::AddPotentialCorrection() {
    std::vector<double> avg_rdf;
    profiler.GetRDFTotal(avg_rdf);

    optimizer->step(avg_rdf, ibi_iteration);
}

void IBI::DerivativeOfPotential() {
    pairs_handler->getForceFunction() = pairs_handler->getPotentialFunction().Derivative(1e+6, 0.0);
    pairs_handler->getForceFunction() *= -1;
    pairs_handler->getForceFunction().write(createFilepath("force"));
}

void IBI::WriteRDF() {
    std::vector<double> avg_rdf;
    profiler.GetRDFTotal(avg_rdf);

    FunctionPL rdfFunction {};
    rdfFunction.SetXValues(profiler.GetRNodes());
    rdfFunction.SetYValues(avg_rdf);
    rdfFunction.write(createFilepath("rdf"));
}

void IBI::CreateSwapComponent() {
    auto* components = _simulation.getEnsemble()->getComponents();
    // we have asserted during init that we only have one active component
    if (components->at(0).numSites() == 1 && components->at(0).numLJcenters() == 1) return; // no need to change

    //=================================================
    // Create new Component
    Component comp {};
    comp.setID(1);
    comp.setName("IBI-Comp");

    LJcenter site {};
    site.setR(0, 0.0);
    site.setR(1, 0.0);
    site.setR(2, 0.0);
    double mass = 0;
    for (int i = 0; i < components->at(0).numLJcenters(); i++) mass += components->at(0).ljcenter(i).m();
    site.setM(mass);
    std::string site_name = "LJ126";
    site.setName(site_name);

    comp.addLJcenter(site);
    components->push_back(comp);
    _simulation.getEnsemble()->setComponentLookUpIDs();

    //=================================================
    // Replace comp pointer with new comp
    Component* new_comp = &components->at(1);
    for (auto it = _simulation.getMoleculeContainer()->iterator(ParticleIterator::ALL_CELLS); it.isValid(); ++it) {
        it->setComponent(new_comp);
    }
}

std::string IBI::createFilepath(const std::string &prefix) const {
    std::stringstream path;
    path << prefix << "_" << ibi_iteration << ".txt";
    return path.str();
}

void IBI::CreateOptimizer() {
    if (optConfig.type == OptConfig::ADAM) {
        optimizer = std::make_unique<AdamOptimizer>(pairs_handler->getPotentialFunction(), reference_potential,
                                                    reference_rdf, T, optConfig.alpha, optConfig.beta1,
                                                    optConfig.beta2, optConfig.eps);
    }
    else {
        optimizer = std::make_unique<DefaultOptimizer>(pairs_handler->getPotentialFunction(),
                                                       reference_potential, reference_rdf, T, optConfig.alpha);
    }
}

//=============================================================
//  CONVERGENCE CHECK
//=============================================================

void IBI::Convergence::init(double th, const std::string &mode_str, const std::string &stop_str, int window, int ignore) {
    if (mode_str == "l2") mode = L2;
    else mode = INTEGRAL;
    threshold = th;
    if (stop_str == "worse") stopping_mode = ON_WORSE;
    else stopping_mode = WINDOW;
    window_size = window;
    ignore_counter = ignore;
}

std::pair<bool, double> IBI::Convergence::integral(const FunctionPL &ref, RDFProfiler& profiler) {
    const std::vector<double>& x = ref.GetXValues();
    const std::vector<double>& g_0 = ref.GetYValues();
    std::vector<double> g_i;
    profiler.GetRDFTotal(g_i);

    // integrate diff
    double diff = 0;
    for (int idx = 0; idx < x.size()-1; idx++) {
        const double dx = x[idx+1] - x[idx];
        const double lower = std::abs(g_i[idx] - g_0[idx]);
        const double upper = std::abs(g_i[idx+1] - g_0[idx+1]);
        diff += dx * (upper + lower) / 2.0;
    }

    // integrate sum
    double sum = 0;
    for (int idx = 0; idx < x.size()-1; idx++) {
        const double dx = x[idx+1] - x[idx];
        const double lower = std::abs(g_i[idx] + g_0[idx]);
        const double upper = std::abs(g_i[idx+1] + g_0[idx+1]);
        sum += dx * (upper + lower) / 2.0;
    }

    if (sum == 0) sum = 1e-15;
    const double conv_value = 1.0 - (diff / sum);
    conv_values.push_back(conv_value);
    return {conv_value <= threshold, conv_value};
}

std::pair<bool, double> IBI::Convergence::l2(const FunctionPL &ref, RDFProfiler& profiler) {
    const std::vector<double>& g_0 = ref.GetYValues();
    std::vector<double> g_i;
    profiler.GetRDFTotal(g_i);

    double vec_norm = 0.0;
    for (int i = 0; i < g_0.size(); ++i) {
        vec_norm += std::pow(g_0[i] - g_i[i], 2.0);
    }

    conv_values.push_back(vec_norm);
    return {vec_norm <= threshold, vec_norm};
}

void IBI::Convergence::LogValues(std::ostream &ostream) {
    for (double conv : conv_values) {
        ostream << conv << " ";
    }
}

std::pair<bool, double> IBI::Convergence::operator()(const FunctionPL &ref, RDFProfiler& profiler) {
    if (mode == L2) return l2(ref, profiler);
    else if (mode == INTEGRAL) return integral(ref, profiler);
    else throw std::runtime_error("Unknown convergence method");
}

bool IBI::Convergence::ShouldStop() {
    if (0 < ignore_counter--) return false;

    const auto steps = conv_values.size();

    if (stopping_mode == ON_WORSE) {
        if (steps < 2) return false;
        if (mode == INTEGRAL) {
            if (conv_values[steps-2] < threshold) return false;
            if (conv_values[steps-2] <= conv_values[steps-1]) return false;
            return true;
        }
        if (mode == L2) {
            if (conv_values[steps-2] > threshold) return false;
            if (conv_values[steps-2] >= conv_values[steps-1]) return false;
            return true;
        }
    }

    else if (stopping_mode == WINDOW) {
        if (steps < window_size) return false;
        bool all_in_range = true;
        for (int i = 0; i < window_size; i++) {
            if (mode == INTEGRAL) all_in_range &= conv_values[steps - 1 - i] >= threshold;
            if (mode == L2) all_in_range &= conv_values[steps - 1 - i] <= threshold;
        }
        return all_in_range;
    }
    throw std::runtime_error("Unknown stopping method");
}
