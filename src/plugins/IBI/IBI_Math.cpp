//
// Created by alex on 9/20/24.
//
#include "IBI_Math.h"
#include "molecules/FullMolecule.h"
#include "molecules/Component.h"

#include <cmath>

FunctionPL::FunctionPL(double def_low, double def_high): default_value_lower(def_low), default_value_upper(def_high), uniform_dx(false), uniform_dx_val(0.0) { }

void FunctionPL::SetXValues(const std::vector<double>& x) {
    x_values.resize(x.size());
    std::copy(x.begin(), x.end(), x_values.begin());

#ifndef NDEBUG
    mardyn_assert(checkUniformity());
#endif
}

void FunctionPL::SetYValues(const std::vector<double>& y) {
    y_values.resize(y.size());
    std::copy(y.begin(), y.end(), y_values.begin());
}

double FunctionPL::EvaluateAt(double x) {
    if (x >= x_values[x_values.size() - 1]) return default_value_upper;
    if (x < x_values[0]) return default_value_lower;

    //find between which 2 values
    const auto idx_upper = FindUpperNode(x);
    mardyn_assert(idx_upper > 0 && idx_upper < x_values.size());
    const auto idx_lower = idx_upper -1;

    //interpolate
    const auto fx = LinearInterpolation(idx_lower, idx_upper, x);
    return fx;
}

double FunctionPL::EvaluateAt(double x, int idx) {
    if (idx >= x_values.size()) return default_value_upper;
    if (idx <= 0) return default_value_lower;

    const auto idx_upper = idx;
    const auto idx_lower = idx_upper -1;

    //interpolate
    const auto fx = LinearInterpolation(idx_lower, idx_upper, x);
    return fx;
}

int FunctionPL::idxOf(double x) {
    if (x >= x_values[x_values.size() - 1]) return x_values.size();
    if (x < x_values[0]) return 0;

    return FindUpperNode(x);
}

FunctionPL FunctionPL::Derivative(double def_low, double def_high) const {
    FunctionPL result(def_low, def_high);
    result.x_values = x_values;
    const auto size = y_values.size();
    result.y_values.resize(size);
    // forward finite diff at start
    result.y_values[0] = (y_values[1] - y_values[0]) / (x_values[1] - x_values[0]);
    // backward finite diff at end
    result.y_values[size - 1] = (y_values[size - 1] - y_values[size - 2]) / (x_values[size - 1] - x_values[size - 2]);
    // central finite diff for rest
    for (int idx = 1; idx < size - 1; idx++) {
        result.y_values[idx] = (y_values[idx + 1] - y_values[idx - 1]) / (2 * (x_values[idx + 1] - x_values[idx - 1]));
    }
    if (uniform_dx) {
        result.uniform_dx = true;
        result.uniform_dx_val = uniform_dx_val;
    }
    return result;
}

double FunctionPL::LinearInterpolation(int a, int b, double x) {
    const auto ya = y_values[a];
    const auto yb = y_values[b];
    const auto xa = x_values[a];
    const auto xb = x_values[b];

    return ya + (yb - ya) / (xb - xa) * (x - xa);
}

int FunctionPL::FindUpperNode(double x) {
    if (uniform_dx) {
        const int bin = (x - x_values[0]) / uniform_dx_val;
        return bin + 1;
    } else {
        int i = 0;
        const auto i_max = x_values.size();
        while(x > x_values[i]) {
            i++;
            if(i == i_max) break;
        }
        return i;
    }
}

void FunctionPL::read(const std::string &path) {
    std::ifstream file{path};
    if (!file) {
        Log::global_log->error() << "Could not read the file" << std::endl;
        throw std::runtime_error("Failed to open file.");
    }

    double n1, n2;
    while (file >> n1 >> n2) {
        x_values.push_back(n1);
        y_values.push_back(n2);
    }
    file.close();

    checkUniformity();
}

void FunctionPL::write(const std::string &path) {
    std::ofstream file{path};
    if (!file) {
        Log::global_log->error() << "Could not read the file" << std::endl;
        throw std::runtime_error("Failed to open file.");
    }

    bool begin = true;
    for (int idx = 0; idx < x_values.size(); idx++) {
        if (begin) begin = false;
        else file << "\n";

        file << x_values[idx] << "\t" << y_values[idx];
    }
    file.close();
}

bool FunctionPL::checkUniformity() {
    bool uniform = true;
    const auto dx = x_values[1] - x_values[0];
    for (int idx = 1; idx < x_values.size() - 1; idx++) {
        const auto delta = x_values[idx + 1] - x_values[idx];
        uniform &= delta >= dx - 1e-10 && delta <= dx + 1e-10;
    }

    uniform_dx = uniform;
    uniform_dx_val = dx;

    return uniform;
}

std::string IBIOptimizer::createFilepath(const std::string &prefix, int step) const {
    std::stringstream path;
    path << prefix << "_" << step << ".txt";
    return path.str();
}

void IBIOptimizer::extrapolate() {
    // extrapolate function values (linearly) for which log was undefined in range of [0, r_min)
    // first find r_min
    int r_min_idx = -1;
    for (int idx = 0; idx < _ref_rdf.GetXValues().size(); idx++) {
        if (_ref_rdf.GetYValues()[idx] != 0.0) { r_min_idx = idx; break; }
    }
    // r_min_idx has the first position with a valid y value
    if (r_min_idx != -1) {
        const double dy = _pot.GetYValues()[r_min_idx + 1] - _pot.GetYValues()[r_min_idx];
        double y_val = _pot.GetYValues()[r_min_idx];
        for (int idx = r_min_idx; idx >= 0; idx--) {
            _pot.GetYValues()[idx] = y_val;
            y_val -= dy;
        }
    }
}

DefaultOptimizer::DefaultOptimizer(FunctionPL &pot, FunctionPL &ref_pot, FunctionPL &ref_rdf, double T, double alpha)
: IBIOptimizer(pot, ref_pot, ref_rdf, T), _updateFunction(), _alpha(alpha) {
    _updateFunction.SetXValues(ref_rdf.GetXValues());
    _updateFunction.GetYValues().resize(ref_pot.GetYValues().size(), 0.0);
}

void DefaultOptimizer::step(const std::vector<double> &rdf, int ibi_step) {
    for (int idx = 0; idx < _pot.GetXValues().size(); idx++) {
        const double update = _alpha * _temp * std::log(rdf[idx] / _ref_rdf.GetYValues()[idx]);
        if (std::isfinite(update)) {
            _pot.GetYValues()[idx] += update;
            _updateFunction.GetYValues()[idx] = update;
        }
    }
    _updateFunction.write(createFilepath("update", ibi_step));

    extrapolate();
}

AdamOptimizer::AdamOptimizer(
        FunctionPL &pot, FunctionPL &ref_pot, FunctionPL &ref_rdf,
        double T, double alpha, double beta1, double beta2, double eps) :
        IBIOptimizer(pot, ref_pot, ref_rdf, T), _m(), _v(), _grad(), _update(),
        _alpha(alpha), _beta1(beta1), _beta2(beta2), _eps(eps) {
    auto& nodes = ref_rdf.GetXValues();
    _m.SetXValues(nodes);
    _v.SetXValues(nodes);
    _grad.SetXValues(nodes);
    _update.SetXValues(nodes);

    const auto size = nodes.size();
    _m.GetYValues().resize(size, 0.0);
    _v.GetYValues().resize(size, 0.0);
    _grad.GetYValues().resize(size, 0.0);
    _update.GetYValues().resize(size, 0.0);
}

void AdamOptimizer::step(const std::vector<double> &rdf, int ibi_step) {
    computeGradient(rdf);
    updateMean();
    updateVar();

    const double m_div = (1 - std::pow(_beta1, ibi_step+1));
    const double v_div = (1 - std::pow(_beta2, ibi_step+1));
    //compute and add update term
    for (int idx = 0; idx < _update.GetYValues().size(); idx++) {
        //add bias correction
        const double m_hat = _m.GetYValues()[idx] / m_div;
        const double v_hat = _v.GetYValues()[idx] / v_div;

        const double update = - _alpha * m_hat / (std::sqrt(v_hat) + _eps);
        _update.GetYValues()[idx] = update;
        _pot.GetYValues()[idx] += update;
    }
    _update.write(createFilepath("update", ibi_step));

    extrapolate();
}

void AdamOptimizer::computeGradient(const std::vector<double> &rdf) {
    for (int idx = 0; idx < rdf.size(); idx++) {
        double grad = 0;
        double tmp = (2.0 / (_temp)) * (std::log(_ref_rdf.GetYValues()[idx]) - std::log(rdf[idx]));
        if (std::isfinite(tmp)) grad = tmp;
        _grad.GetYValues()[idx] = grad;
    }
}

void AdamOptimizer::updateMean() {
    for (int idx = 0; idx < _m.GetYValues().size(); idx++) {
        _m.GetYValues()[idx] = _beta1 * _m.GetYValues()[idx] + (1 - _beta1) * _grad.GetYValues()[idx];
    }
}

void AdamOptimizer::updateVar() {
    for (int idx = 0; idx < _v.GetYValues().size(); idx++) {
        const double grad = _grad.GetYValues()[idx];
        _v.GetYValues()[idx] = _beta2 * _v.GetYValues()[idx] + (1 - _beta2) * (grad * grad);
    }
}

//===============================================================
// UTIL FUNCTION
//===============================================================

double Distance2(const std::array<double, 3> &p1, const std::array<double, 3> &p2) {
    double r = 0.0;
    for(int i = 0; i < 3; i++) {
        r += std::pow(p2[i] - p1[i], 2);
    }
    return r;
}

double Distance(const std::array<double, 3> &p1, const std::array<double, 3> &p2) {
    return std::sqrt(Distance2(p1, p2));
}

void Normalize(std::array<double, 3> &v) {
    const auto length = Distance(v, {0, 0, 0});
    for (int dim = 0; dim < 3; dim++) v[dim] /= length;
}

std::array<double, 3> GetCOM(const Molecule &molecule) {
    std::array<double,3> com{0.0,0.0,0.0};
    Component* comp = molecule.component();
    const double total_mass = comp->m();
    for (int lj = 0; lj < comp->numLJcenters(); lj++) {
        LJcenter& lj_center = comp->ljcenter(lj);
        for (int dim = 0; dim < 3; dim++) {
            com[dim] += lj_center.m() * molecule.ljcenter_d_abs(lj)[dim];
        }
    }

    for (int i = 0; i < 3; i++) com[i] /= total_mass;
    return com;
}
