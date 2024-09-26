//
// Created by alex on 9/20/24.
//
#include "IBI_Math.h"
#include "molecules/FullMolecule.h"
#include "molecules/Component.h"

#include <cmath>

FunctionPL::FunctionPL(double def_low, double def_high): default_value_lower(def_low), default_value_upper(def_high) { }

void FunctionPL::SetXValues(const std::vector<double>& x) {
    x_values.resize(x.size());
    std::copy(x.begin(), x.end(), x_values.begin());

#ifndef NDEBUG
    bool uniform = true;
    const auto dx = x[1] - x[0];
    for (int idx = 1; idx < x.size() - 1; idx++) {
        const auto delta = x[idx + 1] - x[idx];
        uniform &= delta >= dx - 1e-15 || delta <= dx + 1e-15;
    }
    mardyn_assert(uniform);
#endif
}

void FunctionPL::SetYValues(const std::vector<double>& y) {
    y_values.resize(y.size());
    std::copy(y.begin(), y.end(), y_values.begin());
}

double FunctionPL::EvaluateAt(double x) {
    if (x > x_values[x_values.size() - 1]) return default_value_upper;
    if (x < x_values[0]) return default_value_lower;

    //find between which 2 values
    const auto idx_upper = FindUpperNode(x);
    mardyn_assert(idx_upper > 0 && idx_upper < x_values.size());
    const auto idx_lower = idx_upper -1;

    //interpolate
    const auto fx = LinearInterpolation(idx_lower, idx_upper, x);
    return fx;
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
    int i = 0;
    const auto i_max = x_values.size();
    while(x > x_values[i]) {
        i++;
        if(i == i_max) break;
    }
    return i;
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
