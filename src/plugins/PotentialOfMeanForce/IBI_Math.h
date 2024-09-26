//
// Created by alex on 9/20/24.
//

#ifndef MARDYN_IBI_MATH_H
#define MARDYN_IBI_MATH_H

#include <array>
#include <vector>
#include <string>
#include "molecules/MoleculeForwardDeclaration.h"

/**
 * Piecewise linear function
 * */
class FunctionPL {
public:
    explicit FunctionPL(double def_low = 0.0, double def_high = 0.0);
    void SetXValues(const std::vector<double>& x);
    void SetYValues(const std::vector<double>& y);

    std::vector<double>& GetXValues() { return x_values; };
    std::vector<double>& GetYValues() { return y_values; };
    [[nodiscard]] const std::vector<double>& GetXValues() const { return x_values; };
    [[nodiscard]] const std::vector<double>& GetYValues() const { return y_values; };
    [[nodiscard]] double GetLowerDefault() const { return default_value_lower; }
    [[nodiscard]] double GetUpperDefault() const { return default_value_upper; }
    double EvaluateAt(double x);
    [[nodiscard]] FunctionPL Derivative(double def_low = 0.0, double def_high = 0.0) const;

    template<typename T>
    FunctionPL& operator*=(T val) { for (double& y : y_values) y *= val; return *this; }
    template<typename T>
    FunctionPL& operator+=(T val) { for (double& y : y_values) y *= val; return *this; }

    /// read from file
    void read(const std::string& path);
    /// write to file
    void write(const std::string& path);
private:
    /// Interpolates between two points
    double LinearInterpolation(int a, int b, double x);
    /// Returns the index of the next large x node
    int FindUpperNode(double x);
    /// function value nodes
    std::vector<double> x_values;
    /// function values
    std::vector<double> y_values;
    /// default value lower
    double default_value_lower;
    /// default value upper
    double default_value_upper;
};

//===============================================================
// UTIL FUNCTION
//===============================================================

/// Returns squared distance between 2 points
[[maybe_unused]] double Distance2(const std::array<double,3>& p1, const std::array<double,3>& p2);
/// Returns L2 distance between 2 points
[[maybe_unused]] double Distance(const std::array<double,3>& p1, const std::array<double,3>& p2);
/// Normalizes vector to unit length
[[maybe_unused]] void Normalize(std::array<double,3>& v);
/// Returns the center of mass of the molecule
std::array<double, 3> GetCOM(const Molecule& molecule);

#endif //MARDYN_IBI_MATH_H
