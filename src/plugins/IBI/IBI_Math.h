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
    double EvaluateAt(double x, int idx);
    int idxOf(double x);
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
    /// Checks whether the x values have uniform dx
    bool checkUniformity();
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
    /// has uniform x deltas?
    bool uniform_dx;
    /// uniform dx value
    double uniform_dx_val;
};

class IBIOptimizer {
public:
    IBIOptimizer(FunctionPL& pot, FunctionPL& ref_pot, FunctionPL& ref_rdf, double T) : _pot(pot), _ref_pot(ref_pot), _ref_rdf(ref_rdf), _temp(T) {}
    virtual ~IBIOptimizer() = default;
    /**
     * @brief performs the actual optimization step
     * @param rdf current rdf
     * */
    virtual void step(const std::vector<double>& rdf, int ibi_step) = 0;

protected:
    [[nodiscard]] std::string createFilepath(const std::string& prefix, int step) const;
    void extrapolate();

    /// potential function
    FunctionPL& _pot;
    /// reference potential
    FunctionPL& _ref_pot;
    /// reference rdf
    FunctionPL& _ref_rdf;
    /// Target Temperature
    const double _temp;
};

/**
 * Implements the default update process for IBI
 * */
class DefaultOptimizer final : public IBIOptimizer {
public:
    DefaultOptimizer(FunctionPL& pot, FunctionPL& ref_pot, FunctionPL& ref_rdf, double T, double alpha);
    ~DefaultOptimizer() final = default;

    void step(const std::vector<double> &rdf, int ibi_step) override;

private:
    FunctionPL _updateFunction;
    double _alpha;
};

class AdamOptimizer final : public IBIOptimizer {
public:
    AdamOptimizer(FunctionPL& pot, FunctionPL& ref_pot, FunctionPL& ref_rdf, double T, double alpha, double beta1, double beta2, double eps);
    ~AdamOptimizer() final = default;

    void step(const std::vector<double> &rdf, int ibi_step) override;
private:
    void computeGradient(const std::vector<double> &rdf);
    void updateMean();
    void updateVar();

    /// means -> accumulate velocity
    FunctionPL _m;
    /// variances -> stop high frequency oscillation
    FunctionPL _v;
    /// current gradient
    FunctionPL _grad;
    /// update term
    FunctionPL _update;
    /// update rate
    const double _alpha;
    /// retention rate mean
    const double _beta1;
    /// retention rate variance
    const double _beta2;
    /// variance zero offset
    const double _eps;
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
