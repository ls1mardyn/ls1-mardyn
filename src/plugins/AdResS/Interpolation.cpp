//
// Created by alex on 04.01.24.
//

#include "utils/mardyn_assert.h"
#include "Interpolation.h"
#include "Simulation.h"

#include <array>
#include <vector>

using Log::global_log;

/**
 * @param begin starting index of folding
 * @param end end index of folding (exclusive)
 * @param write_offset, writes with offset 0 to same position as begin
 * */
template<std::size_t N>
static inline void fold(const std::array<double, N>& mask, std::vector<double>& input, std::vector<double>& output, int begin, int end, int write_offset) {
    for(int i = begin; i < end; i++) {
        double tmp = 0.0;
        for(int j = 0; j < N; j++) {
            tmp += mask[j] * input[i + j];
        }
        output[i + write_offset] = tmp;
    }
}

[[maybe_unused]] void Interpolation::computeGradient(std::vector<double> &input, std::vector<double> &output) {
    static const std::array<double, 5> c_central = {1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0};
    static const std::array<double, 3> c_forward = {-3.0/2.0, 2.0, -1.0/2.0};
    static const std::array<double, 3> c_backward = {1.0/2.0, -2.0, 3.0/2.0};

    if(input.size() < 5) {
        global_log->fatal() << "[AdResS] gradient computation requires at least 5 sample points" << std::endl;
        Simulation::exit(670);
    }

    output.resize(input.size(), 0.0);
    fold(c_forward, input, output, 0, 2, 0);
    fold(c_backward, input, output, input.size()-1-2-1, input.size()-2, 2);
    fold(c_central, input, output, 0, input.size()-4, 2);
}

[[maybe_unused]] void Interpolation::solveTriDiagonalMatrix(vector<double> &a, vector<double> &b, vector<double> &c, vector<double> &x,
                                           vector<double> &d) {
    std::size_t n = b.size();
    mardyn_assert(n == a.size());
    mardyn_assert(n == c.size());
    mardyn_assert(n == d.size());
    x.resize(n, 0.0);

    for (int i = 1; i < n; i++) {
        double w = a[i] / b[i - 1];
        b[i] = b[i] - w * c[i - 1];
        d[i] = d[i] - w * d[i - 1];
    }
    x[n-1] = d[n-1] / b[n-1];
    for (int i = n - 1 - 1; i >= 0; i--) {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    }
}

static inline double bernstein_3(double t, int k) {
    static const double bc[4] = {1, 3, 3, 1};
    mardyn_assert(t >= 0.0);
    mardyn_assert(t <= 1.0);
    mardyn_assert(k >= 0);
    mardyn_assert(k <= 3);

    double t_k = pow(t, k);
    double inv_t_k = pow((1.0-t), (3-k));
    return bc[k] * t_k * inv_t_k;
}

[[maybe_unused]] double Interpolation::computeHermiteAt(double x, Function &fun) {
    //get active spline and conv x to t
    int c_step = 0;
    double knot_begin = 0;
    {
        double c_pos = x - fun.begin;
        while(true) {
            if(c_step >= fun.n - 2) break;
            if(fun.step_width[c_step] + knot_begin > c_pos) break;
            knot_begin += fun.step_width[c_step];
            c_step++;
        }
    }

    mardyn_assert((c_step >= 0) && (c_step < (fun.n-1)));
    double t = (x - knot_begin - fun.begin) / fun.step_width[c_step];

    //cache bernstein values
    double b0 = bernstein_3(t, 0);
    double b1 = bernstein_3(t, 1);
    double b2 = bernstein_3(t, 2);
    double b3 = bernstein_3(t, 3);

    //compute base functions
    double h00 = b0 + b1;
    double h10 = 1.0/3.0 * b1;
    double h01 = b3 + b2;
    double h11 = -1.0/3.0 * b2;

    double result = h00 * fun.function_values[c_step] +
                    h01 * fun.function_values[c_step + 1] +
                    h10 * fun.step_width[c_step] * fun.gradients[c_step] +
                    h11 * fun.step_width[c_step] * fun.gradients[c_step + 1];

    return result;
}

[[maybe_unused]] void Interpolation::computeHermite(double begin, vector<double> &fVals, vector<double> &steps, int samples,
                                   Function &fun) {
    // the base idea is to create the hermite matrix
    // this here is the simplified version, in which uniform step size is assume
    // in the approach beneath the 141 matrix is replaced by the appropriate step sizes
    // 4 1   | y1'            y2 - y0 - h/3 * y0'
    // 1 4 1 | y2'  =  3/h *  y_i+2 - y_i
    //   1 4 | y3'            y4 - y2 - h/3 * y4'
    std::vector<double> a, b, c, x, d;
    a.resize(samples - 2, 1.0);
    b.resize(samples - 2, 4.0);
    c.resize(samples - 2, 1.0);
    d.resize(samples - 2, 0.0);
    a[0] = 0.0;
    c[samples - 2 - 1] = 0.0;
    for(int i = 0; i < samples - 2; i++) {
        d[i] = 3 * (fVals[i + 2] - fVals[i]);
        b[i] = 2 * (steps[i] + steps[i+1]);
    }
    for(int i = 0; i < samples - 2 - 1; i++) {
        a[i+1] = steps[i+1];
        c[i] = steps[i+1];
    }

    solveTriDiagonalMatrix(a, b, c, x, d);
    fun.begin = begin;
    fun.step_width = std::move(steps);
    fun.n = samples;
    fun.function_values = std::move(fVals);
    fun.gradients = std::move(x);

    //insert gradients x_0 and x_n (equal to 0)
    fun.gradients.resize(fun.gradients.size()+2, 0.0); //make two larger
    std::rotate(fun.gradients.rbegin(), fun.gradients.rbegin() + 1, fun.gradients.rend()); //rotate right for correct position
}

[[maybe_unused]] void Interpolation::computeIntegral(Interpolation::Function &f, Interpolation::Function &F) {
    unsigned long samples = f.n;
    F.begin = f.begin;
    F.n = samples;
    F.step_width = std::vector<double>(f.step_width);
    F.gradients = std::vector<double>(f.function_values);

    //Assign function values
    {
        F.function_values.resize(samples);
        F.function_values[0] = 0.0;

        double running_integral = 0.0;
        for(int i = 0; i < samples-1; i++) {
            running_integral += 0.5 * f.function_values[i]
                                 + 0.5 * f.function_values[i+1]
                                 + (1.0/12.0) * f.gradients[i] * f.step_width[i]
                                 - (1.0/12.0) * f.gradients[i+1] * f.step_width[i];
            F.function_values[i+1] = running_integral;
        }
    }
}

[[maybe_unused]] void Interpolation::computeGradient(Interpolation::Function &F, Interpolation::Function &f) {
    std::vector<double> steps = std::vector<double>(F.step_width);
    std::vector<double> fVals = std::vector<double>(F.gradients);
    computeHermite(F.begin, fVals, steps, F.n, f);
}

[[maybe_unused]] inline double Interpolation::gaussian_kernel(double x, double x_i, double sigma) {
    return std::exp(-(x - x_i)*(x - x_i) / (sigma*sigma));
}

void Interpolation::createGaussianMatrix(double begin, double end, double step_width, double sigma,
                                         Interpolation::Matrix &output) {
    std::vector<double> positions;
    auto count = static_cast<unsigned long>((end-begin)/step_width);
    positions.resize(count);

    {
        double pos = begin;
        for(double & position : positions) {
            position = pos;
            pos += step_width;
        }
    }

    Matrix cols {positions, count};
    Matrix rows {count, positions};
    output._dim0 = count;
    output._dim1 = count;
    output._data.resize(count * count);

    #if defined(_OPENMP)
    #pragma omp parallel for simd collapse(2)
    #endif
    for(unsigned long n = 0; n < count; n++) {
        for(unsigned long m = 0; m < count; m++) {
            output.setAt(n, m, gaussian_kernel(cols.getAt(n, m), rows.getAt(n, m), sigma));
        }
    }

    //normalize dim1 to sum up to 1
    std::vector<double> sums;
    sums.resize(count, 0.0);
    auto* sum_data = std::data(sums);
    #if defined(_OPENMP)
    #pragma omp parallel for reduction(+:sum_data[:count]) collapse(2)
    #endif
    for(unsigned long n = 0; n < count; n++) {
        for(unsigned long m = 0; m < count; m++) {
            sum_data[n] += output.getAt(n, m);
        }
    }

    #if defined(_OPENMP)
    #pragma omp parallel for collapse(2)
    #endif
    for(unsigned long n = 0; n < count; n++) {
        for(unsigned long m = 0; m < count; m++) {
            output.setAt(n, m, output.getAt(n, m) / sums[n]);
        }
    }
}

void Interpolation::resampleFunction(double begin, double end, double step_width, Interpolation::Function &function) {
    double pos = begin;
    std::vector<double> fVals;
    std::vector<double> steps;
    auto count = static_cast<unsigned long>((end - begin) / step_width);
    fVals.resize(count);
    steps.resize(count-1, step_width);

    for(unsigned long i = 0; i < count; i++) {
        fVals[i] = computeHermiteAt(pos, function);
        pos += step_width;
    }

    Function f;
    computeHermite(begin, fVals, steps, count, f);
    function = std::move(f);
}
