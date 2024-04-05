//
// Created by alex on 04.01.24.
//

#ifndef MARDYN_INTERPOLATION_H
#define MARDYN_INTERPOLATION_H

#include <vector>
#include <string>

namespace Interpolation {
    /**
    * Represents one interpolated function using Cubic Hermite Splines
    * */
    struct Function {
        //! @brief number of knots
        unsigned long n;
        //! @brief distance between each knot
        std::vector<double> step_width;
        //! @brief starting point of samples, x_0
        double begin;
        //! @brief samples of f(x)
        std::vector<double> function_values;
        //! @brief samples of f'(x)
        std::vector<double> gradients;
		//! @brief write this function to file in XML format
		void writeXML(const std::string &filename);
    };

    /**
     * Matrix NxM N rows, M columns
     * Row-First Memory alignment
     * */
    struct Matrix {
        Matrix() : _dim0(0), _dim1(0) {}
        Matrix(unsigned long dim0, unsigned long dim1) : _dim0(dim0), _dim1(dim1) {
            _data.resize(dim0 * dim1, 0.0);
        }

        /**
         * @param vec single column
         * @param repeat how often vec should be repeated horizontally
         * */
        Matrix(const std::vector<double>& vec, unsigned long repeat) : Matrix(vec.size(), repeat) {
            #if defined(_OPENMP)
            #pragma omp parallel for collapse(2)
            #endif
            for(unsigned long n = 0; n < _dim0; n++) {
                for(unsigned long m = 0; m < _dim1; m++) {
                    setAt(n, m, vec[n]);
                }
            }
        }

        /**
         * @param vec single row
         * @param repeat how often vec should be repeated vertically
         * */
        Matrix(unsigned long repeat, const std::vector<double>& vec) : Matrix(repeat, vec.size()) {
            #if defined(_OPENMP)
            #pragma omp parallel for collapse(2)
            #endif
            for(unsigned long n = 0; n < _dim0; n++) {
                for(unsigned long m = 0; m < _dim1; m++) {
                    setAt(n, m, vec[m]);
                }
            }
        }

        std::vector<double> operator*(const std::vector<double>& vec) const {
            mardyn_assert((vec.size() == _dim1));
            std::vector<double> result;
            result.resize(_dim0, 0.0);
            auto* raw_result = std::data(result);

            #if defined(_OPENMP)
            #pragma omp parallel for simd reduction(+:raw_result[:_dim0]) collapse(2)
            #endif
            for(unsigned long n = 0; n < _dim0; n++) {
                for(unsigned long m = 0; m < _dim1; m++) {
                    raw_result[n] += vec[m] * _data[n * _dim1 + m];
                }
            }
            return result;
        }

        void setAt(unsigned long n, unsigned long m, double value) {
            _data[n * _dim1 + m] = value;
        }

        double getAt(unsigned long n, unsigned long m) {
            return _data[n * _dim1 + m];
        }

        unsigned long _dim0;
        unsigned long _dim1;
        std::vector<double> _data;
    };

    /**
     * Numerically computes the gradient of the input vector. Uses finite difference coefficients (based on Lagrange Polynomials).
     * Input and output have equal size. The output on the borders use forward or backward differences respectively, the rest central.
     * @param input sample points of f(x)
     * @param output sample points of f'(x)
     * */
    [[maybe_unused]] void computeGradient(std::vector<double>& input, std::vector<double>& output);

    /**
     * Solves the equation system Ax=d, where A is a tridiagonal matrix composed of the diagonals a, b, and c.
     * The lengths of a,b, and c should be equal. Pad with 0.
     * Solves it in O(n).
     * Implementation from: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
     * @param a lower diagonal (0, a_2 to a_n) padded with 0 in front
     * @param b middle diagonal (b_1 to b_n)
     * @param c upper diagonal (c_1 to c_(n-1), 0) padded with 0 in the end
     * @param x output vector
     * @param d right vector
     * */
    [[maybe_unused]] void solveTriDiagonalMatrix(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c,
                                std::vector<double>& x, std::vector<double>& d);

    /**
     * Evaluates the Cubic Hermite interpolation spline at position x.
     * @param x evaluation point
     * @param fun function representation
     * */
    [[maybe_unused]] double computeHermiteAt(double x, Function& fun);

    /**
     * Creates the Cubic Hermite interpolation spline based on the specified knots and stores them in fun.
     * fun will take ownership of knots and steps.
     * Assumes that the boundary gradients are zero.
     * @param begin starting point of knots, x_0
     * @param fVals sample points of f(x)
     * @param steps distance between each knot, i.e. x_i - x_(i+1), must be of size: samples - 1.
     * @param samples total number of samples
     * @param fun output buffer
     * */
    [[maybe_unused]] void computeHermite(double begin, std::vector<double>& fVals, std::vector<double>& steps, int samples, Function& fun);

    /**
     * Integrates f by integrating each spline piece of f symbolically.
     * Constant offset is assumed to be zero, but can be added manually later.
     * @param f function f(x)
     * @param F integral function F(x), with F'(x) = f(x)
     * */
    [[maybe_unused]] void computeIntegral(Function& f, Function& F);

    /**
     * Computes the derivative of F piece-wise.
     * @param F function F(x)
     * @param f derivative of F, with f(x) = F'(x)
     * */
    [[maybe_unused]] void computeGradient(Function& F, Function& f);

    /**
     * Gaussian Kernel for kernel smoothing. The term 2b**2 was replaced with sigma**2
     * */
    [[maybe_unused]] inline double gaussian_kernel(double x, double x_i, double sigma);

    /**
     * Generates as matrix to be used for smoothing sampled data, if sampled data is y=f(x)
     * @param begin begin of x
     * @param end end of x
     * @param step_width sampling step width
     * @param sigma filter strength
     * @param output output
     * */
    [[maybe_unused]] void createGaussianMatrix(double begin, double end, double step_width, double sigma, Matrix& output);

    /**
     * Evaluates the provided function in the range begin:step_width:end to in- or decrease resolution.
     * Results is stored back into same function.
     * @param begin begin inclusive
     * @param end end exclusive
     * @param step_width distance between each step
     * @param function in/out function
     * */
    [[maybe_unused]] void resampleFunction(double begin, double end, double step_width, Function& function);
}
#endif //MARDYN_INTERPOLATION_H
