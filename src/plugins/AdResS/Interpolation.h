//
// Created by alex on 04.01.24.
//

#ifndef MARDYN_INTERPOLATION_H
#define MARDYN_INTERPOLATION_H

#include <vector>

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
}
#endif //MARDYN_INTERPOLATION_H
