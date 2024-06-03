//
// Created by alex on 04.01.24.
//

#include "utils/mardyn_assert.h"
#include "Interpolation.h"
#include "Simulation.h"
#include "parallel/DomainDecompBase.h"

#include <array>
#include <vector>
#include <math.h>

using Log::global_log;

void Interpolation::Function::writeXML(const std::string &filename) const {
#ifdef ENABLE_MPI
	if (_simulation.domainDecomposition().getRank() != 0) return;
#endif
	try {
		std::ofstream file {filename};
		file << "<forceFunction>\n";
		file << "    <startX>" << begin << "</startX>\n";
		for(unsigned long i = 0; i < n; i++) {
			file << "    <samplePoint id=\"" << i + 1 <<"\">\n";
			file << "    " << "    <grad>" << gradients[i] << "</grad>\n";
			file << "    " << "    <func>" << function_values[i] << "</func>\n";
			if(i < n - 1) file << "    " << "    <step>" << step_width[i] << "</step>\n";
			file << "    </samplePoint>\n";
		}
		file << "</forceFunction>\n";
		file.flush();
		file.close();
	} catch (std::ifstream::failure& e) {
		global_log->error() << "[Interpolation] Failed to write Interpolation function.\n" << e.what() << std::endl;
		_simulation.exit(-1);
	}
}

void Interpolation::Function::writeTXT(const std::string &filename) const {
#ifdef ENABLE_MPI
	if (_simulation.domainDecomposition().getRank() != 0) return;
#endif
	try {
		std::ofstream file {filename};
		file << "IF1\n";
		file << "start\t" << begin << "\n";
		file << "num_points\t" << n << "\n";
		file << "#p point_idx\tgrad\tf_val\tstep_width\n";
		for(unsigned long i = 0; i < n; i++) {
			file << "p\t" << i+1 << "\t" << gradients[i] << "\t" << function_values[i];
			if (i < n - 1) file << "\t" << step_width[i] << "\n";
			else file << "\n";
		}
		file.flush();
		file.close();
	} catch (std::ifstream::failure& e) {
		global_log->error() << "[Interpolation] Failed to write Interpolation function.\n" << e.what() << std::endl;
		_simulation.exit(-1);
	}
}

void Interpolation::Function::loadTXT(const std::string &filename) {
	try {
		std::ifstream file {filename};
		std::string s_buf;
		int i_buf;

		if (file >> s_buf; s_buf != "IF1") throw std::runtime_error("wrong format!");
		if (file >> s_buf; s_buf != "start") throw std::runtime_error("wrong format!");
		file >> begin;
		if (file >> s_buf; s_buf != "num_points") throw std::runtime_error("wrong format!");
		file >> n;

		while (!file.eof() && file.good()) {
			file >> s_buf;
			if (s_buf[0] == '#' || s_buf[0] != 'p') { file.ignore(2048, '\n'); continue; }

			file >> i_buf;
			i_buf--;
			if (i_buf < 0 || i_buf >= n) throw std::runtime_error("bad idx!");
			file >> gradients[i_buf] >> function_values[i_buf];
			if (i_buf < n - 1) file >> step_width[i_buf];
		}
		file.close();
	} catch (std::ifstream::failure& e) {
		global_log->error() << "[Interpolation] Failed to load Interpolation function.\n" << e.what() << std::endl;
		_simulation.exit(-1);
	}
}

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

[[maybe_unused]] void Interpolation::solveTriDiagonalMatrix(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, std::vector<double> &x,
                                           std::vector<double> &d) {
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

[[maybe_unused]] void Interpolation::computeHermite(double begin, std::vector<double> &fVals, std::vector<double> &steps, int samples,
                                   Function &fun) {
    // the base idea is to create the hermite matrix
    // this here is the simplified version, in which uniform step size is assumed
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

[[maybe_unused]] static inline double gmm_kernel(double x, double x_i, double xi) {
	// see https://arxiv.org/pdf/1504.07351.pdf and https://pubs.acs.org/doi/10.1021/jp909219k
	double scale = std::pow(2*M_PI*xi*xi, -3.0/2.0);
	double exp = std::exp(-0.5 * std::pow((x-x_i)/xi, 2));
	return scale*exp;
}

void
Interpolation::createGMM(double begin, double end, int samples, double xi, const std::vector<double>& centers, Interpolation::Function &function) {
	std::vector<double> hermite_f_val;
	std::vector<double> hermite_steps;

	auto step_width = static_cast<double>((end-begin)/samples);
	hermite_f_val.resize(samples, 0.0);
	hermite_steps.resize(samples-1, step_width);

	// every center corresponds to one exp function and must be evaluated at all roots
	#if defined(_OPENMP)
	# pragma omp parallel for
	#endif
	for(int step = 0; step < samples; step++)
	{
		double pos = step * step_width;

		for(auto& center : centers) {
			hermite_f_val[step] += gmm_kernel(pos, center, xi);
		}
	}

	computeHermite(begin, hermite_f_val, hermite_steps, samples, function);
}

void Interpolation::realFT(const std::vector<double>& R, unsigned int k_max, double T, std::vector<std::complex<double>>& output) {
	static const std::complex<double> j{0.0, 1.0};
	output.resize(k_max+1);

	const auto fac = -j * 2.0 * M_PI / T;

	#if defined(_OPENMP)
	#pragma omp parallel for
	#endif
	for(unsigned long i_w = 0; i_w < k_max+1; i_w++) {

		std::complex<double> tmp{0.0, 0.0};
		for(unsigned long i_r = 0; i_r < R.size(); i_r++) {
			tmp += std::exp(fac * static_cast<double>(i_w) * R[i_r]);
		}
		output[i_w] = tmp;
	}
}

void Interpolation::ift(const std::vector<std::complex<double>>& F, unsigned int k_max, double T, double begin, double end, int samples, Function& function) {
	std::vector<double> hermite_f_val;
	std::vector<double> hermite_steps;

	auto step_width = static_cast<double>((end-begin)/samples);
	hermite_f_val.resize(samples, 0.0);
	hermite_steps.resize(samples-1, step_width);

	const auto fac = 2.0 * M_PI / T;

	#if defined(_OPENMP)
	#pragma omp parallel for
	#endif
	for(unsigned long i_n = 0; i_n < samples; i_n++) {

		double tmp = F[0].real();
		double n = i_n * step_width;

		for(unsigned long k = 1; k <= k_max; k++) {
			const auto arg = fac * k * n;
			tmp += 2 * F[k].real() * std::cos(arg);
			tmp -= 2 * F[k].imag() * std::sin(arg);
		}

		hermite_f_val[i_n] = tmp;
	}

	computeHermite(begin, hermite_f_val, hermite_steps, samples, function);
}

void Interpolation::filterFT(std::vector<std::complex<double>> &F) {
	/*for(unsigned long i = F.size()/16; i < F.size(); i++) {
		F[i] = 0;
	}*/
	for(unsigned long k = 0; k < F.size(); k++) {
		F[k] *= std::exp(-0.5 * std::pow(k/15.0, 2.0));
	}
}
