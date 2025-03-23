//
// Created by alex on 3/22/25.
//

#ifndef MARDYN_FFT_H
#define MARDYN_FFT_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fftw3.h>
#include <complex>

class FFT {
public:
    static std::vector<double> low_pass_filter(const std::vector<double>& data);

private:
    static std::vector<std::complex<double>> rfft(const std::vector<double>& rho);

    static std::vector<double> irfft(const std::vector<std::complex<double>>& freq_domain, int N);

    static void apply_low_pass_filter(std::vector<std::complex<double>>& freq_domain, int N, double cutoff);

    static void apply_frequency_filter(std::vector<std::complex<double>>& freq_domain, int keep_bins);
};


#endif //MARDYN_FFT_H
