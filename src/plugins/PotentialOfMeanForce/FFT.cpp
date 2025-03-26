//
// Created by alex on 3/22/25.
//

#include "FFT.h"

std::vector<double> FFT::low_pass_filter(const std::vector<double> &data) {
    // Compute FFT
    std::vector<std::complex<double>> freq_domain = rfft(data, 8);

    // Compute inverse FFT
    return irfft(freq_domain, data.size());
}

std::vector<std::complex<double>> FFT::rfft(const std::vector<double> &rho, int f_max) {
    int N = rho.size();
    std::vector<std::complex<double>> freq_domain(f_max, std::complex<double>{0.0, 0.0});

    const std::complex<double> j {0.0, 1.0};
    std::vector<std::complex<double>> omegas (f_max, std::complex<double>{0.0, 0.0});
    for (int k = 0; k < f_max; k++) omegas[k] = std::exp(-j * 2.0 * M_PI * (double) k / (double) N);
    std::vector<std::complex<double>> active_omegas {omegas}; // contains with power 1

    // n = 0
    for (int k = 0; k < f_max; k++) freq_domain[k] += rho[0];

    // n = 1...N
    for (int n = 1; n < N; n++) {
        const auto r_n = rho[n];

        for (int k = 0; k < f_max; k++) freq_domain[k] += r_n * active_omegas[k];
        for (int k = 0; k < f_max; k++) active_omegas[k] *= omegas[k];
    }

    return freq_domain;
}

std::vector<double> FFT::irfft(const std::vector<std::complex<double>> &freq_domain, int N) {
    int f_max = freq_domain.size();
    std::vector<double> rho(N, 0.0);
    std::vector<std::complex<double>> c_rho(N, std::complex<double>());

    const std::complex<double> j {0.0, 1.0};
    std::vector<std::complex<double>> omegas (f_max, std::complex<double>{0.0, 0.0});
    for (int k = 0; k < f_max; k++) omegas[k] = std::exp(j * 2.0 * M_PI * (double) k / (double) N);
    std::vector<std::complex<double>> active_omegas {omegas}; // contains with power 1

    // n = 0
    for (int k = 0; k < f_max; k++) c_rho[0] += freq_domain[k];

    // n = 1...N
    for (int n = 1; n < N; n++) {
        for (int k = 0; k < f_max; k++) c_rho[n] += freq_domain[k] * active_omegas[k];
        for (int k = 0; k < f_max; k++) active_omegas[k] *= omegas[k];
    }

    for (int n = 0; n < N; n++) {
        rho[n] = c_rho[n].real() / N;
    }

    return rho;
}
