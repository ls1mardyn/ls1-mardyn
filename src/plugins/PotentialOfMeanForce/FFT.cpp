//
// Created by alex on 3/22/25.
//

#include "FFT.h"

std::vector<double> FFT::low_pass_filter(const std::vector<double> &data) {
    // Compute FFT
    std::vector<std::complex<double>> freq_domain = rfft(data);

    // Apply low-pass filter with cutoff frequency = 8
    apply_frequency_filter(freq_domain, 8);

    // Compute inverse FFT
    return irfft(freq_domain, data.size());
}

std::vector<std::complex<double>> FFT::rfft(const std::vector<double> &rho) {
    int N = rho.size();
    int N_half = N / 2 + 1;  // Size of complex output

    std::vector<std::complex<double>> freq_domain(N_half);

    // FFTW arrays
    fftw_plan plan;
    double* in = fftw_alloc_real(N);
    fftw_complex* out = fftw_alloc_complex(N_half);

    // Copy input data
    std::copy(rho.begin(), rho.end(), in);

    // Create FFT plan
    plan = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
    fftw_execute(plan);

    // Store result
    for (int i = 0; i < N_half; i++) {
        freq_domain[i] = std::complex<double>(out[i][0], out[i][1]);
    }

    // Cleanup
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return freq_domain;
}

std::vector<double> FFT::irfft(const std::vector<std::complex<double>> &freq_domain, int N) {
    int N_half = N / 2 + 1;
    std::vector<double> rho(N);

    // FFTW arrays
    fftw_plan plan;
    double* out = fftw_alloc_real(N);
    fftw_complex* in = fftw_alloc_complex(N_half);

    // Copy input data
    for (int i = 0; i < N_half; i++) {
        in[i][0] = freq_domain[i].real();
        in[i][1] = freq_domain[i].imag();
    }

    // Create inverse FFT plan
    plan = fftw_plan_dft_c2r_1d(N, in, out, FFTW_ESTIMATE);
    fftw_execute(plan);

    // Normalize and store result
    for (int i = 0; i < N; i++) {
        rho[i] = out[i] / N;  // Normalize by N
    }

    // Cleanup
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return rho;
}

void FFT::apply_low_pass_filter(std::vector<std::complex<double>> &freq_domain, int N, double cutoff) {
    double freq_resolution = 1.0 / N;  // Frequency resolution
    int N_half = N / 2 + 1;

    for (int i = 0; i < N_half; i++) {
        double freq = i * freq_resolution;
        if (freq > cutoff) {
            freq_domain[i] = {0.0, 0.0};  // Zero out high frequencies
        }
    }
}

void FFT::apply_frequency_filter(std::vector<std::complex<double>> &freq_domain, int keep_bins) {
    int N_half = freq_domain.size();

    // Zero out all frequency bins above `keep_bins`
    for (int i = keep_bins; i < N_half; i++) {
        freq_domain[i] = {0.0, 0.0};
    }
}
