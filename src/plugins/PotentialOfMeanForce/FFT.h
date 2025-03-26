//
// Created by alex on 3/22/25.
//

#ifndef MARDYN_FFT_H
#define MARDYN_FFT_H

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>

class FFT {
public:
    static std::vector<double> low_pass_filter(const std::vector<double>& data);

private:
    static std::vector<std::complex<double>> rfft(const std::vector<double>& rho, int f_max);

    static std::vector<double> irfft(const std::vector<std::complex<double>>& freq_domain, int N);
};


#endif //MARDYN_FFT_H
