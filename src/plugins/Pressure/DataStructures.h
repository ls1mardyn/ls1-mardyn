//
// Created by alex on 27.02.24.
//

#ifndef MARDYN_DATASTRUCTURES_H
#define MARDYN_DATASTRUCTURES_H
#include "utils/mardyn_assert.h"
#include "thermostats/TemperatureObserver.h"

#include <array>
#include <vector>

template<typename T>
struct Vec3D {
    std::vector<T> data;
    std::array<std::size_t, 3> dims;
    T& at(std::size_t i, std::size_t j, std::size_t k) {
        return data[k * dims[0] * dims[1] + j * dims[0] + i];
    }
};

struct BinVec3D : public Vec3D<BinData> {
    void init(const std::array<std::size_t, 3>& d, std::size_t samples) {
        dims = d;
        data.resize(d[0]*d[1]*d[2]);
        for(auto& bd : data) bd.init(samples);
    }
};

#endif //MARDYN_DATASTRUCTURES_H
