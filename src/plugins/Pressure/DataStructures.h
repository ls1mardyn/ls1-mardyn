//
// Created by alex on 27.02.24.
//

#ifndef MARDYN_DATASTRUCTURES_H
#define MARDYN_DATASTRUCTURES_H
#include "utils/mardyn_assert.h"

#include <array>
#include <vector>

template<typename T>
struct NBuffer {
    NBuffer() = default;
    explicit NBuffer(std::size_t n) : _n(n), _begin(-1), _end(0) {
        _data.resize(n);
    }
    void init(std::size_t n) {
        _n = n;
        _begin = -1;
        _end = 0;
        _data.resize(n);
    }
    void insert(const T& t) {
        if(_begin == -1) {
            _data[0] = t;
            _begin = 0;
            _end = 1;
        }
        else if (_begin == _end) {
            _data[_end] = t;
            _begin = (_begin+1) % _n;
            _end = _begin;
        }
        else {
            _data[_end] = t;
            _end = (_end+1) % _n;
        }
    }
    const T& get(std::size_t index) {
        mardyn_assert(index < _n);
        return _data[(_begin + index) % _n];
    }
    std::size_t size() {
        if(_begin == -1) {
            return 0;
        }
        else if (_begin == _end) {
            return _n;
        }
        else {
            return _end - _begin;
        }
    }
    std::size_t _n{};
    long _begin{};
    long _end{};
    std::vector<T> _data;
};

struct BinData {
    NBuffer<std::size_t> mol_counts;
    NBuffer<std::array<double,3>> velocities;
    NBuffer<std::array<double,3>> v_cam;

    void init(std::size_t n) {
        mol_counts.init(n);
        velocities.init(n);
        v_cam.init(n);
    }
};

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
