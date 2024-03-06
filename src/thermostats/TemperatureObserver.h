//
// Created by alex on 06.03.24.
//

#ifndef MARDYN_TEMPERATUREOBSERVER_H
#define MARDYN_TEMPERATUREOBSERVER_H


#include "utils/xmlfileUnits.h"
#include "utils/mardyn_assert.h"

#include <array>
#include <vector>


class Domain;
class ParticleContainer;
class Ensemble;

/**
 * Ring buffer that stores up to n entries. Will override oldest entry once n+1 elements are inserted.
 * */
template<typename T>
struct NBuffer {
    NBuffer() = default;
    [[maybe_unused]] explicit NBuffer(std::size_t n) : _n(n), _begin(-1), _end(0) {
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

/**
 * Data struct for computing SAM and CAM averages
 * */
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

/**
 * Measures the temperature in selected regions. See Karimian et al. 2011 and 2014
 * */
class TemperatureObserver {
public:
    /**
     * Init all buffers
     * */
    void init();
    void readXML(XMLfileUnits &xmlconfig);
    /**
     * Measure temp in all regions, updates caches and stores current temp
     * */
    void step(ParticleContainer* particleContainer);

    /**
     * @returns all regions, for which the temp is being measured.
     * */
    void getRegions(std::vector<std::pair<std::array<double, 3>, std::array<double, 3>>>& regions);

    /**
     * Gets the temperature for the selected region.
     * Regions can be queried by accessing getRegions().
     * @returns T * k_b
     * */
    double getTemperature(std::size_t index);
private:
    std::size_t _max_samples;

    using d3 = std::array<double, 3>;
    struct box_t {
        box_t() = default;
        box_t(const d3& l, const d3& h) : low(l), high(h) {}
        d3 low, high;
    };
    struct box_data_t {
        BinData history{ };
        double current_temp{ };
    };

    //! @brief each region is defined by a box in 3d space (box_t) and a data cache (box_data_t)
    std::vector<std::pair<box_t, box_data_t>> _region_data;

    /**
     * Based on: Details about pressure calculation in molecular dynamic analysis by Karimian et al. 2014
     * returns T * k_b
     * */
    double measureTemp(const std::array<double,3>& low, const std::array<double, 3>& high, ParticleContainer *particleContainer, BinData& binData);

    /**
     * Computes the SAM over the CAM average of the mean velocity for the specified region. (SCM method)
     * */
    std::array<double, 3> computeMeanVelocityStep(const std::array<double,3>& low, const std::array<double, 3>& high, ParticleContainer *particleContainer, BinData& binData);
};

#endif //MARDYN_TEMPERATUREOBSERVER_H
