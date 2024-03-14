//
// Created by alex on 23.05.23.
//

#ifndef MARDYN_ANALYSIS_H
#define MARDYN_ANALYSIS_H

#include <array>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <cmath>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "Molecule.h"

/**
 * Performs various kind of analysis on the input data
 * */
class Analysis {
public:
    struct Density {
        std::array<float,3> position;
        std::vector<double> densityValues;
    };

    explicit Analysis(const std::vector<Molecule>& data) : _data(data) , _molTypes(0) {
        std::unordered_set<uint32_t> CIDs;
        for(auto& mol : data) {
            CIDs.insert(mol.cid);
        }
        CIDs.insert(0);
        if(CIDs.size()%3==0) _molTypes = CIDs.size()/3;
        else _molTypes = CIDs.size();
    }

    ~Analysis() {}

    /**
     * Computes density information
     * slices up the domain according to the parameters.
     * @param nSamples how often should a sample be taken per dimension
     * @param boxWidth size of the box of one sample
     * @param bbox bounding box of complete simulation
     * @param check if check[i] is true, will slice domain in that dimension
     * */
    std::vector<Density> computeDensities(int nSamples, float boxWidth, std::array<float,3> bbox, std::array<int,3> check) {
        std::vector<Density> buffer;
        std::array<float,3> begin = { }, end = { }, stepSize = { }, halfBoxSize = { };
        std::array<int,3> steps { };
        double volume = 1;
        for(int d = 0; d < 3; d++) {
            halfBoxSize[d] = check[d] ? boxWidth/2 : bbox[d]/2;

            begin[d] = halfBoxSize[d];
            end[d] = bbox[d] - halfBoxSize[d];
            stepSize[d] = (end[d] - begin[d]) / static_cast<float>(nSamples - 1);
            volume *= check[d] ? boxWidth : bbox[d];
            steps[d] = check[d] ? nSamples : 1;
        }


        #if defined(_OPENMP)
        #pragma omp parallel for collapse(3)
        #endif
        for(int sz = 0; sz < steps[2]; sz++) {
            for(int sy = 0; sy < steps[1]; sy++) {
                for(int sx = 0; sx < steps[0]; sx++) {
                    Density density{};
                    std::array<float,3> center{};
                    center[0] = sx * stepSize[0] + begin[0];
                    center[1] = sy * stepSize[1] + begin[1];
                    center[2] = sz * stepSize[2] + begin[2];
                    density.position = center;

                    std::array<float,3> boxLow {};
                    std::array<float,3> boxHigh {};
                    for(int d = 0; d < 3; d++) {
                        boxLow[d] = center[d] - halfBoxSize[d];
                        boxHigh[d] = center[d] + halfBoxSize[d];
                    }

                    density.densityValues.resize(_molTypes);
                    std::for_each(_data.begin(), _data.end(), [&](const Molecule& m){
                        bool cond = true;
                        for(int d = 0; d < 3; d++) cond &= m.r.at(d) >= boxLow[d];
                        for(int d = 0; d < 3; d++) cond &= m.r.at(d) <= boxHigh[d];
                        if(!cond) return;

                        density.densityValues.at(m.cid/3) += 1;
                    });
                    for(int i = 0; i < _molTypes; i++) density.densityValues[i] /= volume;

                    #if defined(_OPENMP)
                    #pragma omp critical
                    #endif
                    {
                        buffer.emplace_back(std::move(density));
                    }
                }
            }
        }

        return buffer;
    }

private:
    //! @brief ref to data from vtk file
    const std::vector<Molecule>& _data;
    //! @brief total count of different molecule types (i.e. different components)
    uint64_t _molTypes;
};

#endif //MARDYN_ANALYSIS_H
