//
// Created by alex on 06.05.23.
//

#ifndef MARDYN_ADRESSDATA_H
#define MARDYN_ADRESSDATA_H

#include <array>
#include "Domain.h"

/**
 * Implemented AdResS runs on 3 different resolutions.
 * FullParticle -> molecules are represented with full detail
 * Hybrid -> molecules in transition area are FullParticle and CoarseGrain at the same time
 * CoarseGrain -> simplified molecule representation with lower detail
 * */
enum Resolution{
    FullParticle = 0, Hybrid = 1, CoarseGrain = 2, ResolutionCount
};

/**
 * FPRegion is the Full Particle region. Currently it is a box starting in _low and goes until _high.
 * The hybrid area is a shell around the FP box. Each FPRegion can have different dimensions for the hybrid area.
 * */
struct FPRegion {
    /**
     * Defines a boundary region.
     * H_FP is the region between the hybrid and full particle region.
     * CG_H is the region between the hybrid and coarse grain region.
     * */
    enum Intersection { H_FP = 0, CG_H = 1};

    //! @brief front left lower corner of FP
    std::array<double, 3> _low;
    //! @brief rear right upper corner of FP
    std::array<double, 3> _high;
    //! @brief front left lower corner of Hybrid
    std::array<double, 3> _lowHybrid;
    //! @brief rear right upper corner of Hybrid
    std::array<double, 3> _highHybrid;
    //! @brief center point of FP region
    std::array<double, 3> _center;
    //! @brief thickness of hybrid walls
    std::array<double, 3> _hybridDims;
    //! @brief dimensions of the FP region
    std::array<double, 3> _dim;

    /**
     * @brief Constructs a FPRegion
     * @param low is the left lower corner of the FP box
     * @param high is the right higher corner of the FP box
     * @param hDim is the width of the hybrid area shell around the FP box
     * */
    FPRegion(const std::array<double,3>& low = {0,0,0}, const std::array<double, 3>& high = {0,0,0}, const std::array<double,3>& dims = {0,0,0})
            : _low(low), _high(high), _hybridDims(dims), _dim({0,0,0}) {}

    //! @brief gets the lower corner of hybrid region
    [[nodiscard]] const std::array<double, 3>& getLowHybrid() const { return _lowHybrid; }

    //! @brief gets the upper corner of hybrid region
    [[nodiscard]] const std::array<double, 3>& getHighHybrid() const { return _highHybrid; }

    //! @brief computes all fields according to input data low,high,dims
    void init() {
        for(int d = 0; d < 3; d++) {
            _lowHybrid[d] = _low[d] - _hybridDims[d];
            _highHybrid[d] = _high[d] + _hybridDims[d];
            _center[d] = (_high[d] - _low[d])/2 + _low[d];
            _dim[d] = _high[d] - _low[d];
        }
    }

    //! @brief initialized this FPRegion instance from XML data
    void readXML(XMLfileUnits &xmlconfig) {
        xmlconfig.getNodeValue("lowX", _low[0]);
        xmlconfig.getNodeValue("lowY", _low[1]);
        xmlconfig.getNodeValue("lowZ", _low[2]);
        xmlconfig.getNodeValue("highX", _high[0]);
        xmlconfig.getNodeValue("highY", _high[1]);
        xmlconfig.getNodeValue("highZ", _high[2]);
        xmlconfig.getNodeValue("hybridDimX", _hybridDims[0]);
        xmlconfig.getNodeValue("hybridDimY", _hybridDims[1]);
        xmlconfig.getNodeValue("hybridDimZ", _hybridDims[2]);


        init();
    }

    /**
     * @brief checks if the point is in the provided box defined by low and high
     * @param point is the point to be checked
     * @param low is the inclusive lower left front corner of the box
     * @param high is the exclusive upper right rear corner of the box
     * @returns true when point is in the box
     * */
    static bool isInnerPoint(std::array<double,3> point, std::array<double, 3> low, std::array<double, 3> high) {
        bool res = true;
        for(int d = 0; d < 3; d++) res &= point[d] >= low[d];
        for(int d = 0; d < 3; d++) res &= point[d] < high[d];
        return res;
    }

    /**
     * @brief checks if the given point is inside the specified region, which is either the FP or H box defined by this.
     * When this region goes beyond the domain bounds, it should wrap around as the simulation uses periodic bounds.
     * This method also checks for that by using the passed domain.
     * @param domain to check periodic bounds
     * @param region either FullParticle or Hybrid
     * @param point to check
     * */
    bool isInnerPointDomain(Domain* domain, Resolution region, std::array<double, 3> point) const {
        std::array<double,3> globLen{0};
        for(int d = 0; d < 3; d++) {
            globLen[d] = domain->getGlobalLength(d);
        }

        std::array<bool, 3> inDim {false};
        bool inDomain = true;
        for(int d = 0; d < 3; d++) {
            inDim[d] = _lowHybrid[d] >= 0 && _highHybrid[d] <= globLen[d];
            inDomain &= inDim[d];
        }

        //special case for interface situation
        std::array<double,3> borderOffset{0,0,0};
        for(int d = 0; d < 3; d++) {
            borderOffset[d] = (_hybridDims[d] == 0) ? _simulation.getcutoffRadius() : 0;
        }

        //region is within domain bounds
        if(inDomain && region == FullParticle) return isInnerPoint(point, {_low[0] - borderOffset[0], _low[1] - borderOffset[1], _low[2] - borderOffset[2]}, {_high[0] + borderOffset[0], _high[1] + borderOffset[1], _high[2] + borderOffset[2]});
        if(inDomain && region == Hybrid) return isInnerPoint(point, {_lowHybrid[0] - borderOffset[0], _lowHybrid[1] - borderOffset[1], _lowHybrid[2] - borderOffset[2]}, {_highHybrid[0] + borderOffset[0], _highHybrid[1] + borderOffset[1], _highHybrid[2] + borderOffset[2]});

        //region crosses bound
        bool checkOuter = region == Hybrid;
        std::array<std::array<double,2>,3> deltaLowHighDim{};
        for(int d = 0; d < 3; d++) {
            deltaLowHighDim[d][0] = std::max(-(_low[d] - checkOuter * _hybridDims[d]), 0.);
            deltaLowHighDim[d][1] = std::max((_high[d] + checkOuter * _hybridDims[d]) - globLen[d], 0.);
        }

        // check all regions that wrap around due to periodic bounds
        // we only need to check for each dim twice if it actually wraps around in that dimension
        // in the arrays we create the indices of the boxes depending on the viewed case
        // for x y or z: if they are 0 then we do not check a wrap around in that dimension
        bool result = false;
        for(int x = 0; x <= !inDim[0]; x++) {
            for(int y = 0; y <= !inDim[1]; y++) {
                for(int z = 0; z <= !inDim[2]; z++) {
                    std::array<double, 3> checkLow {
                        (1 - x) * std::max(_low[0] - checkOuter * _hybridDims[0], 0.) + x * ((deltaLowHighDim[0][0] != 0) * (globLen[0] - deltaLowHighDim[0][0]) + (deltaLowHighDim[0][1] != 0) * 0),
                        (1 - y) * std::max(_low[1] - checkOuter * _hybridDims[1], 0.) + y * ((deltaLowHighDim[1][0] != 0) * (globLen[1] - deltaLowHighDim[1][0]) + (deltaLowHighDim[1][1] != 0) * 0),
                        (1 - z) * std::max(_low[2] - checkOuter * _hybridDims[2], 0.) + z * ((deltaLowHighDim[2][0] != 0) * (globLen[2] - deltaLowHighDim[2][0]) + (deltaLowHighDim[2][1] != 0) * 0) };

                    std::array<double, 3> checkHigh {
                        (1 - x) * std::min(_high[0] + checkOuter * _hybridDims[0], globLen[0]) + x * ((deltaLowHighDim[0][0] != 0) * (globLen[0]) + (deltaLowHighDim[0][1] != 0) * deltaLowHighDim[0][1]),
                        (1 - y) * std::min(_high[1] + checkOuter * _hybridDims[1], globLen[1]) + y * ((deltaLowHighDim[1][0] != 0) * (globLen[1]) + (deltaLowHighDim[1][1] != 0) * deltaLowHighDim[1][1]),
                        (1 - z) * std::min(_high[2] + checkOuter * _hybridDims[2], globLen[2]) + z * ((deltaLowHighDim[2][0] != 0) * (globLen[2]) + (deltaLowHighDim[2][1] != 0) * deltaLowHighDim[2][1]) };
                    result |= isInnerPoint(point, checkLow, checkHigh);
                }
            }
        }
        return result;
    }
};

#endif //MARDYN_ADRESSDATA_H
