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
    double _hybridDim;
    //! @brief dimensions of the FP region
    std::array<double, 3> _dim;

    /**
     * @brief Constructs a FPRegion
     * @param low is the left lower corner of the FP box
     * @param high is the right higher corner of the FP box
     * @param hDim is the width of the hybrid area shell around the FP box
     * */
    FPRegion(const std::array<double,3>& low = {0,0,0}, const std::array<double, 3>& high = {0,0,0}, double hDim = 0)
            : _low(low), _high(high), _hybridDim(hDim), _dim({0,0,0}) {}

    //! @brief gets the lower corner of hybrid region
    [[nodiscard]] const std::array<double, 3>& getLowHybrid() const { return _lowHybrid; }

    //! @brief gets the upper corner of hybrid region
    [[nodiscard]] const std::array<double, 3>& getHighHybrid() const { return _highHybrid; }

    //! @brief initialized this FPRegion instance from XML data
    void readXML(XMLfileUnits &xmlconfig) {
        xmlconfig.getNodeValue("lowX", _low[0]);
        xmlconfig.getNodeValue("lowY", _low[1]);
        xmlconfig.getNodeValue("lowZ", _low[2]);
        xmlconfig.getNodeValue("highX", _high[0]);
        xmlconfig.getNodeValue("highY", _high[1]);
        xmlconfig.getNodeValue("highZ", _high[2]);
        xmlconfig.getNodeValue("hybridDim", _hybridDim);


        for(int d = 0; d < 3; d++) {
            _lowHybrid[d] = _low[d] - _hybridDim;
            _highHybrid[d] = _high[d] + _hybridDim;
            _center[d] = (_high[d] - _low[d])/2 + _low[d];
            _dim[d] = _high[d] - _low[d];
        }
    }

    /**
     * @brief Computes the intersection point of the line that is created by connection
     * the center of the region and the provided point with the either the inner or outer box of this region
     * depending on the intersection type.
     * @param point outside or inside of the box that is used for the computation
     * @param inter Intersection::H_FP inner box, Intersection::CG_H outer box
     * */
    std::array<double, 3> computeIntersection(std::array<double,3> point, Intersection inter) {
        double scale = std::max(
                            std::max(std::abs(point[0] - _center[0])/(_dim[0]+2*_hybridDim*inter)*2.0,
                                     std::abs(point[1] - _center[1])/(_dim[1]+2*_hybridDim*inter)*2.0),
                            std::abs(point[2] - _center[2])/(_dim[2]+2*_hybridDim*inter)*2.0
                        );
        std::array<double,3> result{0,0,0};
        for(int d = 0; d < 3; d++) result[d] = _center[d] + (point[d] - _center[d]) / scale;
        return result;
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
};

/**
 * Simple container for mesoscopic values.
 * During normal force computation of ls1-mardyn mesoscopic values are computed. In Hybrid regions those are wrong.
 * AdResS needs to recompute those. This struct stores those values.
 * */
struct MesoValues {
    //! @brief variable used to sum the virial contribution of all pairs
    double _virial;
    //! @brief variable used to sum the Upot6LJ contribution of all pairs
    double _upot6LJ;
    //! @brief variable used to sum the UpotXpoles contribution of all pairs
    double _upotXpoles;
    //! @brief variable used to sum the MyRF contribution of all pairs
    double _myRF;

    /**
     * @brief Constructs a MesoValues Container
     * @param v virial
     * @param lj upot6LJ
     * @param pole upotXpoles
     * @param rf myRF
     * */
    explicit MesoValues(double v = 0.f, double lj = 0.f, double pole = 0.f, double rf = 0.f) : _virial(v), _upot6LJ(lj), _upotXpoles(pole), _myRF(rf) {}

    /**
     * @brief sets all values to 0
     * */
    void clear() {
       _virial = 0;
       _upot6LJ = 0;
       _upotXpoles = 0;
       _myRF = 0;
    }

    /**
     * @brief Sets the stores mesoscopic values in the passed @param domain.
     * */
    void setInDomain(Domain* domain) const {
        domain->setLocalUpot(_upot6LJ / 6. + _upotXpoles + _myRF + domain->getLocalUpot());
        domain->setLocalVirial(_virial + 3.0 * _myRF + domain->getLocalVirial());
    }
};

#endif //MARDYN_ADRESSDATA_H
