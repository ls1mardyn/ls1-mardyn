//
// Created by alex on 16.05.23.
//

#ifndef MARDYN_ADRESSREGIONTRAVERSAL_H
#define MARDYN_ADRESSREGIONTRAVERSAL_H

#include <array>
#include <unordered_map>
#include "particleContainer/ParticleContainer.h"
#include "AdResSData.h"
#include "AdResSForceAdapter.h"

/**
 * Handles the traversal of a region defined by two corner points, that do not overlap with the cell structure of any
 * particle container. This is needed during the force calculation in AdResS as FPRegions can be defined anywhere.
 * To increase performance this class also is capable of handling the traversal with multiple threads.
 * For that, there is a simplified C08 traversal implemented alongside with a cell structure consisting of region iterators.
 *
 * Based on C08CellPairTraversal.h
 * */
class AdResSRegionTraversal {
public:
    /**
     * @brief Constructor: Initializes the traversal region defined by checkLow and checkHigh
     * After this the region can be traversed.
     * @param checkLow Lower corner of box to be traversed
     * @param checkHigh Higher corner of box to be traversed
     * @param particleContainer global particle Container
     * @param compMap mapping of component IDs to their actual resolution
     * */
    AdResSRegionTraversal(std::array<double, 3> checkLow, std::array<double, 3> checkHigh,
                          ParticleContainer* particleContainer, std::unordered_map<unsigned long, Resolution>& compMap);

    void traverse(AdResSForceAdapter& forceAdapter, FPRegion& region, bool invert);

private:
    //! @brief Since we are working an cells structure on top of the existing one, this gives offsets to get neighbours
    std::array<std::pair<unsigned long, unsigned long>, 14> _cellPairOffsets8Pack;
    //! @brief global particle container
    ParticleContainer* _particleContainer;
    //! @brief our cell structure, when using this always copy the iterator. cell is defined by region of iterator.
    std::vector<RegionParticleIterator> _cells;
    //! @brief end index of C08 traversal, is for all threads valid
    std::array<unsigned long, 3> _end;
    //! @brief dimensions of cell structure in cells
    std::array<long, 3> _dims;
    //! @brief mapping of component ids to a AdResS resolution
    std::unordered_map<unsigned long, Resolution>& _comp_to_res;
    //! @brief cutoff distance
    double _cutoff;
    //! @brief cutoff squared
    double _cutoff2;
    //! @brief lennard jones cutoff distance squared
    double _ljcutoff2;

    void processCell(RegionParticleIterator& cell, AdResSForceAdapter& forceAdapter, FPRegion& region, bool invert);
    void processCellPair(RegionParticleIterator& cell1, RegionParticleIterator& cell2, AdResSForceAdapter& forceAdapter, FPRegion& region, bool invert);
};


#endif //MARDYN_ADRESSREGIONTRAVERSAL_H
