//
// Created by alex on 09.05.23.
//

#ifndef MARDYN_INTERACTIONLOGPROCESSOR_H
#define MARDYN_INTERACTIONLOGPROCESSOR_H


#include <vector>
#include "particleContainer/adapter/CellProcessor.h"
#include "plugins/AdResS/AdResSData.h"

/**
 * Traverses all cells and maps all cells, that interact with each other while at least
 * one of those cells is part of a FPRegion, to the corresponding region.
 * This is used in AdResS to clear forces that were incorrectly computed during normal cell traversal.
 * Once all cells have been analyzed, results will be cached as FP regions do not move.
 * */
class InteractionLogProcessor : public CellProcessor {
public:
    /**
     * @brief Constructor
     * @param cutoffRadius no clue todo Alex check this
     * @param ljCutoffRadius noclue todo Alex check this
     * @param fpRegions reference to all FP regions from the AdResS plugin instance
     * */
    explicit InteractionLogProcessor(const double cutoffRadius, const double ljCutoffRadius,
                                     const std::vector<FPRegion> &fpRegions);
    //! @brief does nothing
    ~InteractionLogProcessor() override;
    //! @brief does nothing
    void initTraversal() override;
    //! @brief does nothing
    void preprocessCell(ParticleCell &cell) override;
    //! @brief checks if either of the cells is part of a FP region. If so, then those are inserted into the map.
    void processCellPair(ParticleCell &cell1, ParticleCell &cell2, bool sumAll) override;
    //! @brief checks if the cell is part of a FP region. If so, then it is inserted into the map.
    void processCell(ParticleCell &cell) override;
    //! @brief does nothing
    double processSingleMolecule(Molecule *m1, ParticleCell &cell2) override;
    //! @brief does nothing
    void postprocessCell(ParticleCell &cell) override;
    //! @brief set up caching state after first full traversal
    void endTraversal() override;
    //! @returns whether the results are already available or still have to be computed
    bool isCached() {return _processed;}
    //! @returns an immuatable reference to the created map
    const auto& getFPR2CellMap() { return _fprID_to_cells; }
    //! @returns an immuatable reference to the created map
    const auto& getFPR2CellPairMap() { return _fprID_to_cellpairs; }

private:
    //! @brief reference to all FP regions from AdResS plugin instance
    const std::vector<FPRegion>& _fpRegions;

    //! @brief map of index of FP region to all related Particle Cells
    std::unordered_map<unsigned long, std::set<ParticleCell*>> _fprID_to_cells;

    //! @brief map of index of FP region to all related Particle Cell pairs
    std::unordered_map<unsigned long, std::set<std::pair<ParticleCell*, ParticleCell*>>> _fprID_to_cellpairs;

    //! @brief cache state
    bool _processed;
};


#endif //MARDYN_INTERACTIONLOGPROCESSOR_H
