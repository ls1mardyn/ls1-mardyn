//
// Created by alex on 7/15/24.
//

#ifndef MARDYN_VCPADR_WRAPPER_H
#define MARDYN_VCPADR_WRAPPER_H


#include "VCPADR.h"

#include "particleContainer/adapter/CellProcessor.h"
#include "particleContainer/adapter/VectorizedCellProcessor.h"

#include <unordered_map>
#include <memory>

class VCPADR_Wrapper : public CellProcessor {
public:
    VCPADR_Wrapper(Domain & domain, double cutoffRadius, double LJcutoffRadius, const Resolution::Handler& resolutionHandler);

    /**
	 * Delegates to VCPADR::init
	 * */
    void init();

    /**
	 * \brief Delegates to both impls
	 */
    void initTraversal() override;

    /**
     * \brief Not needed
     */
    void preprocessCell(ParticleCell& /*cell*/) override {}

    /**
     * \brief Not needed
     */
    double processSingleMolecule(Molecule* /*m1*/, ParticleCell& /*cell2*/) override { return 0.0; }

    /**
     * \brief Calculate forces between pairs of Molecules in cell.
     */
    void processCell(ParticleCell& cell) override;

    /**
     * \brief Calculate forces between pairs of Molecules of two different cell.
     */
    void processCellPair(ParticleCell& cell1, ParticleCell& cell2, bool sumAll = false) override;

    /**
     * \brief Not needed
     */
    void postprocessCell(ParticleCell& /*cell*/) override {}

    /**
     * \brief Store macroscopic values in the Domain.
     */
    void endTraversal() override;

private:
    /**
     * @brief Checks which impl is required to call, and updates buffers if cell is unknown
     * @returns true if cell requires AdResS implementation, false else
     * */
    bool checkCell(ParticleCell& cell);

    /// ls1 implementation of the vectorized cell processor
    std::unique_ptr<VectorizedCellProcessor> _reference_processor;
    /// AdResS extended implementation of the vectorized cell processor
    std::unique_ptr<VCPADR> _adr_processor;
    /// mapping from cells to whether AdResS should be used (per thread)
    std::vector<std::unordered_map<ParticleCell*, bool>> _cell_maps;
    /// reference to active AdResS resolution handler
    const Resolution::Handler& _resolution_handler;
    /// flag if all cells have been processed
    bool _is_init;
};


#endif //MARDYN_VCPADR_WRAPPER_H
